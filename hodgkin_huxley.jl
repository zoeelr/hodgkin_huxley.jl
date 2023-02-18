module HodgkinHuxley

    export solve_hodgkin_huxley

    using OrdinaryDiffEq
    using Plots

    Ena = 50
    Ek = -77
    El = -54.387
    Gna = 1.20
    Gk = 0.36
    Gl = 0.003
    Cm = 0.01
 
    function solve_hodgkin_huxley(method, time_span, current_injected₀, time₁, duration₁, current_injected₁, time₂, duration₂, current_injected₂, relative_tolerance, absolute_tolerance, show_plot::Bool, save_figure::Union{Bool, AbstractString} = false)

        voltage₀ = -65
        n_gate₀ = 0.3
        m_gate₀ = 0.1
        h_gate₀ = 0.6

        maximum_time_step = max(duration₁, duration₂)/2

        u₀ = [voltage₀, n_gate₀, m_gate₀, h_gate₀]

        p = (current_injected₀, time₁, duration₁, current_injected₁, time₂, duration₂, current_injected₂)

        ode_problem = ODEProblem(hodgkin_huxley_rhs, u₀, time_span, p)

        solution = solve(ode_problem, method, reltol = relative_tolerance, abstol = absolute_tolerance, dtmax = maximum_time_step)

        time = solution.t
        voltage = solution[1, :]
        n_gate = solution[2, :]
        m_gate = solution[3, :]
        h_gate = solution[4, :]

        spike_times, periods = calculate_periods(time, voltage)

        current_Na = Gna .* m_gate .^3 .* h_gate .* (Ena .- voltage);
        current_K = Gk .* n_gate .^4 .* (Ek .- voltage);
        current_l = Gl .* (El .- voltage);


        firing_rate = calculate_firing_rate(time, spike_times)

        if !(isnothing(firing_rate))
            firing_rate_to_plot, time_to_plot = prepare_firing_rate_to_plot(firing_rate, time, spike_times)
        end

        if show_plot || (save_figure != false)
            plot_voltage_v_time = plot(solution, plotdensity=10000, idxs=[1], title="Voltage", label="voltage", xlabel="time (ms)", ylabel="(mV)"; dpi=600)
            plot_gates_v_time = plot(solution, plotdensity=10000, idxs=[2, 3, 4], title="Current gating variables", label=["n" "m" "h"], xlabel="time (ms)"; dpi=600)
            plot_current_v_time = plot(time, plotdensity=10000, [current_Na, current_K, current_l], title="Currents", label=["Na" "K" "leak"], xlabel="time (ms)", ylabel="(mA)"; dpi=600)
            if !(isnothing(firing_rate))
                plot_spike_periods = plot(time_to_plot, plotdensity=10000, firing_rate_to_plot, title="Spike firing rate", label="firing rate", xlabel="time (ms)", dpi=600)
                full_plot = plot(plot_voltage_v_time, plot_gates_v_time, plot_current_v_time, plot_spike_periods, layout=(2, 2); dpi = 600)
            else
                full_plot = plot(plot_voltage_v_time, plot_gates_v_time, plot_current_v_time, layout=(3, 1); dpi = 600)
            end

            if (save_figure != false)
                savefig(plot_voltage_v_time, save_figure*"_volt")
                savefig(plot_gates_v_time, save_figure*"_vars")
                savefig(plot_current_v_time, save_figure*"_currs")
                if !(isnothing(firing_rate))
                savefig(plot_spike_periods, save_figure*"_frrt")
                end
                savefig(full_plot, save_figure*"_full")
            end
            if show_plot
                display(full_plot)
            end
        end

        return (solution, spike_times, periods)
        
    end

    function calculate_periods(time, voltage)
        mask_voltage_sign_changes = (voltage[1:end-1] .* voltage[2:end]) .< 0
    
        mask_voltage_negative = voltage[1:end-1] .< 0
    
        mask_spikes = (mask_voltage_sign_changes .== 1) .& (mask_voltage_negative .== 1)
    
        spike_times = (time[1:end-1])[mask_spikes]
    
        if length(mask_spikes) < 2
            @warn "Less than 2 spikes, no periods calculated"
            return spike_times, []
        end
    
        periods = spike_times[2:end] .- spike_times[1:end-1]
    
        return spike_times, periods
    end


    function calculate_periods(time, voltage)
        mask_voltage_sign_changes = (voltage[1:end-1] .* voltage[2:end]) .< 0
    
        mask_voltage_negative = voltage[1:end-1] .< 0
    
        mask_spikes = (mask_voltage_sign_changes .== 1) .& (mask_voltage_negative .== 1)
    
        spike_times = (time[1:end-1])[mask_spikes]
    
        if length(mask_spikes) < 2
            @warn "Less than 2 spikes, no periods calculated"
            return spike_times, []
        end
    
        periods = spike_times[2:end] .- spike_times[1:end-1]
    
        return spike_times, periods
    end

    function calculate_firing_rate(time, spike_times)
        if length(spike_times) < 2
            @warn "Less than 2 spikes, no periods calculated"
            return nothing
        end
        modified_spike_times = copy(spike_times)

        append!(modified_spike_times, time[end] + 1)
    
        next_spike = zeros((length(spike_times), length(time)))
        most_recent_spike = zeros((length(spike_times), length(time)))
        for spike_index in eachindex(spike_times)
            most_recent_spike[spike_index, :] = modified_spike_times[spike_index] .* ( (time .>= modified_spike_times[spike_index]) .&& (time .< modified_spike_times[spike_index + 1]) )
        end
        extrapolated_next_spike = 2 * spike_times[end] - spike_times[end - 1]
        modified_spike_times[end] = extrapolated_next_spike
        for spike_index in eachindex(spike_times)
            next_spike[spike_index, :] = modified_spike_times[spike_index + 1] .* ( (time .>= modified_spike_times[spike_index]) .&& (time .< modified_spike_times[spike_index + 1]) )
        end
    
        most_recent_spike[end] = most_recent_spike[end - 1]
        next_spike[end] = extrapolated_next_spike
    
        most_recent_spike = sum(most_recent_spike, dims=1)
        next_spike = sum(next_spike, dims=1)
    
        for index in eachindex(next_spike)
            if next_spike[index] != 0
                break
            end
            next_spike[index] = modified_spike_times[1]
        end
    
        interspike_interval = next_spike .- most_recent_spike
        firing_rate = 1 ./ interspike_interval
    
        firing_rate = firing_rate[1:end-1]
    
        return firing_rate
    end

    function prepare_firing_rate_to_plot(firing_rate, time, spike_times)
        time_to_plot = copy(time)[1:end - 1]
        firing_rate_to_plot = copy(firing_rate)
        firing_rate_to_plot = firing_rate_to_plot[time_to_plot .> spike_times[1]]
        time_to_plot = time_to_plot[time_to_plot .> spike_times[1]]
        return firing_rate_to_plot, time_to_plot
    end

    function hodgkin_huxley_rhs(u, p, time)
        (current_injected₀, time₁, duration₁, current_injected₁, time₂, duration₂, current_injected₂) = p

        (voltage, n_gate, m_gate, h_gate) = u

        current_injected = current_injected₀ + current_injected₁ * (time >= time₁) * (time <= time₁ + duration₁) + current_injected₂ * (time >= time₂) * (time <= time₂ + duration₂)

        voltage_threshold = -40

        αm = 0.1 * (voltage - voltage_threshold)/(1 - exp(-0.1 * (voltage - voltage_threshold)))
        βm = 4 * exp(-0.0556 * (voltage + 65))
        m_prime = αm * (1 - m_gate) - βm * m_gate

        αh = 0.07 * exp(-0.05 * (voltage + 65))
        βh = 1 / (1 + exp(-0.1 * (voltage + 35)))
        h∞ = αh / (αh + βh)
        τh = 1 / (αh + βh)
        h_prime = (h∞ - h_gate) / τh

        voltage_threshold = -55
        αn = 0.01 * (voltage - voltage_threshold) / (1 - exp(-0.1 * (voltage - voltage_threshold)))
        βn = 0.125 * exp(-0.0125 * (voltage + 65))
        n∞ = αn / (αn + βn)
        τn = 1 / (αn + βn)
        n_prime = (n∞ - n_gate) / τn

        voltage_prime = 1 / Cm * (Gna * m_gate^3 * h_gate * (Ena - voltage) + Gk * n_gate^4 * (Ek - voltage) + Gl *(El - voltage) + current_injected)

        u_prime = [voltage_prime, n_prime, m_prime, h_prime]

        return u_prime
    end

    function solve_hodgkin_huxley(;method, time_span, current_injected₀, time₁, duration₁, current_injected₁, time₂, duration₂, current_injected₂, relative_tolerance, absolute_tolerance, show_plot::Bool, save_figure::Union{Bool, AbstractString})
        return solve_hodgkin_huxley(method, time_span, current_injected₀, time₁, duration₁, current_injected₁, time₂, duration₂, current_injected₂, relative_tolerance, absolute_tolerance, show_plot::Bool, save_figure::Union{Bool, AbstractString})
    end

end