using OrdinaryDiffEq
using Plots

Ena = 50
Ek = -77
El = -54.387
Gna = 1.20
Gk = 0.36
Gl = 0.003
Cm = 0.01
 
 function hodgkin_huxley(plotting::Bool, method)
    time_span = (0.0, 500.0)

    current_injected₀ = 0.06

    time₁ = 20
    duration₁ = 500
    current_injected₁ = 0.0

    time₂ = 20
    duration₂ = 500
    current_injected₂ = 0.0

    voltage₀ = -65
    n_gate₀ = 0.3
    m_gate₀ = 0.1
    h_gate₀ = 0.6

    relative_tolerance = 1e-11
    absolute_tolerance = 1e-11
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

    current_Na = Gna .* m_gate .^3 .* h_gate .* (Ena .- voltage);
    current_K = Gk .* n_gate .^4 .* (Ek .- voltage);
    current_l = Gl .* (El .- voltage);

    print(size(current_Na))
    print(size(current_K))
    print(size(current_l))

    plot_voltage_v_time = plot(solution, plotdensity=10000, idxs=[1])
    plot_gates_v_time = plot(solution, plotdensity=10000, idxs=[2, 3, 4])
    plot_current_v_time = plot(time, plotdensity=10000, [current_Na, current_K, current_l])

    plot(plot_voltage_v_time, plot_gates_v_time, plot_current_v_time, layout=(3, 1))
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

function hodgkin_huxley()
    hodgkin_huxley(true, QNDF())
end

hodgkin_huxley()