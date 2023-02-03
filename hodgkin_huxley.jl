using OrdinaryDiffEq

hodgkin_huxley = function(plotting::Bool)
    time_span = [0, 500]
    tolerance = 1e-11

    method = QNDF()

    current_injected₀ = 0.06

    time₁ = 20
    duration₁ = 500
    current_injected₁ = 0.0

    volatage₀ = -65
    variable_n₀ = 0.3
    variable_m₀ = 0.1
    variable_h₀ = 0.6

    u₀ = [volatage₀, variable_n₀, variable_m₀, variable_h₀]

    ode_problem = ODEProblem(X, u₀, time_span)

    (time, u) = solve(ode_problem, method)
end

hodgkin_huxley_rhs = function ()
    
end