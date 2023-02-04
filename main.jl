include("hodgkin_huxley.jl")
using .HodgkinHuxley
using OrdinaryDiffEq

hodgkin_huxley(
    method = Rodas4(),
    time_span = (0.0, 500.0),
    current_injected₀ = 0.06,
    time₁ = 20.0,
    duration₁ = 500.0,
    current_injected₁ = 0.0,
    time₂ = 20.0,
    duration₂ = 500.0,
    current_injected₂ = 0.0,
    relative_tolerance = 1e-11,
    absolute_tolerance = 1e-11,
    show_plot = true,
    save_figure = false # or filename of output image
    )