module PoGO

using JuMP, Plots

include("utilities.jl")
include("approximation.jl")
include("bilinear.jl")
include("visualisation.jl")

export bilinear, approximate, exponential, plot_approximation
end
