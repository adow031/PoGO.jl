module PoGO

using JuMP, Plots

include("utilities.jl")
include("approximation.jl")
include("bilinear.jl")
include("visualisation.jl")

export bilinear, approximate, power, plot_approximation
end
