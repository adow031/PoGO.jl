module PoGO

using JuMP, Plots, ForwardDiff, Delaunay

import Base: ==

(==)(x::U, y::T) where {U<:ForwardDiff.Dual,T<:Int} = x.value == y

include("utilities.jl")
include("approximation.jl")
include("bilinear.jl")
include("interpolation.jl")
include("visualisation.jl")

export bilinear,
    approximate, power, plot_approximation, interpolate, interpolate_fn, interpolate_points
end
