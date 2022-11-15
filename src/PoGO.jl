module PoGO

using JuMP, Plots, ForwardDiff

import Base: ==

(==)(x::U, y::T) where {U<:ForwardDiff.Dual,T<:Int} = x.value == y

include("utilities.jl")
include("approximation.jl")
include("bilinear.jl")
include("2dsets.jl")
include("visualisation.jl")

export bilinear, approximate, power, plot_approximation, xy_in_sets
end
