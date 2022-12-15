module PoGO

using JuMP, ForwardDiff, Delaunay

import Base: ==, *

(==)(x::U, y::T) where {U<:ForwardDiff.Dual,T<:Int} = x.value == y

include("utilities.jl")
include("approximation.jl")
include("bilinear.jl")
include("interpolation.jl")
include("visualisation.jl")

function (*)(x::U, y::T) where {U<:Union{VariableRef,AffExpr},T<:Union{VariableRef,AffExpr}}
    return bilinear(x, y)
end

export bilinear,
    approximate, power, plot_approximation, interpolate, interpolate_fn, interpolate_points
end
