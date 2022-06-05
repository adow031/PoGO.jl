module Panco

using JuMP

function get_model(x::Union{VariableRef,AffExpr})
    if typeof(x) == VariableRef
        return x.model
    else
        return collect(keys(x.terms))[1].model
    end
end

function approximate(
    x::Union{VariableRef,AffExpr},
    func::Function,
    n::Int;
    convex::Bool = false,
)
    model = get_model(x)

    fx = @variable(model)

    lx, ux = find_bounds(x)

    xi = []
    fxi = []
    for i in 0:n
        push!(xi, lx + i * (ux - lx) / n)
        push!(fxi, func(xi[end]))
    end

    set_lower_bound(fx, minimum(fxi))
    set_upper_bound(fx, maximum(fxi))

    α = Dict()

    if !convex
        ϵα = Dict()
    end

    for i in 0:n
        if i != n && !convex
            ϵα[i] = @variable(model, binary = true)
        end
        α[i] = @variable(model)
        set_lower_bound(α[i], 0)
        set_upper_bound(α[i], 1)
    end

    @constraint(model, sum(α[i] for i in 0:n) == 1)
    @constraint(model, sum(α[i] * xi[i+1] for i in 0:n) == x)
    @constraint(model, sum(α[i] * fxi[i+1] for i in 0:n) == fx)

    if !convex
        for i in 2:n-1
            @constraint(model, α[i] <= ϵα[i] - ϵα[i-2])
        end

        for i in 0:n-2
            @constraint(model, ϵα[i] <= ϵα[i+1])
        end

        @constraint(model, α[0] <= ϵα[0])
        @constraint(model, α[1] <= ϵα[1])
        @constraint(model, α[n] <= 1 - ϵα[n-2])
    end

    return fx
end

function find_bounds(x::Union{VariableRef,AffExpr})
    if typeof(x) == VariableRef
        lx = lower_bound(x)
        ux = upper_bound(x)
    else
        lx = x.constant
        ux = x.constant
        for (var, coeff) in x.terms
            if coeff > 0
                lx += coeff * lower_bound(var)
                ux += coeff * upper_bound(var)
            elseif coeff < 0
                lx += coeff * upper_bound(var)
                ux += coeff * lower_bound(var)
            end
        end
    end
    return lx, ux
end

function bilinear(x::Union{VariableRef,AffExpr}, y::Union{VariableRef,AffExpr}, n::Int)
    model = get_model(x)

    a² = @variable(model)
    b² = @variable(model)
    xy = @variable(model)

    lx, ux = find_bounds(x)
    ly, uy = find_bounds(y)

    set_lower_bound(xy, min(ux * ly, lx * uy))
    set_upper_bound(xy, max(ux * uy, lx * ly))

    a² = approximate(
        ((x - (lx + ux) / 2) / (ux - lx) + (y - (ly + uy) / 2) / (uy - ly)) / 2,
        x -> x^2,
        2 * n,
    )
    b² = approximate(
        ((x - (lx + ux) / 2) / (ux - lx) - (y - (ly + uy) / 2) / (uy - ly)) / 2,
        x -> x^2,
        2 * n,
    )
    @constraint(
        model,
        xy - x * (ly + uy) / 2 - y * (lx + ux) / 2 + (lx + ux) * (ly + uy) / 4 ==
        (a² - b²) * (ux - lx) * (uy - ly)
    )

    return xy
end

export bilinear, approximate
end
