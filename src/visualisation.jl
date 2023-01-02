"""
    setup_parametric_model(
        optimizer,
        func::Function,
        knots::Union{Nothing,Vector{Float64},Vector{Tuple{Float64,Symbol}}},
        lb::Real,
        ub::Real,
        δ::Real,
        type::Symbol,
    )

Internal function that creates a simple `JuMP.Model` to evaluate an interpolated function.
"""
function setup_parametric_model(
    optimizer,
    func::Function,
    knots::Union{Nothing,Vector{Float64},Vector{Tuple{Float64,Symbol}}},
    lb::Real,
    ub::Real,
    δ::Real,
    type::Symbol,
)
    model = JuMP.Model(optimizer)

    @variable(model, lb <= x <= ub)

    @constraint(model, fixed_x, x == 0)

    @objective(model, Max, approximate(x, func, δ = δ, type = type, knots = knots))
    return model
end

"""
    plot_approximation(
        optimizer,
        func::Function,
        lb::Real,
        ub::Real,
        δ::Real,
        type::Symbol;
        detail::Int = 3,
        knots::Union{Nothing,Vector{Float64},Vector{Tuple{Float64,Symbol}}} = nothing,
    )

Function that returns a set of points corresponding to an approximation of the function `func`.
"""
function plot_approximation(
    optimizer,
    func::Function,
    lb::Real,
    ub::Real,
    δ::Real,
    type::Symbol;
    detail::Int = 3,
    knots::Union{Nothing,Vector{Float64},Vector{Tuple{Float64,Symbol}}} = nothing,
)
    knots2 = nothing
    if typeof(knots) === Nothing
        knots2 = Float64[]
    elseif typeof(knots) == Vector{Float64}
        knots2 = copy(knots)
    else
        knots2 = [knots[i][1] for i in eachindex(knots)]
    end

    insert!(knots2, 1, lb)
    push!(knots2, ub)

    x = []
    step = δ / 2 / detail
    for i in 1:length(knots2)-1
        if knots2[i+1] > lb && knots2[i] < ub
            x = [x; collect(max(lb, knots2[i]):step:min(ub, knots2[i+1]))]
        end
    end

    y = Float64[]
    model = setup_parametric_model(optimizer, func, knots, lb, ub, δ, type)
    for i in x
        set_normalized_rhs(model[:fixed_x], i)
        map(Main.eval, [:(model = $model)])
        optimize!(model)
        obj_fn = objective_value(model)
        push!(y, obj_fn)
    end

    fx = Float64[]
    for i in eachindex(x)
        fxi = func(x[i])[1]
        if typeof(fxi) <: Real
            push!(fx, fxi)
        else
            push!(fx, fxi[1])
        end
    end
    return [x, fx, y]
end
