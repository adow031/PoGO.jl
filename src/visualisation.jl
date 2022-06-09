function setup_parametric_model(
    optimizer,
    func::Function,
    initial::Symbol,
    knots::Union{Nothing,Vector{Float64}},
    lb::Real,
    ub::Real,
    n::Int,
    type::Symbol,
)
    model = JuMP.Model(optimizer)

    @variable(model, lb <= x <= ub)

    @constraint(model, fixed_x, x == 0)

    @objective(
        model,
        Max,
        approximate(x, func, n, type = type, initial = initial, knots = knots)
    )
    return model
end

function plot_approximation(
    optimizer,
    func::Function,
    initial::Symbol,
    lb::Real,
    ub::Real,
    n::Int,
    type::Symbol;
    detail::Int = 3,
    knots::Union{Nothing,Vector{Float64}} = nothing,
)
    knots2 = copy(knots)
    insert!(knots2, 1, lb)
    push!(knots2, ub)

    x = []
    for i in 1:length(knots2)-1
        step = (knots2[i+1] - knots2[i]) / n / 2 / detail
        x = [x; collect(knots2[i]:step:knots2[i+1])]
    end

    y = Float64[]
    model = setup_parametric_model(optimizer, func, initial, knots, lb, ub, n, type)
    for i in x
        set_normalized_rhs(model[:fixed_x], i)
        optimize!(model)
        obj_fn = objective_value(model)
        push!(y, obj_fn)
    end

    Plots.plot(x, [func(x[i])[1] for i in 1:length(x)])
    return Plots.plot!(x, y)
end
