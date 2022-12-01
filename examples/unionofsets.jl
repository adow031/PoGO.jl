using JuMP, PoGO, Plots

# Uncomment appropriate solver

using Gurobi
optimizer =
    optimizer_with_attributes(() -> Gurobi.Optimizer(), "OutputFlag" => 1, "MIPGap" => 0.0)

function in_sets(cx::Real, cy::Real)
    points = [
        (0.0, 0.0),
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (2.0, 2.0),
        (3.0, 2.0),
        (3.0, 3.0),
        (2.0, 3.0),
    ]

    sets = [1, 1, 1, 1, 2, 2, 2, 2]

    model = JuMP.Model(optimizer)
    @variable(model, x)
    @variable(model, y)

    interpolate([x, y], points, sets, update_bounds = true)

    @objective(model, Min, cx * x + cy * y)

    optimize!(model)

    return value(x), value(y)
end

in_sets(1, 1) == (0.0, 0.0)
in_sets(-2, 1) == (3.0, 2.0)
in_sets(-1, -1) == (3.0, 3.0)
in_sets(1, -2) == (2.0, 3.0)
in_sets(2, -1) == (0.0, 1.0)
in_sets(-1, 2) == (1.0, 0.0)

function in_sets2(cx::Real, cy::Real)
    points = [
        (0.0, 0.0),
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (2.0, 2.0),
        (3.0, 2.0),
        (3.0, 3.0),
        (2.0, 3.0),
    ]

    sets = [1, 1, 1, 1, 2, 2, 2, 2]

    model = JuMP.Model(optimizer)

    @variable(model, x)
    @variable(model, y)
    @variable(model, d)

    interpolate([x, y], points, sets, update_bounds = true)

    @constraint(
        model,
        d >= approximate(x - cx, a -> a^2, 21) + approximate(y - cy, a -> a^2, 21)
    )

    @objective(model, Min, d)

    optimize!(model)

    plot([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)], legend = false)
    plot!([(2.0, 2.0), (3.0, 2.0), (3.0, 3.0), (2.0, 3.0), (2.0, 2.0)])
    plot!((cx, cy), seriestype = :scatter)
    return plot!((value(x), value(y)), seriestype = :scatter)
end

in_sets2(-1, -1)
in_sets2(-1, 1)
in_sets2(1, -1)
in_sets2(1.2, 1.2)
in_sets2(1.8, 1.8)
in_sets2(1, 4)
in_sets2(4, 1)
in_sets2(4, 4)

function dynamic_set()
    model = JuMP.Model(optimizer)
    @variable(model, x)
    @variable(model, y)
    @variable(model, 0 <= z <= 1)

    points = [(0.0, z), (0.0, 1.0), (0.5, 1.0), (z, 0.4)]
    sets = [1, 1, 1, 1]
    interpolate([x, y], points, sets, 10, update_bounds = true)

    @objective(model, Min, -2 * x - 1 * y - z)
    optimize!(model)

    plot(
        [(0.0, value(z)), (0.0, 1.0), (0.5, 1.0), (value(z), 0.4), (0.0, value(z))],
        legend = false,
    )
    return plot!((value(x), value(y)), seriestype = :scatter)
end

dynamic_set()
