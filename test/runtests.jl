using PoGO, Test, JuMP, GLPK

function nonlinear(n::Int, method::Symbol)
    model = JuMP.Model(GLPK.Optimizer)
    @variable(model, 0 <= x <= 4)
    @variable(model, 0 <= y <= 4)
    @constraint(
        model,
        approximate(x, a -> a^2, n, method = method) +
        approximate(y, a -> a^2, n, method = method) >= 9
    )

    @objective(
        model,
        Max,
        bilinear(
            approximate(x, a -> exp(-a / 10) * sin(a), n, method = method),
            approximate(y, a -> exp(-a / 5) * sin(2 * a), n, method = method),
            n,
        )
    )

    optimize!(model)

    value(model[:y])

    return [value(model[:x]), value(model[:y]), objective_value(model)]
end

function binary_bilinear(profit_xy, cost_x, cost_y)
    model = JuMP.Model(GLPK.Optimizer)

    @variable(model, x, Bin)
    @variable(model, y, Bin)

    @objective(model, Max, profit_xy * bilinear(x, y) - cost_x * x - cost_y * y)

    optimize!(model)
    return objective_value(model)
end

function fixed_cost_investment(cost)
    model = JuMP.Model(GLPK.Optimizer)

    @variable(model, 0 <= z <= 3, Int)
    @variable(model, 0 <= capacity <= 40)
    @variable(model, x >= 0)
    @variable(model, y >= 0)

    @constraint(model, x <= bilinear(z, capacity))
    @constraint(model, x + y == 50)

    @objective(model, Min, x + 10 * y + 6 * capacity + cost * z)

    optimize!(model)
    return value(z) | return objective_value(model)
end

function integer_bilinear(param::Real)
    model = JuMP.Model(GLPK.Optimizer)

    @variable(model, 0 <= x <= 10, Int)
    @variable(model, 0 <= y <= 10, Int)
    @constraint(model, x + param * y == 10)

    @objective(model, Max, bilinear(x, y))

    optimize!(model)
    return objective_value(model)
end

function approximate_cubic(n::Int, type::Symbol)
    model = JuMP.Model(GLPK.Optimizer)

    @variable(model, -1 <= x <= 2.5)

    @objective(
        model,
        Max,
        approximate(
            x,
            a -> a^3 - 3a^2 + a + 2,
            n,
            type = type,
            knots = [(-1.0, :concave), (1.0, :convex)],
        )
    )
    optimize!(model)
    return objective_value(model)
end

function bilinear_affexpr(cost)
    model = JuMP.Model(GLPK.Optimizer)

    @variable(model, z[1:2], Bin)
    @variable(model, 0 <= capacity <= 40)
    @variable(model, x >= 0)
    @variable(model, y >= 0)

    @constraint(model, x <= bilinear(sum(z), capacity))
    @constraint(model, x + y == 50)

    @objective(model, Min, x + 10 * y + cost * bilinear(z[1], z[2]))

    optimize!(model)

    return objective_value(model)
end

function approximate_cubic_discrete(obj, lb, ub)
    model = JuMP.Model(GLPK.Optimizer)

    @variable(model, lb <= x <= ub, Int)

    if obj == :max
        @objective(model, Max, approximate(x, a -> a^3 - 3a^2 + a + 2))
    else
        @objective(model, Min, approximate(x, a -> a^3 - 3a^2 + a + 2))
    end
    optimize!(model)

    return objective_value(model)
end

function in_sets(cx::Real, cy::Real)
    sets = [
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(2.0, 2.0), (3.0, 2.0), (3.0, 3.0), (2.0, 3.0)],
    ]

    model = JuMP.Model(GLPK.Optimizer)

    @variable(model, x)
    @variable(model, y)
    @variable(model, d)

    z = xy_in_sets(x, y, sets, update_bounds = true)

    @constraint(
        model,
        d >= approximate(x - cx, a -> a^2, 21) + approximate(y - cy, a -> a^2, 21)
    )

    @objective(model, Min, d)

    optimize!(model)

    return [value(x), value(y)]
end

result = nonlinear(10, :binary)
approximate_cubic(10, :upper)
@testset "PoGO.jl" begin
    result = nonlinear(10, :binary)
    @test sum(abs.(result - [1.6, 4.0, 0.37866])) ≈ 0.0 atol = 1e-4

    result = nonlinear(20, :echelon)
    @test sum(abs.(result - [1.4, 3.8, 0.3878])) ≈ 0.0 atol = 1e-4

    result = nonlinear(50, :echelon)
    @test sum(abs.(result - [1.44, 3.84, 0.3923])) ≈ 0.0 atol = 1e-4

    result = nonlinear(50, :echelon)
    @test sum(abs.(result - [1.44, 3.84, 0.3923])) ≈ 0.0 atol = 1e-4

    @test binary_bilinear(10, 4, 5) ≈ 1.0 atol = 1e-4
    @test binary_bilinear(10, 4, 1) ≈ 5.0 atol = 1e-4
    @test binary_bilinear(10, 6, 5) ≈ 0.0 atol = 1e-4

    @test integer_bilinear(1) ≈ 25.0 atol = 1e-4
    @test integer_bilinear(3) ≈ 8.0 atol = 1e-4

    @test bilinear_affexpr(80) ≈ 130.0 atol = 1e-4
    @test bilinear_affexpr(120) ≈ 140.0 atol = 1e-4

    @test fixed_cost_investment(40) ≈ 270.0 atol = 1e-4
    @test fixed_cost_investment(60) ≈ 320.0 atol = 1e-4

    @test approximate_cubic(10, :upper) ≈ 2.1146 atol = 1e-4
    @test approximate_cubic(10, :lower) ≈ 2.088 atol = 1e-4

    @test approximate_cubic_discrete(:max, -1, 4.5) ≈ 22 atol = 1e-4
    @test approximate_cubic_discrete(:max, -1, 2.5) ≈ 2 atol = 1e-4
    @test approximate_cubic_discrete(:min, -1, 2) ≈ -3 atol = 1e-4
    @test approximate_cubic_discrete(:min, -1.5, 2) ≈ -3 atol = 1e-4

    @test isapprox(in_sets(-1, -1), [0.0, 0.0]; atol = 1e-4)
    @test isapprox(in_sets(-1, 1), [0.0, 1.0]; atol = 1e-4)
    @test isapprox(in_sets(1, -1), [1.0, 0.0]; atol = 1e-4)
    @test isapprox(in_sets(1.2, 1.2), [1.0, 1.0]; atol = 1e-4)
    @test isapprox(in_sets(1.8, 1.8), [2.0, 2.0]; atol = 1e-4)
    @test isapprox(in_sets(1, 4), [2.0, 3.0]; atol = 1e-4)
    @test isapprox(in_sets(4, 1), [3.0, 2.0]; atol = 1e-4)
    @test isapprox(in_sets(4, 4), [3.0, 3.0]; atol = 1e-4)
end
