using PoGO, Test, JuMP, GLPK, ForwardDiff

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

function binary_bilinear2()
    model = Model(GLPK.Optimizer)

    @variable(model, x, Bin)
    @variable(model, z, Bin)
    @variable(model, 0 <= y <= 10)

    @constraint(model, x * y * z == 5)

    @objective(model, Min, 10 * x + y)

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
    @variable(model, 0 <= y <= 12, Int)
    @constraint(model, x + param * y == 10)

    @objective(model, Max, (bilinear(x, y) + bilinear(y, x)) / 2)

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
    model = JuMP.Model(GLPK.Optimizer)

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

    return [value(x), value(y)]
end

function interp(X::Real, Y::Real)
    model = JuMP.Model(GLPK.Optimizer)

    @variable(model, x)
    @variable(model, y)
    @variable(model, z)

    PoGO.interpolate(
        [x, y, z],
        [(0, 0, 3), (0, 2, 4), (1, 1, 5), (2, 0, 6), (2, 2, 7)],
        [0, 0, [0, 1], 1, 1],
    )

    @constraint(model, x == X)
    @constraint(model, y == Y)
    @objective(model, Min, z)

    optimize!(model)
    return value(z)
end

function interp_fn(X::Real, Y::Real, method::Symbol)
    model = JuMP.Model(GLPK.Optimizer)
    @variable(model, x)
    @variable(model, y)

    z = PoGO.interpolate_fn(
        (x, y) -> sqrt(x) * sin(y),
        [x, y],
        [1:0.1:3, 1:0.1:3],
        method = method,
    )
    @objective(model, Min, z)

    @constraint(model, x == X)
    @constraint(model, y == Y)

    optimize!(model)
    return value(z)
end

function interp_fn_pts(X, Y)
    points = [
        0.1345547259221871 0.6674206709045125
        0.5785505908199026 0.5799491307751639
        0.6937923537615268 0.00041410756224624645
        0.5546120304490526 0.444346423889622
        0.9736473421517288 0.36056976267756746
        0.3619967054023656 0.4045676376116868
        0.6859856822126517 0.5307393395868238
        0.693778523223462 0.12042346184441766
        0.5668391239777307 0.2965015701072703
        0.4549256807005051 0.28118851814924295
        0.521495486631126 0.38234635199801137
        0.40803859883694327 0.00038724533209688605
        0.7302285907371966 0.34618641979564346
        0.8789394232667177 0.5166243597119882
        0.19098457328317175 0.03307732223905446
        0.4488428264808636 0.2776225273989721
        0.6495747197271683 0.4629169708005898
        0.4486429598022251 0.8448075421920468
        0.4928407997818922 0.6830491066528377
        0.1520047909926191 0.9313189680471927
    ]

    model = JuMP.Model(GLPK.Optimizer)
    @variable(model, x)
    @variable(model, y)

    z = PoGO.interpolate_fn((x, y) -> [x * y, x^2 * y], [x, y], points)
    @objective(model, Min, z[1] + z[2])

    @constraint(model, x == X)
    @constraint(model, y == Y)

    optimize!(model)
    return objective_value(model)
end

function interp_pts(X, Y)
    points = [
        0.1345547259221871 0.6674206709045125 0.101888239521013
        0.5785505908199026 0.5799491307751639 0.5296509412286065
        0.6937923537615268 0.00041410756224624645 0.00048663443685222204
        0.5546120304490526 0.444346423889622 0.3831183903783676
        0.9736473421517288 0.36056976267756746 0.6928840128024092
        0.3619967054023656 0.4045676376116868 0.19946734842481678
        0.6859856822126517 0.5307393395868238 0.6138329724587234
        0.693778523223462 0.12042346184441766 0.14151047254757548
        0.5668391239777307 0.2965015701072703 0.2633365994113446
        0.4549256807005051 0.28118851814924295 0.1861139156095003
        0.521495486631126 0.38234635199801137 0.30337387119935894
        0.40803859883694327 0.00038724533209688605 0.00022248564718513688
        0.7302285907371966 0.34618641979564346 0.4373935197713632
        0.8789394232667177 0.5166243597119882 0.8531916632373886
        0.19098457328317175 0.03307732223905446 0.007523757148797848
        0.4488428264808636 0.2776225273989721 0.1805386817480868
        0.6495747197271683 0.4629169708005898 0.4960257351603552
        0.4486429598022251 0.8448075421920468 0.5490602452336192
        0.4928407997818922 0.6830491066528377 0.5025416684628153
        0.1520047909926191 0.9313189680471927 0.1630834949750745
    ]
    model = JuMP.Model(GLPK.Optimizer)
    @variable(model, x)
    @variable(model, y)

    z = PoGO.interpolate_points([x, y], points)
    @objective(model, Min, z)

    @constraint(model, x == X)
    @constraint(model, y == Y)

    optimize!(model)
    return objective_value(model)
end

function negative_power(rhs::Real)
    model = Model(GLPK.Optimizer)

    @variable(model, -6 <= x <= -1)
    @variable(model, -3 <= y <= 5, Int)

    z = power(x, y, 20)

    @constraint(model, z == rhs)

    @objective(model, Min, x)

    optimize!(model)

    return [value(model[:x]), value(model[:y]), value(z)]
end

function positive_power(rhs::Real)
    model = Model(GLPK.Optimizer)

    @variable(model, 1 <= x <= 6)
    @variable(model, -3 <= y <= 5, Int)

    z = PoGO.power(x, y, 20)

    @constraint(model, z == rhs)

    @objective(model, Min, x)

    optimize!(model)

    return [value(model[:x]), value(model[:y]), value(z)]
end

function domain_test(d::Real)
    model = JuMP.Model(GLPK.Optimizer)

    @variable(model, x)
    @variable(model, y >= 0)
    z = x ∈ [0, [2, 3]]

    @objective(model, Min, 10 * z[(2, 3)] + x + 5 * y)
    @constraint(model, x + y >= d)
    optimize!(model)

    return [value(x), value(y)]
end

function bounds_test()
    model = JuMP.Model()
    @variable(model, a[1:2], Bin)
    @variable(model, 0 <= b[1:2] <= 5, Int)

    return PoGO.find_bounds(a[1] - a[2] + b[1] - b[2])
end

function general_fn(; δ = 0.1, type = :interior)
    fn = function (x)
        if typeof(x) <: ForwardDiff.Dual && x == 0
            [-sqrt(abs(x) + 0.1), sqrt(abs(x) + 0.1)]
        elseif x < 1
            return sqrt(abs(x) + 0.1)
        elseif x == 1
            return [sqrt(abs(x) + 0.1), x - 1]
        elseif x < 2
            return x - 1
        elseif x == 2
            return [x - 1, 0.25 * x^2 + 1]
        else
            return 0.25 * x^2 + 1
        end
    end
    knots = [0.0, 1.0, 2.0]
    return x, fx, y = plot_approximation(
        GLPK.Optimizer,
        fn,
        -1.0,
        3.0,
        δ,
        type,
        detail = 1,
        knots = knots,
    )
end

@testset "PoGO.jl" begin
    result = nonlinear(10, :binary)
    @test sum(abs.(result - [1.6, 4.0, 0.37866])) ≈ 0.0 atol = 1e-4

    result = nonlinear(10, :bisection)
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

    @test binary_bilinear2() ≈ 15.0 atol = 1e-4

    @test integer_bilinear(1) ≈ 25.0 atol = 1e-4
    @test integer_bilinear(3) ≈ 8.0 atol = 1e-4

    @test bilinear_affexpr(80) ≈ 130.0 atol = 1e-4
    @test bilinear_affexpr(120) ≈ 140.0 atol = 1e-4

    @test bilinear_affexpr(120) ≈ 140.0 atol = 1e-4

    @test fixed_cost_investment(40) ≈ 270.0 atol = 1e-4
    @test fixed_cost_investment(60) ≈ 320.0 atol = 1e-4

    @test approximate_cubic(10, :upper) ≈ 2.1146 atol = 1e-4
    @test approximate_cubic(10, :lower) ≈ 2.088 atol = 1e-4
    @test approximate_cubic(10, :combined) ≈ 2.1003 atol = 1e-4

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

    @test interp(1, 1) ≈ 5.0 atol = 1e-4
    @test interp(1.5, 1.5) ≈ 6.0 atol = 1e-4
    @test interp(1.75, 1.75) ≈ 6.5 atol = 1e-4
    @test interp(2, 2) ≈ 7.0 atol = 1e-4
    @test interp(1.5, 0.5) ≈ 5.5 atol = 1e-4

    @test interp_fn(1.25, 2.25, :binary) ≈ 0.868 atol = 1e-3
    @test interp_fn(1, 2.25, :binary) ≈ 0.777 atol = 1e-3
    @test interp_fn(2.42, 1.21, :binary) ≈ 1.455 atol = 1e-3

    @test interp_fn(1, 2.25, :bisection) ≈ 0.777 atol = 1e-3
    @test interp_fn(2.42, 1.21, :bisection) ≈ 1.455 atol = 1e-3

    @test interp_fn_pts(0.5, 0.5) ≈ 0.394 atol = 1e-3
    @test interp_fn_pts(0.4, 0.4) ≈ 0.225 atol = 1e-3

    @test interp_pts(0.5, 0.5) ≈ 0.394 atol = 1e-3
    @test interp_pts(0.4, 0.4) ≈ 0.225 atol = 1e-3

    result = negative_power(-8)
    @test sum(abs.(result - [-1.9628, 3.0, -8.0])) ≈ 0.0 atol = 1e-2
    result = negative_power(-0.25)
    @test sum(abs.(result - [-4.2622, -1.0, -0.25])) ≈ 0.0 atol = 1e-2
    result = negative_power(16)
    @test sum(abs.(result - [-3.8890, 2.0, 16.0])) ≈ 0.0 atol = 1e-2
    result = negative_power(4)
    @test sum(abs.(result - [-1.9425, 2.0, 4.0])) ≈ 0.0 atol = 1e-2

    result = positive_power(2)
    @test sum(abs.(result - [1.1413, 5.0, 2.0])) ≈ 0.0 atol = 1e-2

    result = domain_test(2)
    @test sum(abs.(result - [0.0, 2.0])) ≈ 0.0 atol = 1e-2
    result = domain_test(3)
    @test sum(abs.(result - [3.0, 0.0])) ≈ 0.0 atol = 1e-2
    result = domain_test(4)
    @test sum(abs.(result - [3.0, 1.0])) ≈ 0.0 atol = 1e-2

    @test PoGO.process_knots([(1.0, :convex), (4.0, :concave)], 0, 5) ==
          ([0.0, 1.0, 4.0, 5.0], [:convex, :convex, :concave, :concave])
    @test PoGO.process_knots([(1.0, :convex), (4.0, :concave)], 2, 5) ==
          ([2.0, 4.0, 5.0], [:convex, :concave, :concave])
    @test PoGO.process_knots([(1.0, :convex), (4.0, :concave)], 0, 3) ==
          ([0.0, 1.0, 3.0], [:convex, :convex, :convex])
    @test PoGO.process_knots([(1.0, :convex), (4.0, :concave)], 2, 3) ==
          ([2.0, 3.0], [:convex, :convex])

    @test bounds_test() == (-6.0, 6.0)

    result = general_fn()
    @test sum(
        sum.(
            abs.(
                vcat(result...) - vcat(
                    [
                        [
                            -1.0,
                            -0.95,
                            -0.9,
                            -0.85,
                            -0.8,
                            -0.75,
                            -0.7,
                            -0.65,
                            -0.6,
                            -0.55,
                            -0.5,
                            -0.45,
                            -0.4,
                            -0.35,
                            -0.3,
                            -0.25,
                            -0.2,
                            -0.15,
                            -0.1,
                            -0.05,
                            0.0,
                            0.0,
                            0.05,
                            0.1,
                            0.15,
                            0.2,
                            0.25,
                            0.3,
                            0.35,
                            0.4,
                            0.45,
                            0.5,
                            0.55,
                            0.6,
                            0.65,
                            0.7,
                            0.75,
                            0.8,
                            0.85,
                            0.9,
                            0.95,
                            1.0,
                            1.0,
                            1.05,
                            1.1,
                            1.15,
                            1.2,
                            1.25,
                            1.3,
                            1.35,
                            1.4,
                            1.45,
                            1.5,
                            1.55,
                            1.6,
                            1.65,
                            1.7,
                            1.75,
                            1.8,
                            1.85,
                            1.9,
                            1.95,
                            2.0,
                            2.0,
                            2.05,
                            2.1,
                            2.15,
                            2.2,
                            2.25,
                            2.3,
                            2.35,
                            2.4,
                            2.45,
                            2.5,
                            2.55,
                            2.6,
                            2.65,
                            2.7,
                            2.75,
                            2.8,
                            2.85,
                            2.9,
                            2.95,
                            3.0,
                        ],
                        [
                            1.0488088481701516,
                            1.02469507659596,
                            1.0,
                            0.9746794344808963,
                            0.9486832980505138,
                            0.9219544457292888,
                            0.8944271909999159,
                            0.8660254037844386,
                            0.8366600265340756,
                            0.806225774829855,
                            0.7745966692414834,
                            0.7416198487095663,
                            0.7071067811865476,
                            0.6708203932499369,
                            0.6324555320336759,
                            0.5916079783099616,
                            0.5477225575051662,
                            0.5,
                            0.4472135954999579,
                            0.3872983346207417,
                            0.31622776601683794,
                            0.31622776601683794,
                            0.3872983346207417,
                            0.4472135954999579,
                            0.5,
                            0.5477225575051662,
                            0.5916079783099616,
                            0.6324555320336759,
                            0.6708203932499369,
                            0.7071067811865476,
                            0.7416198487095663,
                            0.7745966692414834,
                            0.806225774829855,
                            0.8366600265340756,
                            0.8660254037844386,
                            0.8944271909999159,
                            0.9219544457292888,
                            0.9486832980505138,
                            0.9746794344808963,
                            1.0,
                            1.02469507659596,
                            1.0488088481701516,
                            1.0488088481701516,
                            0.050000000000000044,
                            0.10000000000000009,
                            0.1499999999999999,
                            0.19999999999999996,
                            0.25,
                            0.30000000000000004,
                            0.3500000000000001,
                            0.3999999999999999,
                            0.44999999999999996,
                            0.5,
                            0.55,
                            0.6000000000000001,
                            0.6499999999999999,
                            0.7,
                            0.75,
                            0.8,
                            0.8500000000000001,
                            0.8999999999999999,
                            0.95,
                            1.0,
                            1.0,
                            2.050625,
                            2.1025,
                            2.1556249999999997,
                            2.21,
                            2.265625,
                            2.3225,
                            2.380625,
                            2.44,
                            2.5006250000000003,
                            2.5625,
                            2.625625,
                            2.6900000000000004,
                            2.755625,
                            2.8225000000000002,
                            2.890625,
                            2.96,
                            3.030625,
                            3.1025,
                            3.175625,
                            3.25,
                        ],
                        [
                            1.0488088481701514,
                            1.0244044240850756,
                            1.0,
                            0.9743416490252568,
                            0.9486832980505138,
                            0.9215552445252149,
                            0.8944271909999157,
                            0.8655436087669957,
                            0.8366600265340757,
                            0.8056283478877795,
                            0.7745966692414834,
                            0.7408517252140154,
                            0.7071067811865476,
                            0.6697811566101117,
                            0.6324555320336759,
                            0.590089044769421,
                            0.5477225575051661,
                            0.49746807650256203,
                            0.447213595499958,
                            0.38172068075839793,
                            0.31622776601683794,
                            0.31622776601683794,
                            0.38172068075839793,
                            0.4472135954999579,
                            0.497468076502562,
                            0.5477225575051662,
                            0.590089044769421,
                            0.6324555320336759,
                            0.6697811566101117,
                            0.7071067811865476,
                            0.7408517252140154,
                            0.7745966692414834,
                            0.8056283478877796,
                            0.8366600265340756,
                            0.8655436087669957,
                            0.8944271909999159,
                            0.9215552445252149,
                            0.9486832980505138,
                            0.9743416490252568,
                            1.0,
                            1.0244044240850758,
                            1.0488088481701516,
                            1.0488088481701516,
                            0.050000000000000044,
                            0.10000000000000009,
                            0.1499999999999999,
                            0.19999999999999996,
                            0.25,
                            0.30000000000000004,
                            0.3500000000000001,
                            0.3999999999999999,
                            0.44999999999999996,
                            0.5,
                            0.55,
                            0.6000000000000001,
                            0.6499999999999999,
                            0.7,
                            0.75,
                            0.8,
                            0.8500000000000001,
                            0.8999999999999999,
                            0.95,
                            2.0,
                            2.0,
                            2.0512499999999996,
                            2.1025,
                            2.1562499999999996,
                            2.21,
                            2.2662499999999994,
                            2.3225,
                            2.3812499999999996,
                            2.4399999999999995,
                            2.5012499999999998,
                            2.5624999999999996,
                            2.6262500000000006,
                            2.690000000000001,
                            2.75625,
                            2.8225000000000002,
                            2.8912499999999994,
                            2.96,
                            3.03125,
                            3.1024999999999996,
                            3.1762500000000005,
                            3.25,
                        ],
                    ]...,
                ),
            ),
        ),
    ) ≈ 0.0 atol = 1e-2
end
