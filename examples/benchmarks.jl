using JuMP, PoGO, Gurobi, GLPK, CPLEX, HiGHS, BenchmarkTools, PrettyTables

function benchmark_methods(test_fn::Function)
    env = Gurobi.Env()
    bench = Dict{Tuple{Symbol,Symbol},Any}()
    optimizers = Dict(
        zip(
            [:Gurobi, :CPLEX, :GLPK, :HiGHS],
            [() -> Gurobi.Optimizer(env), CPLEX.Optimizer, GLPK.Optimizer, HiGHS.Optimizer],
        ),
    )
    methods = [:bisection, :binary, :echelon, :SOS1, :SOS2]

    global a , b , c , d
    d = test_fn
    for (sym, opt) in optimizers
        for m in methods
            a = sym
            b = opt
            c = m
            try
                bench[(a, c)] = @benchmark d(b, c)
            catch
                bench[(a, c)] = nothing
            end
        end
    end
    header = ["Optimizer"; string.(methods)]

    mean_(a) = typeof(a) === Nothing ? "N/A" : mean(a.times / 1000000)
    data = hcat(
        string.(keys(optimizers)),
        [mean_(bench[(a, b)]) for a in keys(optimizers), b in methods],
    )

    return pretty_table(data, header = header)
end

function quadratic_optimization(optimizer, method)
    model = JuMP.Model(optimizer)
    set_silent(model)

    @variable(model, -1 <= x <= 1)
    @variable(model, -1 <= y <= 1)

    @constraint(model, x + y == 0.5)

    @objective(
        model,
        Min,
        approximate(x, a -> a^2, method = method) +
        approximate(y, a -> a^2, method = method)
    )

    optimize!(model)
    return nothing
end

function interp_pts(optimizer, method)
    X = 0.5
    Y = 0.5
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
    model = JuMP.Model(optimizer)
    set_silent(model)

    @variable(model, x)
    @variable(model, y)

    z = PoGO.interpolate_points([x, y], points, method = method)
    @objective(model, Min, z)

    @constraint(model, x == X)
    @constraint(model, y == Y)

    optimize!(model)
    return objective_value(model)
end

# set number of points between knots to 4, and test methods and solvers
PoGO.set_default(:n, 4);
benchmark_methods(quadratic_optimization)

# set number of points between knots to 20, and test methods and solvers
PoGO.set_default(:n, 20);
benchmark_methods(quadratic_optimization)

# test methods and solvers for interpolation
benchmark_methods(interp_pts)
