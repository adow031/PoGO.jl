using JuMP, PoGO, Gurobi, Plots

function globaloptimal(n::Int, type::Symbol)
    model = JuMP.Model(Gurobi.Optimizer)
    @variable(model, -1 <= x <= 2.5)
    @variable(model, 0 <= y <= 2π)

    fx = approximate(x, a -> a^3 - 3a^2 + a + 6, n, type = type, knots = [1.0])

    fy = approximate(y, a -> sin(a) + 2, n, type = type, knots = [float(π)])

    @objective(model, Max, bilinear(fx, fy, n, type = type))
    optimize!(model)

    return objective_value(model), objective_bound(model)
end

solution1 = []
solution2 = []

for n in 2:4:100
    push!(solution1, globaloptimal(n, :lower))
    push!(solution2, globaloptimal(n, :upper))
end

Plots.plot(1:length(solution1), [s[2] for s in solution1])
Plots.plot!(1:length(solution2), [s[1] for s in solution2])
