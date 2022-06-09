using JuMP, PoGO, Gurobi, Plots

function power_test(n::Int, type::Symbol)
    model = JuMP.Model(Gurobi.Optimizer)
    @variable(model, 0.1 <= x <= 3)
    @variable(model, 0.1 <= y <= 3)
    @constraint(model, x + y == 4)

    @objective(model, Max, PoGO.power(x, y, n; type = type))
    optimize!(model)

    println("x: $(value(x))")
    println("y: $(value(y))")
    println("approximate xy: $(objective_value(model))")
    println("actual xy: $(value(x)^value(y))")
    return objective_value(model)
end

gap = Float64[]

for i in 4:2:50
    push!(gap, power_test(i, :upper) - power_test(i, :lower))
end

Plots.plot(4:2:50, gap)
