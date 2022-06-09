using JuMP, PoGO, Gurobi, Plots

model = JuMP.Model(Gurobi.Optimizer)

@variable(model, z, Bin)
@variable(model, 0 <= capacity <= 40)
@variable(model, x >= 0)
@variable(model, y >= 0)

@constraint(model, x <= bilinear(z, capacity))
@constraint(model, x + y == 50)

@objective(model, Min, x + 10 * y + 5 * capacity)

optimize!(model)

value(x)

println(model)

model = JuMP.Model(Gurobi.Optimizer)

@variable(model, 0 <= z <= 3, Int)
@variable(model, 0 <= capacity <= 40)
@variable(model, x >= 0)
@variable(model, y >= 0)

@constraint(model, x <= bilinear(z, capacity))
@constraint(model, x + y == 50)

@objective(model, Min, x + 10 * y + 5 * capacity + 100 * z)

optimize!(model)

value(x)
value(y)
value(z)
value(capacity)

objective_value(model)

println(model)

model = JuMP.Model(Gurobi.Optimizer)

@variable(model, 0 <= x <= 10, Int)
@variable(model, 0 <= y <= 10, Int)
@constraint(model, x + 2 * y == 10)

@objective(model, Max, bilinear(x, y))

optimize!(model)

value(x)
value(y)

model = JuMP.Model(Gurobi.Optimizer)

@variable(model, x, Bin)
@variable(model, y, Bin)

@objective(model, Max, 10 * bilinear(x, y) - 5 * x - 4 * y)

println(model)
optimize!(model)

value(x)
value(y)
