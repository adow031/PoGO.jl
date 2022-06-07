using JuMP, PoGO, Gurobi, Plots

function particle(p::Int, n::Int, method::Symbol)
    model = JuMP.Model(Gurobi.Optimizer)

    @variable(model, 0 <= x[1:p] <= 10)
    @variable(model, 0 <= y[1:p] <= 10)
    @variable(model, mindist)

    for i in 1:p-1
        for j in i+1:p
            temp =
                approximate(x[i] - x[j], a -> a^2, n, method = method) +
                approximate(y[i] - y[j], a -> a^2, n, method = method)
            @constraint(model, temp >= mindist)
            @constraint(model, temp <= 25)
        end
        @constraint(model, x[i] <= x[i+1])
    end

    @constraint(model, y[2] - y[1] == 3 * (x[2] - x[1]))
    @objective(model, Max, mindist)
    optimize!(model)

    println(value.(model[:x]))
    println(value.(model[:y]))

    return Plots.plot(
        value.(model[:x]),
        value.(model[:y]),
        seriestype = :scatter,
        title = "Particle Equilibrium",
    )
end

particle(5, 10, :SOS1)
particle(7, 10, :SOS2)
particle(5, 10, :binary)
particle(7, 10, :echelon)
