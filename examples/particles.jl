using JuMP, Panco, Gurobi, Plots

function particle(p::Int, n::Int)
    model = JuMP.Model(Gurobi.Optimizer)

    @variable(model, 0 <= x[1:p] <= 10)
    @variable(model, 0 <= y[1:p] <= 10)
    @variable(model, mindist)

    for i in 1:p-1
        for j in i+1:p
            @constraint(
                model,
                bilinear(x[i] - x[j], x[i] - x[j], n) +
                bilinear(y[i] - y[j], y[i] - y[j], n) >= mindist
            )
            @constraint(
                model,
                approximate(x[i] - x[j], a -> a^2, n) +
                approximate(y[i] - y[j], a -> a^2, n) <= 25
            )
        end
        @constraint(model, x[i] <= x[i+1])
    end

    @constraint(model, y[2] - y[1] == 3 * (x[2] - x[1]))
    #@constraint(model,y[1]<=y[2])
    @objective(model, Max, mindist)
    optimize!(model)
    value.(model[:x])
    value.(model[:y])
    return plot(
        value.(model[:x]),
        value.(model[:y]),
        seriestype = :scatter,
        title = "Particle Equilibrium",
    )
end

particle(7, 10)
