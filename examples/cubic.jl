using JuMP, PoGO, Gurobi, Plots

function cubic(rhs::Real, n::Int, type::Symbol)
    model = JuMP.Model(Gurobi.Optimizer)
    @variable(model, -1 <= x <= 2.5)

    @constraint(model, x == rhs)

    @objective(
        model,
        Max,
        approximate(
            x,
            a -> [a^3 - 3a^2 + a + 2, 3a^2 - 6a + 1],
            n,
            type = type,
            initial = :concave,
            knots = [1.0],
        )
    )
    optimize!(model)
    return objective_value(model)
end

function approximate_cubic(; step = 0.05, n = 10, type = :interior)
    x = collect(-1:step:2.5)
    y = Float64[]
    for i in x
        println(i)
        push!(y, cubic(i, n, type))
    end

    Plots.plot(x, y)
    return Plots.plot!(x, [i^3 - 3i^2 + i + 2 for i in x])
end

approximate_cubic(step = 0.05, n = 4, type = :interior)
approximate_cubic(step = 0.01, n = 4, type = :tangent_cuts)
approximate_cubic(step = 0.01, n = 4, type = :lower)
approximate_cubic(step = 0.01, n = 4, type = :upper)
