using Panco, Test, JuMP, GLPK

function nonlinear(n::Int)
    model = JuMP.Model(GLPK.Optimizer)
    @variable(model, 0 <= x <= 4)
    @variable(model, 0 <= y <= 4)
    @constraint(model, approximate(x, a -> a^2, n) + approximate(y, a -> a^2, n) >= 9)

    @objective(
        model,
        Max,
        bilinear(
            approximate(x, a -> exp(-a / 10) * sin(a), n),
            approximate(y, a -> exp(-a / 5) * sin(2 * a), n),
            n,
        )
    )

    optimize!(model)

    value(model[:y])

    return [value(model[:x]), value(model[:y]), objective_value(model)]
end

@testset "Panco.jl" begin
    result = nonlinear(10)
    @test sum(abs.(result - [1.6, 4.0, 0.37866])) ≈ 0.0 atol = 1e-4

    result = nonlinear(20)
    @test sum(abs.(result - [1.4, 3.8, 0.3878])) ≈ 0.0 atol = 1e-4

    result = nonlinear(50)
    @test sum(abs.(result - [1.44, 3.84, 0.3922])) ≈ 0.0 atol = 1e-4
end
