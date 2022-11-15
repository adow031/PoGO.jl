using JuMP, PoGO

# Uncomment appropriate solver

using Gurobi
optimizer =
    optimizer_with_attributes(() -> Gurobi.Optimizer(), "OutputFlag" => 0, "MIPGap" => 0.0)

# using CPLEX
# optimizer = optimizer_with_attributes(
#     CPLEX.Optimizer,
#     "CPXPARAM_MIP_Tolerances_MIPGap" => 0.0,
#     "CPX_PARAM_SCRIND" => 0,
# )
#
# using GLPK
# optimizer = optimizer_with_attributes(
#     (method = GLPK.INTERIOR) -> GLPK.Optimizer(),
#     "msg_lev" => 0,
#     "mip_gap" => 0.0,
# )
#
# using Cbc
# optimizer = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0, "ratioGap" => 0.0)

function cubic_plot(; δ = 2, type = :interior)
    return plot_approximation(
        optimizer,
        a -> a^3 - 3a^2 + a + 2,
        -1.0,
        2.5,
        δ,
        type,
        detail = 40,
        knots = [(-1.0, :concave), (1.0, :convex)],
    )
end

cubic_plot(δ = 0.5, type = :interior)
cubic_plot(δ = 0.5, type = :tangent_cuts)
cubic_plot(δ = 0.5, type = :lower)
cubic_plot(δ = 0.5, type = :upper)

function waves(; δ = 0.1, type = :interior, wave = :sine)
    knots = nothing

    if wave == :square
        wave_fn = function (x)
            if abs(x % 2) == 1
                return [-1, 1]
            elseif x % 2 == 0
                return [1, -1]
            elseif abs(floor(x)) % 2 == 1
                return 1
            elseif floor(x) % 2 == 0
                return -1
            else
                return 0
            end
        end
        knots = [(i, :linear) for i in -4.0:4.0]
    elseif wave == :saw
        wave_fn = function (x)
            if abs(x % 1) == 0
                return [1, 0]
            else
                return x - floor(x)
            end
        end
        knots = [(i, :linear) for i in -4.0:4.0]
    elseif wave == :sine
        wave_fn = x -> sin(x)
        knots = [(float(-π), :convex), (0.0, :concave), (float(π), :convex)]
    elseif wave == :other
        wave_fn = function (x)
            if abs(x % 2) == 1
                return [-x, x]
            elseif x % 2 == 0
                return [x, -x]
            elseif abs(floor(x)) % 2 == 1
                return x
            elseif floor(x) % 2 == 0
                return -x
            else
                return 0
            end
        end
        knots = [(i, :linear) for i in -4.0:4.0]
    elseif wave == :test
        wave_fn = function (x)
            if x < 1
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
        knots = [
            (-4.0, :concave),
            (-0.01, :linear),
            (0.01, :concave),
            (1.0, :linear),
            (2.0, :convex),
        ]
    end

    return plot_approximation(
        optimizer,
        wave_fn,
        -4.5,
        4.5,
        δ,
        type,
        detail = 50,
        knots = knots,
    )
end

waves(δ = 1, type = :tangent_cuts, wave = :square)
waves(δ = 1, type = :tangent_cuts, wave = :saw)
waves(δ = 1, type = :tangent_cuts, wave = :other)
waves(δ = 0.8, type = :tangent_cuts, wave = :sine)
waves(δ = 0.8, type = :interior, wave = :sine)
waves(δ = 0.2, type = :tangent_cuts, wave = :sine)
waves(δ = 0.2, type = :interior, wave = :sine)
waves(δ = 0.8, type = :tangent_cuts, wave = :test)

function general_fn(; δ = 0.1, type = :interior)
    fn = function (x)
        if x < 1
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
    knots = [
        (-1.0, :concave),
        (-0.01, :linear),
        (0.01, :concave),
        (1.0, :linear),
        (2.0, :convex),
    ]

    return plot_approximation(optimizer, fn, -1.0, 3.0, δ, type, detail = 30, knots = knots)
end

general_fn(δ = 0.5, type = :interior)
general_fn(δ = 0.5, type = :lower)
general_fn(δ = 0.5, type = :upper)
general_fn(δ = 0.5, type = :tangent_cuts)
