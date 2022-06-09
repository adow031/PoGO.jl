using JuMP, PoGO

# Uncomment appropriate solver

# using Gurobi
# optimizer =
#     optimizer_with_attributes(() -> Gurobi.Optimizer(), "OutputFlag" => 0, "MIPGap" => 0.0)
#
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

function cubic_plot(; n = 2, type = :interior)
    return plot_approximation(
        optimizer,
        a -> [a^3 - 3a^2 + a + 2, 3a^2 - 6a + 1],
        :concave,
        -1.0,
        2.5,
        n,
        type,
        detail = 20,
        knots = [1.0],
    )
end

cubic_plot(n = 4, type = :interior)
cubic_plot(n = 4, type = :tangent_cuts)
cubic_plot(n = 4, type = :lower)
cubic_plot(n = 4, type = :upper)
