using JuMP, PoGO, Gurobi

struct Tranche
    q::Float64
    p::Float64
end

function define_offerstack(q::Vector{<:Real}, p::Vector{<:Real})
    @assert(length(q) == length(p))
    return [Tranche(float(q[t]), float(p[t])) for t in 1:length(q)]
end

struct Line
    capacity::Float64
    reactance::Float64
    resistance::Float64
end

function Line(capacity::Real, reactance::Real, resistance::Real)
    return Line(float(capacity), float(reactance), float(resistance))
end

function offercost(x::Real, tranches::Vector{Tranche})
    totalq = 0.0
    totalc = 0.0
    for t in tranches
        if x <= totalq + t.q
            return totalc + t.p * (x - totalq)
        end
        totalq += t.q
        totalc += t.q * t.p
    end
    return totalc
end

function dispatch(nodes, G, offers, demand, lines, n_loss_tranches)
    model = JuMP.Model(Gurobi.Optimizer)

    @variable(model, 0 <= x[g in keys(G)] <= sum([o.q for o in offers[g]]))
    @variable(model, -lines[l].capacity <= flow[l in keys(lines)] <= lines[l].capacity)
    @variable(model, θ[nodes])

    @constraint(
        model,
        nodebalance[n in nodes],
        sum(x[g] for (g, m) in G if m == n; init = 0) + sum(
            flow[(i, j)] - approximate(
                flow[(i, j)],
                x -> 0.5 * lines[(i, j)].resistance * x^2,
                n_loss_tranches,
                knots = [0.0],
            ) for (i, j) in keys(lines) if n == j;
            init = 0,
        ) - sum(
            flow[(i, j)] + approximate(
                flow[(i, j)],
                x -> 0.5 * lines[(i, j)].resistance * x^2,
                n_loss_tranches,
                knots = [0.0],
            ) for (i, j) in keys(lines) if n == i;
            init = 0,
        ) == demand[n]
    )

    @constraint(model, θ[nodes[1]] == 0)
    @constraint(
        model,
        voltageangles[(i, j) in keys(lines)],
        lines[(i, j)].reactance * (θ[i] - θ[j]) == flow[(i, j)]
    )

    @objective(
        model,
        Min,
        sum(
            approximate(
                x[g],
                a -> offercost(a, offers[g]),
                knots = [sum(offers[g][τ].q for τ in 1:t) for t in 1:length(offers[g])],
                method = :convex,
            ) for g in keys(G)
        )
    )

    optimize!(model)

    vars = all_variables(model)
    solution = Dict(zip(vars, value.(vars)))

    for v in vars
        if is_binary(v)
            unset_binary(v)
            fix(v, solution[v])
        end
    end

    optimize!(model)
    for g in keys(G)
        println("Node $(g) dispatch: $(value(x[g]))MW")
    end
    for n in nodes
        println("Node $(n) price: \$$(dual(nodebalance[n]))/MWh")
    end
    for l in keys(lines)
        println("Line $(l) flow: $(value(flow[l]))MW")
    end
    #value.(flow)
end

nodes = [:a, :b, :c]
offers = Dict{Symbol,Vector{Tranche}}()
demand = Dict{Symbol,Real}()

lines = Dict{Tuple{Symbol,Symbol},Line}()

G = Dict{Symbol,Symbol}()
G[:x] = :a
G[:y] = :b

offers[:x] = define_offerstack([2, 1, 1], [1, 2, 3])
offers[:y] = define_offerstack([1, 1, 2], [4, 10, 8])
#offers[:c] = define_offerstack([0,0,0],[3,6,9])

demand[:a] = 2
demand[:b] = 2
demand[:c] = 0.1

lines[(:a, :b)] = Line(10, 1, 0.01)
lines[(:b, :c)] = Line(0.2, 1, 0.01)
lines[(:a, :c)] = Line(10, 1, 0.01)

dispatch(nodes, G, offers, demand, lines, 1)
