using JuMP, PoGO, Gurobi, Plots

struct Point
    x::Real
    y::Real
end

points = Point[]
push!(points, Point(-0.75, 0.5))
push!(points, Point(0, 1))
push!(points, Point(0.9, 0.8))
push!(points, Point(1.2, 0.5))
push!(points, Point(0.5, 0))

roughs = Point[]
push!(roughs, Point(0, 0))
push!(roughs, Point(0, 1))
push!(roughs, Point(1.5, 2))
push!(roughs, Point(1, 0))

model = JuMP.Model(Gurobi.Optimizer)
@variable(model, sx)
@variable(model, sy)
@variable(model, 0 <= m <= 2)
@variable(model, 0 <= θ <= 2π)
@variable(model, 0 <= α[1:length(points), 1:length(roughs)])

@variable(model, 0 <= px[1:length(points)] <= 2)
@variable(model, 0 <= py[1:length(points)] <= 2)

sinθ = approximate(θ, a -> sin(a), 50)
cosθ = approximate(θ, a -> cos(a), 50)

msinθ = bilinear(m, sinθ, 50)
mcosθ = bilinear(m, cosθ, 50)

for i in 1:length(points)
    @constraint(model, px[i] == sx + mcosθ * points[i].x - msinθ * points[i].y)
    @constraint(model, py[i] == sy + msinθ * points[i].x + mcosθ * points[i].y)
    @constraint(model, sum(α[i, j] for j in 1:length(roughs)) == 1)
    for j in 1:length(roughs)
        @constraint(model, px[i] == sum(α[i, j] * roughs[j].x for j in 1:length(roughs)))
        @constraint(model, py[i] == sum(α[i, j] * roughs[j].y for j in 1:length(roughs)))
    end
end

@objective(model, Max, m)
optimize!(model)

value.(θ)

Plots.plot(
    [value.(model[:px])[i] for i in [1:length(points); 1]],
    [value.(model[:py])[i] for i in [1:length(points); 1]],
    aspect_ratio = :equal,
)
Plots.plot!(
    [roughs[i].x for i in [1:length(roughs); 1]],
    [roughs[i].y for i in [1:length(roughs); 1]],
    aspect_ratio = :equal,
)
Plots.plot!(
    [points[i].x for i in [1:length(points); 1]],
    [points[i].y for i in [1:length(points); 1]],
    aspect_ratio = :equal,
)
