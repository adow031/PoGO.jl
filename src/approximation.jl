function approximate(
    x::Union{VariableRef,AffExpr},
    func::Function,
    n::Int;
    method::Symbol = :echelon,
    type::Symbol = :interior,
    initial::Union{Nothing,Symbol} = nothing,
    knots::Union{Nothing,Vector{Float64}} = nothing,
)
    methods = [:convex, :SOS1, :SOS2, :binary, :echelon]

    if method ∉ methods
        error("Invalid method.")
    end

    if type != nothing && type ∉ [:lower, :upper, :interior, :tangent_cuts]
        error("Invalid approximation bound type.")
    elseif type ∈ [:lower, :upper] && initial == nothing
        error("Must specify initial curvature if type is :lower or :upper.")
    elseif type ∈ [:interior, :cutting_planes] && initial != nothing
        @warn(
            "Bound type is :interior or :tangent_cuts, so initial curvature will be ignored."
        )
    end

    model = get_model(x)

    fx = @variable(model)

    lx, ux = find_bounds(x)

    if knots != nothing
        knots = copy(knots)
    end
    xi, fxi = find_points(lx, ux, func, n, knots, initial == :convex, type)

    n = length(xi) - 1

    set_lower_bound(fx, minimum(fxi))
    set_upper_bound(fx, maximum(fxi))

    α = Dict()

    for i in 0:n
        α[i] = @variable(model)
        set_lower_bound(α[i], 0)
        set_upper_bound(α[i], 1)
    end

    @constraint(model, sum(α[i] for i in 0:n) == 1)
    @constraint(model, sum(α[i] * xi[i+1] for i in 0:n) == x)
    @constraint(model, sum(α[i] * fxi[i+1] for i in 0:n) == fx)

    if method == :echelon
        ϵ = Dict()
        for i in 0:n-1
            ϵ[i] = @variable(model, binary = true)
        end
        for i in 2:n-1
            @constraint(model, α[i] <= ϵ[i] - ϵ[i-2])
        end
        for i in 0:n-2
            @constraint(model, ϵ[i] <= ϵ[i+1])
        end

        @constraint(model, α[0] <= ϵ[0])
        @constraint(model, α[1] <= ϵ[1])
        @constraint(model, α[n] <= 1 - ϵ[n-2])

    elseif method == :binary || method == :SOS1
        z = Dict()

        if method == :binary
            for i in 0:n-1
                z[i] = @variable(model, binary = true)
            end
        else
            for i in 0:n-1
                z[i] = @variable(model)
                set_lower_bound(z[i], 0)
            end
            @constraint(model, [z[i] for i in 0:n-1] in SOS1())
        end

        @constraint(model, sum(z[i] for i in 0:n-1) == 1)

        for i in 1:n-1
            @constraint(model, α[i] <= z[i-1] + z[i])
        end

        @constraint(model, α[0] <= z[0])
        @constraint(model, α[n] <= z[n-1])
    elseif method == :SOS2
        @constraint(model, [α[i] for i in 0:n] in SOS2())
    end

    return fx
end
