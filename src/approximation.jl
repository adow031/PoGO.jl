function approximate(
    x::Union{VariableRef,AffExpr},
    func::Function,
    n::Int = 0;
    δ::Real = 0.0,
    method::Symbol = :echelon,
    type::Symbol = :interior,
    knots::Union{Nothing,Vector{Float64},Vector{Tuple{Float64,Symbol}}} = nothing,
)
    methods = [:convex, :SOS1, :SOS2, :binary, :echelon]

    if method ∉ methods
        error("Invalid method.")
    end

    if type !== nothing && type ∉ [:lower, :upper, :interior, :tangent_cuts, :combined]
        error("Invalid approximation bound type.")
    end

    lx, ux = find_bounds(x)
    model = get_model(x)
    fx = @variable(model)
    x_type, x_values = get_type(x)
    method2 = :continuous

    if x_type == :binary || x_type == :integer || x_type == :discrete
        method2 = :discrete
        if x_values == nothing
            x_values = collect(ceil(lx - 0.01):floor(ux + 0.01))
        end
        knots = [float(i) for i in x_values]
        knots_shape = [:linear for i in x_values]
    else
        knots = infer_curvature(func, knots, lx, ux)
        knots, knots_shape = process_knots(knots, lx, ux)
    end

    xi, fxi = find_points(lx, ux, func, n, δ, knots, knots_shape, type)
    n = length(xi) - 1
    set_lower_bound(fx, minimum(fxi))
    set_upper_bound(fx, maximum(fxi))

    α = Dict()

    for i in 0:n
        α[i] = @variable(model)
    end

    @constraint(model, sum(α[i] for i in 0:n) == 1)
    @constraint(model, sum(α[i] * xi[i+1] for i in 0:n) == x)
    @constraint(model, sum(α[i] * fxi[i+1] for i in 0:n) == fx)

    if method2 == :discrete
        if method == :SOS1
            @constraint(model, [α[i] for i in 0:n] in SOS1())
        else
            for i in 0:n
                set_binary(α[i])
            end
        end
    else
        for i in 0:n
            set_lower_bound(α[i], 0)
            set_upper_bound(α[i], 1)
        end

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
            if n > 1
                @constraint(model, α[1] <= ϵ[1])
                @constraint(model, α[n] <= 1 - ϵ[n-2])
            end
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
    end
    return fx
end
