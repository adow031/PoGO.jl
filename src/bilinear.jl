function bilinear(
    x::Union{VariableRef,AffExpr},
    y::Union{VariableRef,AffExpr},
    n::Int = 0;
    method::Symbol = :echelon,
    type::Symbol = :interior,
)
    if method == :convex
        error("Invalid method.")
    end

    type2 = nothing

    model = get_model(x)

    lx, ux = find_bounds(x)
    ly, uy = find_bounds(y)

    x_type, x_values = get_type(x)
    y_type, y_values = get_type(y)

    if n == 0
        if x_type == :binary && y_type == :binary
            xy = @variable(model, binary = true) # doesn't need to be binary, but it may help with branch-and-bound.
            @constraint(model, xy >= x + y - 1)
            @constraint(model, xy <= x)
            @constraint(model, xy <= y)
            return xy
        end

        if x_type == :binary
            xy = @variable(model)
            if ly >= 0
                @constraint(model, xy <= y)
            else
                @constraint(model, xy <= y - ly * (1 - x))
            end
            if uy <= 0
                @constraint(model, xy >= y)
            else
                @constraint(model, xy >= y - uy * (1 - x))
            end
            @constraint(model, xy <= uy * x)
            @constraint(model, xy >= ly * x)
            return xy
        end

        if y_type == :binary
            xy = @variable(model)
            if lx >= 0
                @constraint(model, xy <= x)
            else
                @constraint(model, xy <= x - lx * (1 - y))
            end
            if ux <= 0
                @constraint(model, xy >= x)
            else
                @constraint(model, xy >= x - ux * (1 - y))
            end
            @constraint(model, xy <= ux * y)
            @constraint(model, xy >= lx * y)
            return xy
        end

        if x_values == nothing
            x_values = collect(lx:ux)
        end
        if y_values == nothing
            y_values = collect(ly:uy)
        end

        if (x_type ∈ [:integer,:discrete] || y_type ∈ [:integer,:discrete])
            xy = @variable(model)

            mx = max(
                maximum(y_values) * maximum(x_values),
                minimum(y_values) * minimum(x_values),
                maximum(y_values) * minimum(x_values),
                minimum(y_values) * maximum(x_values),
            )
            mn = min(
                maximum(y_values) * maximum(x_values),
                minimum(y_values) * minimum(x_values),
                maximum(y_values) * minimum(x_values),
                minimum(y_values) * maximum(x_values),
            )

            set_lower_bound(xy, mn)
            set_upper_bound(xy, mx)

            ϵ = Dict()

            if x_type ∈ [:integer,:discrete] &&
               (y_type ∉ [:integer,:discrete] || length(y_values) >= length(x_values))
                if length(x_values) >= 50
                   @warn(
                       "The number of discrete levels being modelled ($(length(x_values))) within the bilinear function is large, consider using an approximation"
                   )
                end
                ϵ[0] = 0.0
                ϵ[length(x_values)] = 1.0
                for i in 1:length(x_values)-1
                    ϵ[i] = @variable(model, binary = true)
                end

                for i in 2:length(x_values)-1
                    @constraint(model, ϵ[i-1] <= ϵ[i])
                end
                @constraint(
                    model,
                    x == sum(x_values[i] * (ϵ[i] - ϵ[i-1]) for i in 1:length(x_values))
                )

                for i in 1:length(x_values)
                    @constraint(
                        model,
                        xy <= (mx - mn) * (1 - ϵ[i] + ϵ[i-1]) + x_values[i] * y
                    )
                    @constraint(
                        model,
                        xy >= -(mx - mn) * (1 - ϵ[i] + ϵ[i-1]) + x_values[i] * y
                    )
                end

                return xy
            else
                if length(y_values) >= 50
                   @warn(
                       "The number of discrete levels being modelled ($(length(y_values))) within the bilinear function is large, consider using an approximation"
                   )
                end
                ϵ[0] = 0.0
                ϵ[length(y_values)] = 1.0
                for i in 1:length(y_values)-1
                    ϵ[i] = @variable(model, binary = true)
                end

                for i in 2:length(y_values)-1
                    @constraint(model, ϵ[i-1] <= ϵ[i])
                end
                @constraint(
                    model,
                    y == sum(y_values[i] * (ϵ[i] - ϵ[i-1]) for i in 1:length(y_values))
                )

                for i in 1:length(y_values)
                    @constraint(
                        model,
                        xy <= (mx - mn) * (1 - ϵ[i] + ϵ[i-1]) + y_values[i] * x
                    )
                    @constraint(
                        model,
                        xy >= -(mx - mn) * (1 - ϵ[i] + ϵ[i-1]) + y_values[i] * x
                    )
                end

                return xy
            end
        end
    end

    if n <= 1
        error("n must be at least 2.")
    end

    if type == :interior || type == :tangent_cuts
        type2 = type
    elseif type == :lower
        type = :tangent_cuts
        type2 = :interior
    elseif type == :upper
        type = :interior
        type2 = :tangent_cuts
    else
        error("Invalid bound type.")
    end

    a² = @variable(model)
    b² = @variable(model)
    xy = @variable(model)

    set_lower_bound(xy, min(ux * ly, lx * uy))
    set_upper_bound(xy, max(ux * uy, lx * ly))

    a² = approximate(
        ((x - (lx + ux) / 2) / (ux - lx) + (y - (ly + uy) / 2) / (uy - ly)) / 2,
        x -> x^2,
        2 * n,
        method = method,
        type = type,
    )

    b² = approximate(
        ((x - (lx + ux) / 2) / (ux - lx) - (y - (ly + uy) / 2) / (uy - ly)) / 2,
        x -> x^2,
        2 * n,
        method = method,
        type = type2,
    )
    @constraint(
        model,
        xy - x * (ly + uy) / 2 - y * (lx + ux) / 2 + (lx + ux) * (ly + uy) / 4 ==
        (a² - b²) * (ux - lx) * (uy - ly)
    )

    return xy
end

function power(
    x::Union{VariableRef,AffExpr},
    y::Union{VariableRef,AffExpr},
    n::Int;
    method::Symbol = :echelon,
    type::Symbol = :interior,
)
    if method == :convex
        error("Invalid method.")
    end

    lx, ux = find_bounds(x,ignore_errors=true)

    if lx >= 1e-6
        z = approximate(x, a -> log(a), n, type = type, method = method)
        w = bilinear(y, z, n, method = method, type = type)
        v = approximate(w, a -> exp(a), n, type = type, method = method)
        return v
    elseif ux <= -1e-6
        y_type, y_values = get_type(y)
        if y_type ∈ [:binary,:integer]
            z = approximate(x, a -> log(-a), n, type = type, method = method)
            w = bilinear(y, z, n, method = method, type = type)
            v = approximate(w, a -> exp(a), n, type = type, method = method)
            u = approximate(y, a -> (-1)^round(Int,a))
            set_integer(u)
            set_lower_bound(u,-1.0)
            set_upper_bound(u,1.0)
            p = bilinear(u, v, n, method = method, type = type)
            return p
        else
            error("A negative number must have an integer or discrete exponent.")
        end
    else
        error("Currently PoGO.jl cannot model exponents of near-zero numbers.")
    end
end
