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

    if n == 0
        if typeof(x) == VariableRef &&
           is_binary(x) &&
           typeof(y) == VariableRef &&
           is_binary(y)
            xy = @variable(model, binary = true) # doesn't need to be binary, but it may help with branch-and-bound.
            @constraint(model, xy >= x + y - 1)
            @constraint(model, xy <= x)
            @constraint(model, xy <= y)
            return xy
        end

        if typeof(x) == VariableRef && is_binary(x)
            xy = @variable(model)
            @constraint(model, xy <= (uy - ly) * x + y)
            @constraint(model, xy >= -(uy - ly) * x + y)
            @constraint(model, xy <= uy * (1 - x))
            @constraint(model, xy >= -uy * (1 - x))
            return xy
        end

        if typeof(y) == VariableRef && is_binary(y)
            xy = @variable(model)
            @constraint(model, xy <= (ux - lx) * y + x)
            @constraint(model, xy >= -(ux - lx) * y + x)
            @constraint(model, xy <= ux * (1 - y))
            @constraint(model, xy >= -ux * (1 - y))
            return xy
        end

        if typeof(x) == VariableRef && is_integer(x)
            if typeof(y) != VariableRef || !is_integer(y) || uy - ly >= ux - lx
                xy = @variable(model)

                ϵ = Dict()
                for i in lx:ux
                    ϵ[i] = @variable(model, binary = true)
                end
                for i in (lx+1):ux
                    @constraint(model, ϵ[i-1] <= ϵ[i])
                end
                @constraint(model, ϵ[ux] == 1)
                @constraint(
                    model,
                    x == lx * ϵ[lx] + sum(i * (ϵ[i] - ϵ[i-1]) for i in (lx+1):ux)
                )

                mx = max(uy * ux, ly * lx, uy * ux, ly * lx)
                mn = min(uy * ux, ly * lx, uy * ux, ly * lx)

                for i in (lx+1):ux
                    @constraint(model, xy <= (mx - mn) * (1 - ϵ[i] + ϵ[i-1]) + i * y)
                    @constraint(model, xy >= -(mx - mn) * (1 - ϵ[i] + ϵ[i-1]) + i * y)
                end

                @constraint(model, xy <= (mx - mn) * (1 - ϵ[lx]) + lx * y)

                return xy
            end
        end

        if typeof(y) == VariableRef && is_integer(y)
            xy = @variable(model)

            ϵ = Dict()
            for i in ly:uy
                ϵ[i] = @variable(model, binary = true)
            end
            for i in (ly+1):uy
                @constraint(model, ϵ[i-1] <= ϵ[i])
            end

            @constraint(
                model,
                y == ly * ϵ[lx] + sum(i * (ϵ[i] - ϵ[i-1]) for i in (ly+1):uy)
            )

            mx = max(uy * ux, ly * lx, uy * ux, ly * lx)
            mn = min(uy * ux, ly * lx, uy * ux, ly * lx)

            for i in (ly+1):uy
                @constraint(model, xy <= (mx - mn) * (1 - ϵ[i] + ϵ[i-1]) + i * x)
                @constraint(model, xy >= -(mx - mn) * (1 - ϵ[i] + ϵ[i-1]) + i * x)
            end
            return xy
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
        x -> [x^2, 2x],
        2 * n,
        method = method,
        type = type,
    )

    b² = approximate(
        ((x - (lx + ux) / 2) / (ux - lx) - (y - (ly + uy) / 2) / (uy - ly)) / 2,
        x -> [x^2, 2x],
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

function exponential(
    x::Union{VariableRef,AffExpr},
    y::Union{VariableRef,AffExpr},
    n::Int;
    method::Symbol = :echelon,
    type::Symbol = :interior,
)
    if method == :convex
        error("Invalid method.")
    end

    if lower_bound(x) < 10^-6
        error("Currently PoGO.jl cannot model exponents of zero or negative numbers.")
    end

    z = approximate(
        x,
        a -> [log(a), 1 / a],
        n,
        type = type,
        initial = :concave,
        method = method,
    )
    w = bilinear(y, z, n, method = method, type = type)
    v = approximate(
        w,
        a -> [exp(a), exp(a)],
        n,
        type = type,
        initial = :convex,
        method = method,
    )
    return v
end
