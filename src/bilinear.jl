function one_binary_bilinear_formulation(
    model::JuMP.Model,
    xy::VariableRef,
    x::Union{VariableRef,AffExpr},
    y::Union{VariableRef,AffExpr},
)
    ly, uy = find_bounds(y)
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
    set_lower_bound(xy, ly)
    set_upper_bound(xy, uy)
    return xy
end

function discrete_bilinear_formulation(
    model::JuMP.Model,
    xy::VariableRef,
    x::Union{VariableRef,AffExpr},
    y::Union{VariableRef,AffExpr},
    x_values::Vector,
    mx::Float64,
    mn::Float64,
)
    ϵ = Dict()

    if length(x_values) >= 50
        @warn(
            "The number of discrete levels being modelled ($(length(x_values))) within the bilinear function is large, consider using an approximation"
        )
    end
    ϵ[0] = 0.0
    ϵ[length(x_values)] = 1.0
    for i in 1:length(x_values)-1
        ϵ[i] = @variable(model, binary = true)
        if name == ""
            set_name(ϵ[i], "ϵ$(i)_[$(index(ϵ[i]).value)]")
        else
            set_name(ϵ[i], "ϵ$(i)_$name")
        end
    end

    for i in 2:length(x_values)-1
        @constraint(model, ϵ[i-1] <= ϵ[i])
    end
    @constraint(model, x == sum(x_values[i] * (ϵ[i] - ϵ[i-1]) for i in eachindex(x_values)))

    for i in eachindex(x_values)
        @constraint(model, xy <= (mx - mn) * (1 - ϵ[i] + ϵ[i-1]) + x_values[i] * y)
        @constraint(model, xy >= -(mx - mn) * (1 - ϵ[i] + ϵ[i-1]) + x_values[i] * y)
    end

    return xy
end

function bilinear(
    x::Union{VariableRef,AffExpr},
    y::Union{VariableRef,AffExpr},
    n::Int = 0;
    method::Symbol = :default,
    type::Symbol = :interior,
    name::String = "",
)
    if method == :convex
        error("Invalid method.")
    end

    type2 = nothing

    model = get_model(x)

    x_type, x_values = get_type(x)
    y_type, y_values = get_type(y)

    if x_type == :binary && y_type == :binary
        xy = @variable(model, binary = true) # doesn't need to be binary, but it may help with branch-and-bound.
        if name == ""
            set_name(xy, "xy_[$(index(xy).value)]")
        else
            set_name(xy, name)
        end
        @constraint(model, xy >= x + y - 1)
        @constraint(model, xy <= x)
        @constraint(model, xy <= y)
        return xy
    end

    if x_type == :binary || y_type == :binary
        xy = @variable(model)
        if name == ""
            set_name(xy, "xy_[$(index(xy).value)]")
        else
            set_name(xy, name)
        end
        if x_type == :binary
            return one_binary_bilinear_formulation(model, xy, x, y)
        else
            return one_binary_bilinear_formulation(model, xy, y, x)
        end
    end

    lx, ux = find_bounds(x)
    ly, uy = find_bounds(y)

    if x_values === nothing
        if x_type == :integer
            x_values = collect(lx:ux)
        else
            x_values = [lx, ux]
        end
    end
    if y_values === nothing
        if y_type == :integer
            y_values = collect(ly:uy)
        else
            y_values = [ly, uy]
        end
    end

    if n == 0
        if (x_type ∈ [:integer, :discrete] || y_type ∈ [:integer, :discrete])
            xy = @variable(model)
            if name == ""
                set_name(xy, "xy_[$(index(xy).value)]")
            else
                set_name(xy, name)
            end
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

            if x_type ∈ [:integer, :discrete] &&
               (y_type ∉ [:integer, :discrete] || length(y_values) >= length(x_values))
                return discrete_bilinear_formulation(model, xy, x, y, x_values, mx, mn)
            else
                return discrete_bilinear_formulation(model, xy, y, x, y_values, mx, mn)
            end
        end

        n = parse(Int, get(ENV, "POGO_N", "10"))
    elseif n == 1
        error("n must be at least 2.")
    end

    if type ∈ [:interior, :tangent_cuts, :combined]
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

    xy = @variable(model)
    if name == ""
        set_name(xy, "xy_[$(index(xy).value)]")
    else
        set_name(xy, name)
    end

    set_lower_bound(xy, min(ux * ly, lx * uy, ux * uy, lx * ly))
    set_upper_bound(xy, max(ux * uy, lx * ly, ux * ly, lx * uy))

    a² = approximate(
        ((x - (lx + ux) / 2) / (ux - lx) + (y - (ly + uy) / 2) / (uy - ly)) / 2,
        x -> x^2,
        n,
        knots = [0.0],
        method = method,
        type = type,
        name = name,
    )

    b² = approximate(
        ((x - (lx + ux) / 2) / (ux - lx) - (y - (ly + uy) / 2) / (uy - ly)) / 2,
        x -> x^2,
        n,
        knots = [0.0],
        method = method,
        type = type2,
        name = name,
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
    n::Int = 0;
    method::Symbol = :default,
    type::Symbol = :interior,
    name::String = "",
)
    if method == :convex
        error("Invalid method.")
    end

    lx, ux = find_bounds(x, ignore_errors = true)

    if lx >= 1e-6
        z = approximate(x, a -> log(a), n, type = type, method = method)
        w = bilinear(y, z, n, method = method, type = type)
        v = approximate(w, a -> exp(a), n, type = type, method = method)
        return v
    elseif ux <= -1e-6
        y_type, y_values = get_type(y)
        if y_type ∈ [:binary, :integer]
            z = approximate(x, a -> log(-a), n, type = type, method = method)
            w = bilinear(y, z, 0, method = method, type = type)
            v = approximate(w, a -> exp(a), n, type = type, method = method)
            u = approximate(y, a -> (-1)^round(Int, a))
            set_integer(u)
            set_lower_bound(u, -1.0)
            set_upper_bound(u, 1.0)
            p = bilinear(u, v, 0, method = method, type = type, name = name)
            return p
        else
            error("A negative variable must have an integer or discrete exponent.")
        end
    else
        error("PoGO.jl cannot model exponents of near-zero numbers.")
    end
end
