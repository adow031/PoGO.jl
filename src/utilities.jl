function get_model(x::Union{VariableRef,AffExpr})
    if typeof(x) == VariableRef
        return x.model
    else
        return collect(keys(x.terms))[1].model
    end
end

function get_type(x::Union{VariableRef,AffExpr})
    if typeof(x) == VariableRef
        if is_binary(x)
            return :binary, nothing
        elseif is_integer(x)
            return :integer, nothing
        else
            return :continuous, nothing
        end
    else
        for (var, coeff) in x.terms
            if !(is_binary(var) || is_integer(var))
                return :continuous, nothing
            end
        end
        #x_type = abs(x.constant - round(x.constant)) < 10^-10 ? :integer : :continuous
        vals = Vector{Float64}[]
        count = 1
        for (var, coeff) in x.terms
            lx, ux = find_bounds(var)
            temp = collect(lx:ux)
            push!(vals, coeff * temp)
            count *= length(temp)
        end
        if count > 5000
            @warn(
                "Number of discrete levels for AffExpr exceeded 5000, skipping enumeration."
            )
            return :continuous, nothing
        end

        col = 1
        indices = ones(Int, length(vals))
        final = Float64[]
        while true
            temp = sum(vals[i][indices[i]] for i in eachindex(vals)) + x.constant
            if temp ∉ final
                push!(final, temp)
            end
            indices[col] += 1
            while indices[col] > length(vals[col])
                indices[col] = 1
                col += 1
                if col > length(vals)
                    return :discrete, final
                end
                indices[col] += 1
            end
            col = 1
        end
    end
end

function find_bounds(x::Union{VariableRef,AffExpr}; ignore_errors = false)
    if typeof(x) == VariableRef
        if is_binary(x)
            lx = 0
            ux = 1
        else
            if ignore_errors
                try
                    lx = lower_bound(x)
                catch
                    lx = -Inf
                end
                try
                    ux = upper_bound(x)
                catch
                    ux = Inf
                end
            else
                lx = lower_bound(x)
                ux = upper_bound(x)
            end
        end
    else
        lx = x.constant
        ux = x.constant
        for (var, coeff) in x.terms
            if is_binary(var)
                if coeff > 0
                    ux += coeff
                elseif coeff < 0
                    lx += coeff
                end
            else
                if ignore_errors
                    if coeff > 0
                        try
                            lx += coeff * lower_bound(var)
                        catch
                            lx -= Inf
                        end
                        try
                            ux += coeff * upper_bound(var)
                        catch
                            ux += Inf
                        end
                    elseif coeff < 0
                        try
                            lx += coeff * upper_bound(var)
                        catch
                            lx -= Inf
                        end
                        try
                            ux += coeff * lower_bound(var)
                        catch
                            ux += Inf
                        end
                    end
                else
                    if coeff > 0
                        lx += coeff * lower_bound(var)
                        ux += coeff * upper_bound(var)
                    elseif coeff < 0
                        lx += coeff * upper_bound(var)
                        ux += coeff * lower_bound(var)
                    end
                end
            end
        end
    end
    return lx, ux
end

function find_points(
    lx::Real,
    ux::Real,
    func::Function,
    n::Integer,
    δ::Real,
    knots::Vector{Float64},
    knots_shape::Vector{Symbol},
    type::Symbol,
)
    xi = Float64[]
    fxi = Float64[]

    x = knots[1]
    f = func(x)
    g = ForwardDiff.derivative(func, x)

    if typeof(f) <: Real
        push!(xi, x)
        push!(fxi, f)
    else
        push!(xi, x)
        push!(fxi, f[1])
        push!(xi, x)
        push!(fxi, f[2])
    end
    n1 = n
    for j in 1:length(knots)-1
        shape = knots_shape[j]
        if n1 == 0 && δ > 0.0
            n = max(2, ceil(Int, (knots[j+1] - knots[j]) / δ))
        end
        if (shape == :concave && type == :upper) ||
           (shape == :convex && type == :lower) ||
           (shape != :linear && type ∈ [:tangent_cuts, :combined])
            for i in 1:n-1
                xold = x
                if length(f) == 2
                    fold = f[2]
                else
                    fold = f
                end
                if length(g) == 2
                    gold = g[2]
                else
                    gold = g
                end

                x = knots[j] + i * (knots[j+1] - knots[j]) / (n - 1)

                f = func(x)
                g = ForwardDiff.derivative(func, x)
                if length(f) == 2
                    if length(g) == 2
                        xreal = (f[1] - fold + gold * xold - g[1] * x) / (gold - g[1])
                        freal = (xreal - x) * g[1] + f[1]
                    else
                        xreal = (f[1] - fold + gold * xold - g * x) / (gold - g)
                        freal = (xreal - x) * g + f[1]
                    end
                else
                    if length(g) == 2
                        xreal = (f - fold + gold * xold - g[1] * x) / (gold - g[1])
                        freal = (xreal - x) * g[1] + f
                    else
                        xreal = (f - fold + gold * xold - g * x) / (gold - g)
                        freal = (xreal - x) * g + f
                    end
                end
                push!(xi, xreal)
                push!(fxi, freal)
            end
            if type == :combined
                for i in 0:n-2
                    fxi[end-i] = (fxi[end-i] * 0.5 + 0.5 * func(xi[end-i]))
                end
            end
        elseif (shape == :concave && type == :lower) ||
               (shape == :convex && type == :upper) ||
               (shape != :linear && type == :interior)
            for i in 1:n-1
                x = knots[j] + i * (knots[j+1] - knots[j]) / n
                f = func(x)
                g = ForwardDiff.derivative(func, x)
                if typeof(f) <: Real
                    push!(xi, x)
                    push!(fxi, f)
                else
                    push!(xi, x)
                    push!(fxi, f[1])
                    push!(xi, x)
                    push!(fxi, f[2])
                end
            end
        end
        x = knots[j+1]
        f = func(x)
        g = ForwardDiff.derivative(func, x)
        if typeof(f) <: Real
            push!(xi, x)
            push!(fxi, f)
        else
            push!(xi, x)
            push!(fxi, f[1])
            push!(xi, x)
            push!(fxi, f[2])
        end
    end
    return xi, fxi
end

function process_knots(
    knots::Union{Vector{Tuple{Float64,Symbol}}},
    lx::Real,
    ux::Real,
)::Tuple{Vector{Float64},Vector{Symbol}}
    knots_shape = nothing
    i = 1

    prev_x = -Inf
    prev_shape = knots[1][2]
    while knots[i][1] < lx
        prev_x = knots[i][1]
        prev_shape = knots[i][2]
        deleteat!(knots, i)
    end
    i = length(knots)
    while knots[i][1] > ux
        deleteat!(knots, i)
        i -= 1
    end
    knots_shape = [knots[i][2] for i in eachindex(knots)]
    knots = [knots[i][1] for i in eachindex(knots)]
    if knots[1] > lx
        insert!(knots, 1, lx)
        insert!(knots_shape, 1, prev_shape)
    end
    if knots[end] < ux
        push!(knots, ux)
        push!(knots_shape, knots_shape[end])
    end
    return knots, knots_shape
end

function infer_curvature(
    f::Function,
    knots::Union{Nothing,Vector{Float64},Vector{Tuple{Float64,Symbol}}},
    lx::Real,
    ux::Real,
)
    function get_curve(f::Function, m::Real)
        curvature = ForwardDiff.derivative(x -> ForwardDiff.derivative(f, x), m)
        if curvature > 0
            return :convex
        elseif curvature == 0
            return :linear
        else
            return :concave
        end
    end

    c = Symbol[]

    if knots === nothing
        knots = [float(lx), float(ux)]
    end

    if typeof(knots) == Vector{Float64}
        if lx < knots[1]
            insert!(knots, 1, lx)
        end
        if ux > knots[end]
            push!(knots, ux)
        end
        for i in eachindex(knots)
            if i == length(knots)
                m = knots[i] + 1e-5
            else
                m = (knots[i] + knots[i+1]) / 2
            end
            push!(c, get_curve(f, m))
        end
        return [(knots[i], c[i]) for i in eachindex(knots)]
    end

    knots = copy(knots)

    if lx < knots[1][1]
        m = (lx + knots[1][1]) / 2
        insert!(knots, 1, (lx, get_curve(f, m)))
    end
    if ux > knots[end][1]
        m = ux + 1e-5
        push!(knots, (ux, get_curve(f, m)))
    end
end
