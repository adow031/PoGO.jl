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
            return :integer, nothing
        end

        col = 1
        indices = ones(Int, length(vals))
        final = Float64[]
        while true
            temp = sum(vals[i][indices[i]] for i in 1:length(vals)) + x.constant
            if temp ∉ final
                push!(final, temp)
            end
            indices[col] += 1
            while indices[col] > length(vals[col])
                indices[col] = 1
                col += 1
                if col > length(vals)
                    return :integer, final
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
    knots_shape::Union{Vector{Nothing},Vector{Symbol}},
    type::Symbol,
)
    f = 0.0
    g = 0.0

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
           (shape != :linear && type == :tangent_cuts)
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
                    xreal = (f - fold + gold * xold - g * x) / (gold - g)
                    freal = (xreal - x) * g + f
                end
                push!(xi, xreal)
                push!(fxi, freal)
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
    knots::Union{Vector{Float64},Vector{Tuple{Float64,Symbol}}},
    lx::Real,
    ux::Real,
    initial::Union{Nothing,Symbol} = nothing,
)
    knots_shape = nothing
    i = 1
    if typeof(knots) == Vector{Tuple{Float64,Symbol}}
        prev_x = -Inf
        prev_shape = knots[1][2]
        while i <= length(knots)
            if knots[i][1] <= lx || knots[i][1] >= ux
                if knots[i][1] <= lx && knots[i][1] > prev_x
                    prev_x = knots[i][1]
                    prev_shape = knots[i][2]
                end
                deleteat!(knots, i)
            else
                i += 1
            end
        end
        knots_shape = [knots[i][2] for i in 1:length(knots)]
        knots = [knots[i][1] for i in 1:length(knots)]
        insert!(knots, 1, lx)
        push!(knots, ux)
        insert!(knots_shape, 1, prev_shape)
        push!(knots_shape, knots_shape[end])

    else
        while i <= length(knots)
            if knots[i] <= lx || knots[i] >= ux
                deleteat!(knots, i)
            else
                i += 1
            end
        end
        insert!(knots, 1, lx)
        push!(knots, ux)

        knots_shape = [initial]
        for i in 2:length(knots)
            if knots_shape[i-1] == :concave
                push!(knots_shape, :convex)
            else
                push!(knots_shape, :concave)
            end
        end
    end

    return knots, knots_shape
end
