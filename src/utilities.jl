function get_model(x::Union{VariableRef,AffExpr})
    if typeof(x) == VariableRef
        return x.model
    else
        return collect(keys(x.terms))[1].model
    end
end

function find_bounds(x::Union{VariableRef,AffExpr})
    if typeof(x) == VariableRef
        if is_binary(x)
            lx = 0
            ux = 1
        else
            lx = lower_bound(x)
            ux = upper_bound(x)
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
    return lx, ux
end

function find_points(
    lx::Real,
    ux::Real,
    func::Function,
    n::Integer,
    knots::Union{Nothing,Vector{Float64}},
    convex::Bool,
    type::Symbol,
)
    i = 1
    f = 0.0
    g = 0.0
    single_output = false
    try
        f, g = func(lx)
    catch
        single_output = true
        if type ∈ [:lower, :upper, :tangent_cuts]
            error(
                "Function must return both a value and its derivative, if tangent cuts may be required.",
            )
        end
    end

    if knots != nothing
        while i <= length(knots)
            if knots[i] <= lx || knots[i] >= ux
                if knots[i] <= lx
                    convex = !convex
                end
                deleteat!(knots, i)
            else
                i += 1
            end
        end
    else
        knots = Float64[]
    end

    insert!(knots, 1, lx)
    push!(knots, ux)

    xi = Float64[]
    fxi = Float64[]

    x = lx
    if single_output
        f = func(x)
    else
        f, g = func(x)
    end

    if typeof(f) <: Real
        push!(xi, lx)
        push!(fxi, f)
    else
        push!(xi, lx)
        push!(fxi, f[1])
        push!(xi, lx)
        push!(fxi, f[2])
    end

    for j in 1:length(knots)-1
        if (!convex && type == :upper) ||
           (convex && type == :lower) ||
           type == :tangent_cuts
            for i in 1:n-1
                xold = x
                fold = f
                gold = g
                x = knots[j] + i * (knots[j+1] - knots[j]) / (n - 1)
                (f, g) = func(x)
                xreal = (f - fold + gold * xold - g * x) / (gold - g)
                freal = (xreal - x) * g + f
                push!(xi, xreal)
                push!(fxi, freal)
            end
        elseif (!convex && type == :lower) ||
               (convex && type == :upper) ||
               type == :interior
            for i in 1:n-1
                x = knots[j] + i * (knots[j+1] - knots[j]) / n
                if single_output
                    f = func(x)
                else
                    f, g = func(x)
                end
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
        if single_output
            f = func(x)
        else
            f, g = func(x)
        end
        if typeof(f) <: Real
            push!(xi, x)
            push!(fxi, f)
        else
            push!(xi, x)
            push!(fxi, f[1])
            push!(xi, x)
            push!(fxi, f[2])
        end
        convex = !convex
    end
    return xi, fxi
end
