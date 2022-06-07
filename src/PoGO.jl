module PoGO

using JuMP

function get_model(x::Union{VariableRef,AffExpr})
    if typeof(x) == VariableRef
        return x.model
    else
        return collect(keys(x.terms))[1].model
    end
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

    # xi = Float64[]
    # fxi = Float64[]
    # for i in 0:n
    #     push!(xi, lx + i * (ux - lx) / n)
    #     push!(fxi, func(xi[end]))
    # end

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
        # elseif method==:echelon2
        #         ϵ = Dict()
        #         z = Dict()
        #         for i in 0:n-1
        #             ϵ[i] = @variable(model, binary = true)
        #             z[i+1] = @variable(model)
        #             set_lower_bound(z[i+1],0)
        #         end
        #
        #         if target != nothing
        #             order = PoGO.find_order(fxi[1:end-1],target)
        #         else
        #             order = 1:n
        #         end
        #
        #         for i in 0:n-2
        #             @constraint(model, z[order[i+1]] == ϵ[i+1]-ϵ[i])
        #         end
        #
        #         @constraint(model,z[order[n]]==1-ϵ[n-1])
        #
        #         for i in 1:n-1
        #             @constraint(model,α[i]<=z[i]+z[i+1])
        #         end
        #
        #         @constraint(model,sum(z[i] for i in 1:n)==1)
        #
        #         @constraint(model,α[0]<=z[1])
        #         @constraint(model,α[n]<=z[n])
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
            #@constraint(model, [z[i] for i in 0:n-1] in SOS1([i+1 for i in reorder(n)]))
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

function find_bounds(x::Union{VariableRef,AffExpr})
    if typeof(x) == VariableRef
        lx = lower_bound(x)
        ux = upper_bound(x)
    else
        lx = x.constant
        ux = x.constant
        for (var, coeff) in x.terms
            if coeff > 0
                lx += coeff * lower_bound(var)
                ux += coeff * upper_bound(var)
            elseif coeff < 0
                lx += coeff * upper_bound(var)
                ux += coeff * lower_bound(var)
            end
        end
    end
    return lx, ux
end

function bilinear(
    x::Union{VariableRef,AffExpr},
    y::Union{VariableRef,AffExpr},
    n::Int;
    method::Symbol = :echelon,
    type::Symbol = :interior,
)
    if method == :convex
        error("Invalid method.")
    end

    type2 = nothing

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

    model = get_model(x)

    a² = @variable(model)
    b² = @variable(model)
    xy = @variable(model)

    lx, ux = find_bounds(x)
    ly, uy = find_bounds(y)

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

export bilinear, approximate
end

# n=5
# (x1,y1)=find_points(0,2π,a->(cos(a),-sin(a)),n,[float(π/2),float(3π/2)],false,:lower)
#
# (x2,y2)=find_points(0,2π,a->(cos(a),-sin(a)),n,[float(π/2),float(3π/2)],false,:upper)
#
#     Plots.plot(
#         x1,
#         y1,
#         seriestype = :line,
#         title = "Particle Equilibrium",
#     )
#     Plots.plot!(
#         x2,
#         y2,
#         seriestype = :line,
#         title = "Particle Equilibrium",
#     )
#     Plots.plot!(
#         xx,
#         yy,
#         seriestype = :line,
#         title = "Particle Equilibrium",
#     )
#
# xx=[2π/100*i for i in 0:100]
# yy=[cos(2π/100*i) for i in 0:100]
#
# func = a->a^2
#
# func(3)
#
# try
#     f,g = func(5)
# catch
#     func = a -> [func(a),0.0]
# end
# [a,b]=func(3)
# function find_order(a::Vector{Float64},target::Union{Real,Symbol})
#     if typeof(target)==Symbol
#         if target == :max
#             order=[i for (v,i) in sort([(maximum(a)-a[i],i) for i in 1:length(a)])]
#         elseif target == :min
#             order=[i for (v,i) in sort([(a[i]-minimum(a),i) for i in 1:length(a)])]
#         end
#     else
#         order=[i for (v,i) in sort([(abs(target-a[i]),i) for i in 1:length(a)])]
#     end
#     order2=Int[]
#     toggle=false
#     for i in 1:length(order)
#         if toggle
#             push!(order2,order[i])
#         else
#             insert!(order2,1,order[i])
#         end
#         toggle=!toggle
#     end
#     return order2
# end
# function reorder(n)
#     a=collect(n ÷ 2:-2:1)
#     b=collect(n ÷ 2+2:2:n)
#
#     c=reverse(collect(n ÷ 2+1:2:n))
#     d=reverse(collect(n ÷ 2-1:-2:1))
#
#     order = [a[1]]
#     for i in 1:length(b)
#         push!(order,b[i])
#         try
#             push!(order,a[i+1])
#         catch
#         end
#     end
#     push!(order,c[1])
#     for i in 1:length(d)
#         push!(order,d[i])
#         try
#             push!(order,c[i+1])
#         catch
#         end
#     end
#     return order
# end
