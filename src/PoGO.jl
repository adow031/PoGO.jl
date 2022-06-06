module PoGO

using JuMP

function get_model(x::Union{VariableRef,AffExpr})
    if typeof(x) == VariableRef
        return x.model
    else
        return collect(keys(x.terms))[1].model
    end
end

function approximate(
    x::Union{VariableRef,AffExpr},
    func::Function,
    n::Int;
    method::Symbol = :echelon,
)
    methods = [:convex, :SOS1, :SOS2, :binary, :echelon, :echelon2, :bisection]

    if method ∉ methods
        error("Invalid method.")
    end

    model = get_model(x)

    fx = @variable(model)

    lx, ux = find_bounds(x)

    xi = []
    fxi = []
    for i in 0:n
        push!(xi, lx + i * (ux - lx) / n)
        push!(fxi, func(xi[end]))
    end

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
        #         order = reorder(n)
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
)
    if method == :convex
        error("Invalid method.")
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
    )
    b² = approximate(
        ((x - (lx + ux) / 2) / (ux - lx) - (y - (ly + uy) / 2) / (uy - ly)) / 2,
        x -> x^2,
        2 * n,
        method = method,
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
