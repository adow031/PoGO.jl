"""
    interpolate(
        x_vector::Vector{<:Union{VariableRef,AffExpr}},
        points::Vector{<:Tuple},
        sets::Vector;
        update_bounds::Bool = false,
        method::Symbol = :default,
        name::String = "",
    )

Function that constrains the variables to lie inside the convex hull of one of the sets of points
(defined by `sets`).

### Required arguments
`x_vector` vector of variables or expressions that will be used as the domain of the interpolation

`points` vector of points that define the interpolation; these points can contain real number and as
well as variables and expressions.

`sets` for each point in the vector `points`, `sets` contains the name(s) of the sets that point is a
member of. 

### Optional arguments
`update_bounds` is set to `true` if the upper and lower bounds of `x_vector` components should be
updated based on the values in `points`.

`method` if the formulation method and can be set to `:convex`, `:SOS1`, `:binary` or `:bisection`.

`name` can be set to give the variables created meaningful names.
"""
function interpolate(
    x_vector::Vector{<:Union{VariableRef,AffExpr}},
    points::Vector{<:Tuple},
    sets::Vector;
    update_bounds::Bool = false,
    method::Symbol = :default,
    name::String = "",
)
    if length(x_vector) != length(points[1])
        error(
            "interpolate() function must be provided with a vector of variables ('x_vector') which is the same length as the dimensions of the each of the points in 'points'.",
        )
    end

    methods = [:convex, :SOS1, :binary, :bisection]

    if method == :default
        method = Symbol(get(ENV, "POGO_METHOD", "binary"))
    end

    if method ∉ methods
        error("Invalid method.")
    end

    model = get_model(x_vector[1])

    if update_bounds
        l = Dict()
        u = Dict()
        for x in x_vector
            l[x] = Inf
            u[x] = -Inf
        end

        for p in points
            for i in eachindex(p)
                if typeof(p[i]) <: Real
                    if p[i] < l[x_vector[i]]
                        l[x_vector[i]] = p[i]
                    end
                    if p[i] > u[x_vector[i]]
                        u[x_vector[i]] = p[i]
                    end
                else
                    lp, up = find_bounds(p[i])
                    if lp < l[x_vector[i]]
                        l[x_vector[i]] = lp
                    end
                    if up > u[x_vector[i]]
                        u[x_vector[i]] = up
                    end
                end
            end
        end

        for x in x_vector
            ll, uu = find_bounds(x, ignore_errors = true)
            l[x] = max(l[x], ll)
            u[x] = min(u[x], uu)
            set_upper_bound(x, u[x])
            set_lower_bound(x, l[x])
        end
    end

    all_sets = Vector{Any}()
    for s in sets
        if typeof(s) <: Vector
            for i in s
                if i ∉ all_sets
                    push!(all_sets, i)
                end
            end
        else
            if s ∉ all_sets
                push!(all_sets, s)
            end
        end
    end

    y = Dict()
    z = Dict()
    α = Dict()

    for p in eachindex(points)
        α[p] = @variable(model)
        set_name(α[p], "α#$(p)")
        set_lower_bound(α[p], 0)
        set_upper_bound(α[p], 1) #perhaps only include if corresponding point is a variable.
    end

    if method == :binary
        for set in all_sets
            z[set] = @variable(model, binary = true)
            set_name(z[set], "z_set#$(set)")
        end
    else
        for set in all_sets
            z[set] = @variable(model)
            set_name(z[set], "z_set#$(set)")
            set_lower_bound(z[set], 0)
        end
    end

    if method == :SOS1
        @constraint(model, [z[set] for set in all_sets] in SOS1())
    elseif method == :bisection
        k = ceil(Int, log2(length(all_sets)))

        for j in 0:k-1
            y[j] = @variable(model, binary = true)
            if name == ""
                set_name(y[j], "y$(j)_[$(index(y[j]).value)]")
            else
                set_name(y[j], "y$(j)_$(name)")
            end
        end

        for i in 0:length(all_sets)-1
            total = AffExpr(1 - k)
            for j in 0:k-1
                if (i & 2^j) != 0
                    @constraint(model, z[all_sets[i+1]] <= y[j])
                    total += y[j]

                else
                    @constraint(model, z[all_sets[i+1]] <= 1 - y[j])
                    total += 1 - y[j]
                end
            end
            @constraint(model, z[all_sets[i+1]] >= total)
        end
        for i in length(all_sets):(2^k-1)
            total = AffExpr(1 - k)
            for j in 0:k-1
                total += (i & 2^j) != 0 ? y[j] : 1 - y[j]
            end
            @constraint(model, 0 >= total)
        end
    end

    for p in eachindex(points)
        if typeof(sets[p]) <: Vector
            @constraint(model, α[p] <= sum(z[j] for j in sets[p]))
        else
            @constraint(model, α[p] <= z[sets[p]])
        end
    end
    @constraint(model, sum(values(z)) == 1)

    for i in eachindex(x_vector)
        sum_p = sum(
            α[p] * points[p][i] for p in eachindex(points) if typeof(points[p][i]) <: Real;
            init = 0.0,
        )
        sum_p += sum(
            bilinear(α[p], points[p][i], n) for
            p in eachindex(points) if typeof(points[p][i]) <: Union{VariableRef,AffExpr};
            init = 0.0,
        )
        @constraint(model, x_vector[i] == sum_p)
    end

    @constraint(model, sum(α[p] for p in eachindex(points)) == 1)

    return z
end

"""
    interpolate_fn(
        f::Function,
        x_vector::Vector{<:Union{VariableRef,AffExpr}},
        grid::Vector{<:Union{AbstractRange,Vector{<:Real}}};
        method::Symbol = :default,
    )

Function that interpolates a function `f` over a set of variables (or affine expressions). The `grid`
vector should be the same length as the `x_vector`. 

### Required arguments
`f` a possibly multidimensional function that takes a vector of arguments the same length as `x_vector`.

`x_vector` a vector of variables or expressions that are the input vector for `f`.

`grid` is a vector of ranges or vectors that specify the sample points (which are the Cartesian product
of the vectors) for the function `f`.

### Optional arguments
`method` if the formulation method and can be set to `:convex`, `:SOS1`, `:binary` or `:bisection`.
"""
function interpolate_fn(
    f::Function,
    x_vector::Vector{<:Union{VariableRef,AffExpr}},
    grid::Vector{<:Union{AbstractRange,Vector{<:Real}}};
    method::Symbol = :default,
)
    points = collect(Iterators.product(grid...))
    if !(typeof(points) <: Matrix)
        points = collect(points)
    end
    points2 = convert(Matrix, transpose(hcat(collect.(vcat(points...))...)))
    return interpolate_fn(f, x_vector, points2; method = method)
end

"""
    interpolate_fn(
        f::Function,
        x_vector::Vector{<:Union{VariableRef,AffExpr}},
        points::T where {T<:Matrix};
        method::Symbol = :default,
    )

Function that interpolates a function `f` over a set of variables (or affine expressions). The number
of columns of the `points` matrix should be the same as the length of the `x_vector`. 

### Required arguments
`f` a possibly multidimensional function that takes a vector of arguments the same length as `x_vector`.

`x_vector` a vector of variables or expressions that are the input vector for `f`.

`points` is a `Matrix` that specifies the sample points for the function.

### Optional arguments
`method` if the formulation method and can be set to `:convex`, `:SOS1`, `:binary` or `:bisection`.
"""
function interpolate_fn(
    f::Function,
    x_vector::Vector{<:Union{VariableRef,AffExpr}},
    points::T where {T<:Matrix};
    method::Symbol = :default,
)
    y_matrix = nothing
    for i in 1:size(points)[1]
        y = f(points[i, :]...)
        if i == 1
            y_matrix = y'
        else
            y_matrix = [y_matrix; y']
        end
    end

    points2 = [points y_matrix]

    return interpolate_points(x_vector, points2; method = method)
end

"""
    interpolate_points(
        x_vector::Vector{<:Union{VariableRef,AffExpr}},
        points::T where {T<:Matrix};
        method::Symbol = :default,
    )

Function that interpolates over a set of variables (or affine expressions). The number of columns of
the `points` matrix should be greater than or equal to the length of the `x_vector`. For any extra
columns, new variables will be defined, and their values will be interpolated. 

### Required arguments
`x_vector` a vector of variables or expressions around which the triangulation will be formed.

`points` is a `Matrix` that specifies the points to be interpolated.

### Optional arguments
`method` if the formulation method and can be set to `:convex`, `:SOS1`, `:binary` or `:bisection`.
"""
function interpolate_points(
    x_vector::Vector{<:Union{VariableRef,AffExpr}},
    points::T where {T<:Matrix};
    method::Symbol = :default,
)
    mesh = delaunay(points[:, 1:length(x_vector)])
    points2 = Tuple[]
    for i in 1:size(points)[1]
        push!(points2, Tuple(points[i, :]))
    end

    sets = []
    for p in eachindex(points2)
        set = []
        for i in 1:size(mesh.simplices)[1]
            if p in mesh.simplices[i, :]
                push!(set, i)
            end
        end
        push!(sets, set)
    end

    z = []
    for i in length(x_vector):(size(points)[2]-1)
        var = @variable(get_model(x_vector[1]))
        push!(z, var)
        push!(x_vector, var)
    end

    interpolate(x_vector, points2, sets; method = method)
    return length(z) > 1 ? z : length(z) == 1 ? z[1] : nothing
end

function set_var_domain(x::Union{VariableRef,AffExpr}, domain::Vector)
    sets = []
    vals = Tuple[]

    for d in domain
        if typeof(d) == Int
            push!(vals, (d,))
            push!(sets, d)
        elseif typeof(d) == Vector{<:Real} || length(d) == 2
            push!(vals, (d[1],))
            push!(vals, (d[2],))
            push!(sets, Tuple(d))
            push!(sets, Tuple(d))
        else
            error("Invalid domain specified")
        end
    end

    return interpolate([x], vals, sets)
end
