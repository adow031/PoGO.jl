function xy_in_sets(
    x::Union{VariableRef,AffExpr},
    y::Union{VariableRef,AffExpr},
    sets::Any,
    n = 0;
    update_bounds = false,
)
    model = get_model(x)

    lx, ux = find_bounds(x, ignore_errors = true)
    ly, uy = find_bounds(y, ignore_errors = true)

    lx2 = Inf
    ux2 = -Inf
    ly2 = Inf
    uy2 = -Inf

    for set in sets
        for p in set
            if typeof(p[1]) <: Real
                if p[1] < lx2
                    lx2 = p[1]
                end
                if p[1] > ux2
                    ux2 = p[1]
                end
            else
                lp, up = find_bounds(p[1])
                if lp < lx2
                    lx2 = lp
                end
                if up > ux2
                    ux2 = up
                end
            end

            if typeof(p[2]) <: Real
                if p[2] < ly2
                    ly2 = p[2]
                end
                if p[2] > uy2
                    uy2 = p[2]
                end
            else
                lp, up = find_bounds(p[2])
                if lp < ly2
                    ly2 = lp
                end
                if up > uy2
                    uy2 = up
                end
            end
        end
    end

    lx = max(lx, lx2)
    ly = max(ly, ly2)
    ux = min(ux, ux2)
    uy = min(uy, uy2)

    if update_bounds
        set_upper_bound(x, ux)
        set_upper_bound(y, uy)
        set_lower_bound(x, lx)
        set_lower_bound(y, ly)
    end

    z = Dict()
    for set in sets
        M = Inf
        alpha = Dict()
        z[set] = @variable(model, binary = true)
        for p in set
            alpha[p] = @variable(model)
            set_lower_bound(alpha[p], 0.0)
            set_upper_bound(alpha[p], 1.0)
            if typeof(p[1]) <: Real
                lp1 = p[1]
                up1 = p[1]
            else
                lp1, up1 = find_bounds(p[1])
            end
            if typeof(p[2]) <: Real
                lp2 = p[2]
                up2 = p[2]
            else
                lp2, up2 = find_bounds(p[2])
            end

            newM = max(up1 - lx, ux - lp1, up2 - ly, uy - lp2)
            if newM < M
                M = newM
            end
        end

        sum_p = sum(alpha[p] * p[1] for p in set if typeof(p[1]) <: Real; init = 0.0)
        sum_p += sum(
            bilinear(alpha[p], p[1], n) for
            p in values(set) if typeof(p[1]) <: Union{VariableRef,AffExpr};
            init = 0.0,
        )
        @constraint(model, x <= sum_p + M * (1 - z[set]))
        @constraint(model, x >= sum_p - M * (1 - z[set]))

        sum_p = sum(alpha[p] * p[2] for p in set if typeof(p[2]) <: Real; init = 0.0)
        sum_p += sum(
            bilinear(alpha[p], p[2], n) for
            p in values(set) if typeof(p[2]) <: Union{VariableRef,AffExpr};
            init = 0.0,
        )
        @constraint(model, y <= sum_p + M * (1 - z[set]))
        @constraint(model, y >= sum_p - M * (1 - z[set]))

        @constraint(model, sum(alpha[p] for p in set) == 1)
    end
    @constraint(model, sum(values(z)) == 1)

    return z
end
