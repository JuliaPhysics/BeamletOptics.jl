function render!(
        axis::_RenderEnv,
        acyl::BMO.AbstractAcylindricalSurfaceSDF;
        # Makie kwargs
        color=:white,
        transparency=true,
        kwargs...
    )
    r = BMO.diameter(acyl) / 2
    w_vals = LinRange(-r, r, 100)           # aperture coordinate
    z_vals = LinRange(-BMO.height(acyl) / 2, BMO.height(acyl) / 2, 20)  # extrusion coordinate

    # Parameterize lateral surface as (w, aspheric_equation(w), z)
    # but swap the w and z directions to fix the 90° rotation.
    # That is, let X_local come from z_vals and Z_local from w_vals.
    X_local = [z for w in w_vals, z in z_vals]
    Y_local = [BMO.aspheric_equation(w, acyl) for w in w_vals, z in z_vals]
    Z_local = [w for w in w_vals, z in z_vals]

    # Transform to global coordinates:
    R = BMO.orientation(acyl)
    P = BMO.position(acyl)
    Xt = R[1, 1] .* X_local .+ R[1, 2] .* Y_local .+ R[1, 3] .* Z_local .+ P[1]
    Yt = R[2, 1] .* X_local .+ R[2, 2] .* Y_local .+ R[2, 3] .* Z_local .+ P[2]
    Zt = R[3, 1] .* X_local .+ R[3, 2] .* Y_local .+ R[3, 3] .* Z_local .+ P[3]

    surface!(axis, Xt, Yt, Zt; transparency, colormap=[color, color], kwargs...)
    render_acylindric_caps!(axis, acyl; color)
    return nothing
end

function render_acylindric_cap!(
        axis::_RenderEnv,
        acyl::BMO.AbstractAcylindricalSurfaceSDF,
        top::Bool;
        # Makie kwargs
        color = :white,
        transparency = true,
        kwargs...
    )
    r = BMO.diameter(acyl) / 2
    h = BMO.height(acyl)
    Xval = top ? (h / 2) : -(h / 2)   # local X coordinate for top or bottom cap
    Nw = 60                     # resolution along w
    Ny = 60                     # resolution along y

    w_vals = LinRange(-r, r, Nw)
    X_local = fill(Xval, Nw, Ny)
    Y_local = similar(X_local)
    Z_local = similar(X_local)
    # aspheric value at the edge
    awe = BMO.aspheric_equation(r, acyl)
    for i in 1:Nw
        w = w_vals[i]
        aw = BMO.aspheric_equation(w, acyl)
        # If aw>0, sweep from awe to aw; if aw<0, sweep from aw to awe
        y_range = aw >= 0 ? LinRange(awe, aw, Ny) : LinRange(aw, 0, Ny)
        for j in 1:Ny
            X_local[i, j] = Xval
            Y_local[i, j] = y_range[j]
            Z_local[i, j] = w
        end
    end

    # Transform local -> global
    R = BMO.orientation(acyl)
    P = BMO.position(acyl)
    Xt = R[1, 1] .* X_local .+ R[1, 2] .* Y_local .+ R[1, 3] .* Z_local .+ P[1]
    Yt = R[2, 1] .* X_local .+ R[2, 2] .* Y_local .+ R[2, 3] .* Z_local .+ P[2]
    Zt = R[3, 1] .* X_local .+ R[3, 2] .* Y_local .+ R[3, 3] .* Z_local .+ P[3]

    surface!(axis, Xt, Yt, Zt; transparency, colormap=[color, color], kwargs...)
    return nothing
end

function render_acylindric_caps!(axis::_RenderEnv, acyl::BMO.AbstractAcylindricalSurfaceSDF; kwargs...)
    render_acylindric_cap!(axis, acyl, true; kwargs...)  # top
    render_acylindric_cap!(axis, acyl, false; kwargs...)  # bottom
    return nothing
end