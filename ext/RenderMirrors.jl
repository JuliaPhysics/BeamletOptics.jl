"""
    render!(ax, s::ConicSDF; color=:silver, kwargs...)

Analytical surface renderer for a [`ConicSDF`](@ref).
Renders the concave/convex conic front reflective face, the cylindrical substrate side wall,
and the flat rear substrate base into `ax`.
"""
function render!(
        ax::_RenderEnv,
        s::BMO.ConicSDF;
        color = :silver,
        kwargs...
)
    r_max = s.diameter / 2
    N_r, N_theta, N_u = 50, 100, 30

    r_grid = LinRange(0, r_max, N_r)
    theta_grid = LinRange(0, 2π, N_theta)
    u_grid = LinRange(0, 1, N_u)

    R_rot = BMO.orientation(s)
    P_pos = BMO.position(s)

    R = 2 * s.f
    Z_off = BMO._conic_sag(s.x_off, R, s.k)
    y_surf(xl, zl) = -(BMO._conic_sag(sqrt((xl + s.x_off)^2 + zl^2), R, s.k) - Z_off)

    function draw_surface!(X_loc, Y_loc, Z_loc; surf_color = color)
        Xt = Float32.(R_rot[1, 1] .* X_loc .+ R_rot[1, 2] .* Y_loc .+
                      R_rot[1, 3] .* Z_loc .+ P_pos[1])
        Yt = Float32.(R_rot[2, 1] .* X_loc .+ R_rot[2, 2] .* Y_loc .+
                      R_rot[2, 3] .* Z_loc .+ P_pos[2])
        Zt = Float32.(R_rot[3, 1] .* X_loc .+ R_rot[3, 2] .* Y_loc .+
                      R_rot[3, 3] .* Z_loc .+ P_pos[3])
        surface!(ax, Xt, Yt, Zt; colormap = [surf_color, surf_color], kwargs...)
    end

    # Front Conic Surface
    X_front = [r * cos(t) for r in r_grid, t in theta_grid]
    Z_front = [r * sin(t) for r in r_grid, t in theta_grid]
    Y_front = [min(s.thickness, y_surf(r * cos(t), r * sin(t)))
               for r in r_grid, t in theta_grid]
    draw_surface!(X_front, Y_front, Z_front; surf_color = color)

    # Substrate Cylindrical Side Wall
    X_wall = [r_max * cos(t) for u in u_grid, t in theta_grid]
    Z_wall = [r_max * sin(t) for u in u_grid, t in theta_grid]
    Y_rim = [min(s.thickness, y_surf(r_max * cos(t), r_max * sin(t))) for t in theta_grid]
    Y_wall = [(1 - u) * Y_rim[j] + u * s.thickness
              for (i, u) in enumerate(u_grid), (j, t) in enumerate(theta_grid)]
    draw_surface!(X_wall, Y_wall, Z_wall; surf_color = :grey)

    # Substrate Rear Flat Base
    X_back = [r * cos(t) for r in r_grid, t in theta_grid]
    Z_back = [r * sin(t) for r in r_grid, t in theta_grid]
    Y_back = [s.thickness for r in r_grid, t in theta_grid]
    draw_surface!(X_back, Y_back, Z_back; surf_color = :grey)

    return nothing
end
