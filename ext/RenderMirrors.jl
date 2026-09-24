"""
    _tessellate_conic(s::ConicSDF, r_hole)

Analytic mesh of a [`ConicSDF`](@ref) substrate, optionally pierced by an axial cylindrical bore
of radius `r_hole` (`r_hole = 0` for the unpierced case). Consists of the concave/convex conic
front reflective face (annular if `r_hole > 0`), the cylindrical substrate side wall, the flat
rear substrate base, and, if `r_hole > 0`, the inner bore wall. All faces but the front face are
tagged as substrate (grey), see `_vertex_colors`.
"""
function _tessellate_conic(s::BMO.ConicSDF, r_hole)
    r_max = s.diameter / 2
    N_r = 50
    θs = [2π * (j - 1) / _N_THETA for j in 1:_N_THETA]

    R = 2 * s.f
    Z_off = BMO._conic_sag(s.x_off, R, s.k)
    y_surf(xl, zl) = min(s.thickness, -(BMO._conic_sag(sqrt((xl + s.x_off)^2 + zl^2), R, s.k) - Z_off))
    # outward normal of the front face y = y_surf(x, z), the substrate lies at y ≥ y_surf
    function n_surf(xl, zl)
        r_p = sqrt((xl + s.x_off)^2 + zl^2)
        slope = r_p > 0 ? BMO._conic_slope(r_p, R, s.k) / r_p : zero(r_p)
        return (-slope * (xl + s.x_off), -1, -slope * zl)
    end
    substrate = 0x01
    m = _TriMesh()
    ring(r, y, n; tag = 0x00) = r ≤ 0 ? [_add_vertex!(m, (0, y(0, 0), 0), n(0, 0, 0), tag)] :
        [_add_vertex!(m, (r * cos(θ), y(r * cos(θ), r * sin(θ)), r * sin(θ)),
            n(r * cos(θ), r * sin(θ), θ), tag) for θ in θs]

    # 1. Front conic surface (from r_hole to r_max)
    _add_rows!(m, [ring(r, y_surf, (x, z, _) -> n_surf(x, z)) for r in LinRange(r_hole, r_max, N_r)])
    # 2. Substrate cylindrical side wall
    radial(x, z, θ) = (cos(θ), 0, sin(θ))
    _add_rows!(m, [ring(r_max, y_surf, radial; tag = substrate),
        ring(r_max, (_, _) -> s.thickness, radial; tag = substrate)])
    # 3. Inner bore wall (only if pierced)
    if r_hole > 0
        inward(x, z, θ) = (-cos(θ), 0, -sin(θ))
        _add_rows!(m, [ring(r_hole, y_surf, inward; tag = substrate),
            ring(r_hole, (_, _) -> s.thickness, inward; tag = substrate)])
    end
    # 4. Substrate rear (annular) base
    back(x, z, θ) = (0, 1, 0)
    _add_rows!(m, [ring(r_hole, (_, _) -> s.thickness, back; tag = substrate),
        ring(r_max, (_, _) -> s.thickness, back; tag = substrate)])

    return _transform!(m, s)
end

_has_mesh(::BMO.ConicSDF) = true

_tessellate(s::BMO.ConicSDF) = _tessellate_conic(s, zero(s.diameter))

_default_color(::BMO.ConicSDF) = :silver

"""
    render!(ax, s::ConicSDF; color=:silver, kwargs...)

Analytical mesh renderer for a [`ConicSDF`](@ref).
Renders the concave/convex conic front reflective face in `color`, the cylindrical substrate side
wall and the flat rear substrate base in grey into `ax`.
"""
function render!(
        ax::_RenderEnv,
        s::BMO.ConicSDF;
        color = :silver,
        kwargs...
)
    return _render_mesh!(ax, s; color, kwargs...)
end

"""
    _axial_bore_radius(d::DifferenceSDF{T, <:ConicSDF})

Returns the radius of `d`'s tool if it is a single [`CylinderSDF`](@ref) coaxial with the
[`ConicSDF`](@ref) base and centered on its axis, or `nothing` otherwise.
"""
function _axial_bore_radius(d::BMO.DifferenceSDF{T, <:BMO.ConicSDF}) where {T}
    s = d.base
    if length(d.tools) == 1 && d.tools[1] isa BMO.CylinderSDF
        tool = d.tools[1]
        dot_axis = abs(dot(BMO.orientation(s)[:, 2], BMO.orientation(tool)[:, 2]))
        rel_pos = BMO.position(tool) - BMO.position(s)
        x_proj = dot(rel_pos, BMO.orientation(s)[:, 1])
        z_proj = dot(rel_pos, BMO.orientation(s)[:, 3])
        if isapprox(dot_axis, 1; atol = 1e-3) && abs(x_proj) < 1e-4 && abs(z_proj) < 1e-4
            return tool.radius
        end
    end
    return nothing
end

# An axial bore is part of the conic mesh, other tools are cut via `_difference`
function _tessellate_difference(::BMO.ConicSDF, d::BMO.DifferenceSDF)
    r_hole = _axial_bore_radius(d)
    r_hole === nothing || return _tessellate_conic(d.base, r_hole)
    return _difference((_tessellate(d.base), d.base), [(_tessellate(t), t) for t in d.tools])
end

_default_color(::BMO.DifferenceSDF{T, <:BMO.ConicSDF}) where {T} = :silver

"""
    render!(ax, d::DifferenceSDF{T, <:ConicSDF}; color=:silver, kwargs...)

Analytical mesh renderer for a [`ConicSDF`](@ref) with an axial cylindrical bore
(e.g. Cassegrain, Ritchey-Chrétien, or OAP with `:collimated` through-hole).
Other tools are cut from the mesh of the [`ConicSDF`](@ref), see `_difference`.
"""
function render!(
        ax::_RenderEnv,
        d::BMO.DifferenceSDF{T, <:BMO.ConicSDF};
        color = :silver,
        kwargs...
) where {T}
    return _render_mesh!(ax, d; color, kwargs...)
end
