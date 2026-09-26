#=
Meshes of the cylindric lens surfaces. The cross section lies in a local plane and is extruded
along the local x-axis, see `_extrude!`.
=#

"""Returns the angles of an arc from `-φ_max` to `φ_max` with steps of at most 3°."""
_arc_angles(φ_max) = LinRange(-φ_max, φ_max, max(17, ceil(Int, 2φ_max / _MAX_STEP) + 1))

# Cut disk in the local y-z-plane: the arc of `radius` above the chord z = h, extruded along x
function _tessellate(s::BMO.ConvexCylinderSDF)
    r, w = s.radius, s.diameter / 2
    h = sqrt(r^2 - w^2)
    φs = _arc_angles(asin(w / r))
    ws = [r * sin(φ) for φ in φs]
    lo = fill(h, length(φs))
    hi = [r * cos(φ) for φ in φs]
    lo[1] = lo[end] = hi[1] = hi[end] = h
    m = _TriMesh()
    _extrude!(m, ws, lo, hi, fill((0.0, -1.0), length(φs)), [(sin(φ), cos(φ)) for φ in φs],
        s.height / 2, (2, 3, 1))
    return _transform!(m, s)
end

# Box of the sagitta in the local z-y-plane minus the cylinder of `radius` along x at y = radius
function _tessellate(s::BMO.ConcaveCylinderSDF)
    R, w = s.radius, s.diameter / 2
    φs = _arc_angles(asin(w / abs(R)))
    ws = [abs(R) * sin(φ) for φ in φs]
    arc = [R - sign(R) * abs(R) * cos(φ) for φ in φs]
    flat = zeros(length(φs))
    m = _TriMesh()
    if R > 0
        _extrude!(m, ws, flat, arc, fill((0.0, -1.0), length(φs)), [(-sin(φ), cos(φ)) for φ in φs],
            s.height / 2, (3, 2, 1))
    else
        _extrude!(m, ws, arc, flat, [(-sin(φ), -cos(φ)) for φ in φs], fill((0.0, 1.0), length(φs)),
            s.height / 2, (3, 2, 1))
    end
    return _transform!(m, s)
end

# Aspheric cross section in the local z-y-plane, see `_aspheric_profile`, extruded along x
function _tessellate(acyl::BMO.AbstractAcylindricalSurfaceSDF)
    r = BMO.diameter(acyl) / 2
    ws = LinRange(-r, r, 100)
    lo, hi, nlo, nhi = _aspheric_profile(acyl, ws, acyl isa BMO.AconvexCylinderSDF)
    m = _TriMesh()
    _extrude!(m, ws, lo, hi, nlo, nhi, BMO.height(acyl) / 2, (3, 2, 1))
    return _transform!(m, acyl)
end

const _CylinderLensSDF = Union{BMO.ConvexCylinderSDF, BMO.ConcaveCylinderSDF,
    BMO.AbstractAcylindricalSurfaceSDF}

_has_mesh(::_CylinderLensSDF) = true

_default_color(::BMO.AbstractAcylindricalSurfaceSDF) = :white

function render!(
        axis::_RenderEnv,
        acyl::BMO.AbstractAcylindricalSurfaceSDF;
        # Makie kwargs
        color=:white,
        kwargs...
    )
    _render_mesh!(axis, acyl; color, kwargs...)
    return nothing
end
