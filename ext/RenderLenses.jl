#=
Meshes of the rotationally symmetric lens surfaces. The profiles lie in the local r-y-plane and
are revolved around the local y-axis, see `_revolve!`.
=#

"""Radii of the rings of a curved lens surface with outer radius `r`, equal area per ring."""
_lens_radii(r) = sqrt.(LinRange(0, 1, _N_RADIAL + 1)) .* r

# Convex cap between the spherical face (y = 0 on axis) and the flat base at y = sag
function _tessellate(css::BMO.ConvexSphericalSurfaceSDF)
    R = BMO.radius(css)
    face = [(r, R - sqrt(R^2 - r^2), r / R, -sqrt(R^2 - r^2) / R) for r in _lens_radii(BMO.diameter(css) / 2)]
    m = _TriMesh()
    _revolve!(m, face)
    _revolve!(m, _disc(face[end][2], face[end][1], 1))
    return _transform!(m, css)
end

# Concave face (y = 0 on axis, y = -sag at the rim), rim wall and flat base at y = 0
function _tessellate(css::BMO.ConcaveSphericalSurfaceSDF)
    R = BMO.radius(css)
    face = [(r, sqrt(R^2 - r^2) - R, -r / R, -sqrt(R^2 - r^2) / R) for r in _lens_radii(BMO.diameter(css) / 2)]
    r_rim, y_rim = face[end]
    m = _TriMesh()
    _revolve!(m, face)
    _revolve!(m, _wall(r_rim, y_rim, 0.0, 1))
    _revolve!(m, _disc(0.0, r_rim, 1))
    return _transform!(m, css)
end

"""
    _aspheric_bounds(s, convex)

Returns the lower and upper boundary `(lo, hi)` of the solid of the aspheric (or acylindric)
surface SDF `s` along its local y-axis, each either `:face` for the aspheric face or the constant
height of the closing flat face, following the cases of `convex_aspheric_surface_distance` and
`concave_aspheric_surface_distance`.
"""
function _aspheric_bounds(s, convex::Bool)
    c = 1 / s.radius
    zb = BMO.aspheric_equation(BMO.diameter(s) / 2, s)
    ms = s.max_sag[1]
    if convex
        c > 0 && zb < 0 && return (:face, ms)
        return c > 0 ? (:face, zb) : (zb, :face)
    else
        ms > 0 && zb < 0 && return (zb, :face)
        return zb < 0 ? (:face, 0.0) : (0.0, :face)
    end
end

"""
    _aspheric_profile(s, ws, convex)

Returns the lower and upper boundary heights and their outward 2D normals at the samples `ws`, see
`_aspheric_bounds`.
"""
function _aspheric_profile(s, ws, convex::Bool)
    f(w) = BMO.aspheric_equation(w, s)
    δ = 1e-7 * BMO.diameter(s)
    function slope(w)
        d = (f(w + δ) - f(w - δ)) / 2δ
        isnan(d) && (d = (f(w) - f(w - δ)) / δ)
        return d
    end
    bounds = _aspheric_bounds(s, convex)
    lo, hi = [[b === :face ? f(w) : b for w in ws] for b in bounds]
    face_normal(w, up) = (up ? (-slope(w), 1.0) : (slope(w), -1.0)) ./ sqrt(1 + slope(w)^2)
    nlo = [bounds[1] === :face ? face_normal(w, false) : (0.0, -1.0) for w in ws]
    nhi = [bounds[2] === :face ? face_normal(w, true) : (0.0, 1.0) for w in ws]
    return lo, hi, nlo, nhi
end

function _tessellate(asp::BMO.AbstractAsphericalSurfaceSDF)
    r = BMO.diameter(asp) / 2
    rs = LinRange(0, r, 50)
    convex = asp isa BMO.ConvexAsphericalSurfaceSDF
    lo, hi, nlo, nhi = _aspheric_profile(asp, rs, convex)
    blo, bhi = _aspheric_bounds(asp, convex)
    m = _TriMesh()
    # flat faces with the minimal number of rings
    _revolve!(m, blo === :face ? [(rs[i], lo[i], nlo[i]...) for i in eachindex(rs)] : _disc(blo, r, -1))
    _revolve!(m, bhi === :face ? [(rs[i], hi[i], nhi[i]...) for i in eachindex(rs)] : _disc(bhi, r, 1))
    abs(hi[end] - lo[end]) > 1e-12 * r && _revolve!(m, _wall(r, lo[end], hi[end], 1))
    return _transform!(m, asp)
end

const _LensSurfaceSDF = Union{BMO.ConvexSphericalSurfaceSDF, BMO.ConcaveSphericalSurfaceSDF,
    BMO.AbstractAsphericalSurfaceSDF}

_has_mesh(::_LensSurfaceSDF) = true

_default_color(::_LensSurfaceSDF) = :white

function render!(
        axis::_RenderEnv,
        css::BMO.ConcaveSphericalSurfaceSDF;
        # Makie kwargs
        color = :white,
        kwargs...
    )
    _render_mesh!(axis, css; color, kwargs...)
    return nothing
end

function render!(
        axis::_RenderEnv,
        css::BMO.ConvexSphericalSurfaceSDF;
        # Makie kwargs
        color = :white,
        kwargs...
    )
    _render_mesh!(axis, css; color, kwargs...)
    return nothing
end

function render!(
        axis::_RenderEnv,
        asp::BMO.AbstractAsphericalSurfaceSDF;
        color = :white,
        kwargs...
    )
    _render_mesh!(axis, asp; color, kwargs...)
    return nothing
end
