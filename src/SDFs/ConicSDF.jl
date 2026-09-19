_conic_sag(r, R, k)   = r^2 / (R * (1 + sqrt(1 - (1 + k) * r^2 / R^2)))
_conic_slope(r, R, k) = (r / R) / sqrt(1 - (1 + k) * r^2 / R^2)

"""
    ConicSDF{T} <: AbstractSDF{T}

Signed distance function representation of a segment of a conic of revolution
(sphere, paraboloid, prolate/oblate ellipsoid or hyperboloid), used as the substrate of
[`ParabolicMirror`](@ref), [`ConicMirror`](@ref), [`EllipsoidalMirror`](@ref),
[`HyperbolicMirror`](@ref) and their `OffAxis*` twins.

# Frame convention

- The origin is the segment centre, lying **on** the reflecting surface.
- The parent optical axis is parallel to the local y-axis, offset by `x_off` along the
  negative x-axis; there is no tilt between the segment and its parent axis.
- A concave surface (`R > 0`) opens towards the negative y-axis; a convex surface
  (`R < 0`) opens towards the positive y-axis.
- The parent vertex lies at `(-x_off, +Z(x_off), 0)`, where `Z` is the sag function of the
  parent conic.
- For `x_off = 0` and `k = -1`, `R = 2f`: the vertex lies at the origin and the (real)
  focus lies at `(0, -f, 0)`. This on-axis case needs no separate type: `x_off = 0` is an
  ordinary argument value, not a different `ConicSDF` variant.

# Fields

- `f`: paraxial focal length of the parent conic, `R/2` \\[m\\]
- `k`: conic constant (`k = -1` paraboloid, `k = 0` sphere, `-1 < k < 0` prolate ellipsoid,
  `k > 0` oblate ellipsoid, `k < -1` hyperboloid)
- `x_off`: off-axis distance from the parent vertex to the segment centre \\[m\\]
- `diameter`: segment aperture diameter \\[m\\]
- `thickness`: substrate thickness in y-direction \\[m\\]
- `pos`: position point in world coordinates \\[m\\]
- `dir`: rotation matrix in world coordinates
- `transposed_dir`: transposed rotation matrix

# Surface equation

With `r_p = sqrt((x + x_off)^2 + z^2)` the parent radius and

```math
Z(r) = \\frac{r^2}{R\\left(1 + \\sqrt{1 - (1+k)\\,r^2/R^2}\\right)}
```

the reflecting surface is `y_surf(x, z) = -(Z(r_p) - Z(x_off))`. The substrate is the solid
`{y_surf <= y <= thickness, sqrt(x^2 + z^2) <= diameter/2}`.
"""
mutable struct ConicSDF{T} <: AbstractSDF{T}
    f::T            # paraxial focal length, R/2
    k::T            # conic constant
    x_off::T
    diameter::T
    thickness::T
    pos::Point3{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
end

"""
    ConicSDF(R, k, x_off, diameter, thickness)

Constructs a [`ConicSDF`](@ref) representing a segment of a conic of revolution.

# Inputs

- `R`: radius of curvature at the parent vertex \\[m\\]; `R > 0` is concave (opens towards
  `-y`), `R < 0` is convex (opens towards `+y`). Must be non-zero.
- `k`: conic constant. `k = -1` is a paraboloid, `k = 0` a sphere, `-1 < k <= 0` a prolate
  ellipsoid, `k > 0` an oblate ellipsoid, `k < -1` a hyperboloid.
- `x_off`: off-axis distance from the parent vertex to the segment centre \\[m\\]
- `diameter`: segment aperture diameter \\[m\\]. Must be positive.
- `thickness`: substrate thickness in y-direction \\[m\\]

For `k > -1` the parent conic is only defined for `r < abs(R)/sqrt(1+k)`; the constructor
throws an `ArgumentError` if the aperture (`abs(x_off) + diameter/2`) reaches or exceeds this
limit. For `k <= -1` the surface is defined for all `r`, so no such limit applies.
"""
function ConicSDF(R::RR, k::K, x_off::X, diameter::D,
        thickness::Th) where {RR<:Real, K<:Real, X<:Real, D<:Real, Th<:Real}
    T = float(promote_type(RR, K, X, D, Th))
    Rt, kt, x_offt, dt, tt = T(R), T(k), T(x_off), T(diameter), T(thickness)

    iszero(Rt) && throw(ArgumentError("R must be non-zero; use a plano mirror for a flat surface"))
    dt <= 0 && throw(ArgumentError("diameter must be positive"))

    if kt > -1
        r_hi = abs(x_offt) + dt / 2
        r_lim = abs(Rt) / sqrt(1 + kt)
        if r_hi >= r_lim
            throw(ArgumentError(
                "aperture reaches r = $(r_hi) m but the conic (R = $R m, k = $k) is only " *
                "defined for r < $(r_lim) m; reduce diameter or x_off, or use k <= -1"))
        end
    end

    return ConicSDF{T}(
        Rt / 2, kt, x_offt, dt, tt,
        Point3{T}(0),
        SMatrix{3, 3, T, 9}(I),
        SMatrix{3, 3, T, 9}(I)
    )
end

"""
    radius(s::ConicSDF)

Returns the radius of curvature `R = 2f` of the parent conic at its vertex \\[m\\].
"""
radius(s::ConicSDF) = 2s.f

function bounding_sphere(s::ConicSDF{T}) where {T}
    R = 2s.f
    Z_off = _conic_sag(s.x_off, R, s.k)
    r_max = s.diameter / 2

    r_lo = max(zero(T), abs(s.x_off) - r_max)
    r_hi = abs(s.x_off) + r_max
    y_lo = -(_conic_sag(r_hi, R, s.k) - Z_off)
    y_hi = -(_conic_sag(r_lo, R, s.k) - Z_off)

    y_min = min(y_lo, y_hi, zero(T))
    y_max = max(y_lo, y_hi, s.thickness)

    y_center = (y_min + y_max) / 2
    r_bound = sqrt(r_max^2 + ((y_max - y_min) / 2)^2) + T(0.05)
    return Point3{T}(0, y_center, 0), r_bound
end

function sdf(s::ConicSDF{T}, point) where {T}
    p = _world_to_sdf(s, point)
    x, y, z = p[1], p[2], p[3]

    R = 2s.f
    Z_off = _conic_sag(s.x_off, R, s.k)

    r_hi = abs(s.x_off) + s.diameter / 2
    r_p = sqrt((x + s.x_off)^2 + z^2)
    r_c = min(r_p, r_hi)

    y_surf = -(_conic_sag(r_c, R, s.k) - Z_off)
    # `g` uses the *global* (aperture-rim) slope rather than the local one at `r_c`. `Z'`
    # is monotone increasing in `r` for every conic family here, so `r_hi` is where the
    # true gradient magnitude is largest over the whole segment; dividing by this
    # upper bound can only shrink `|d_front|` relative to the locally-linearized estimate
    # `(y_surf - y) / sqrt(1 + Z'(r_c)^2)`, i.e. it can only make `d_front` *more*
    # conservative. This matters because the local estimate is only a first-order
    # (graph/gradient) approximation of the true distance to a curved surface and can
    # slightly *overestimate* it away from the vertex (unlike a plane, a conic's true
    # distance-to-surface has a curvature-dependent second-order correction the local
    # estimate ignores) — this is exactly the failure mode the global-Lipschitz fallback
    # in the plan's Rationale anticipates. It does not change the surface itself: `d_front`
    # is still exactly zero on `y = y_surf(x, z)` regardless of `g`.
    g = sqrt(one(T) + _conic_slope(r_hi, R, s.k)^2)

    d_front = (y_surf - y) / g
    d_cyl = sqrt(x^2 + z^2) - s.diameter / 2
    d_back = y - s.thickness

    return max(d_front, d_cyl, d_back)
end

"""
    OffAxisParaboloidSDF(f, x_off, diameter, thickness)

Constructs a [`ConicSDF`](@ref) with `R = 2f` and `k = -1`, i.e. a segment of a paraboloid of
revolution. Kept for backwards compatibility; prefer `ConicSDF`.
"""
OffAxisParaboloidSDF(f::Real, x_off::Real, diameter::Real, thickness::Real) =
    ConicSDF(2f, -1, x_off, diameter, thickness)
