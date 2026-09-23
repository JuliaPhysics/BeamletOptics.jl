function _conic_from_conjugates(s, s′)
    iszero(s + s′) && throw(ArgumentError("s + s′ must be non-zero (s = $s, s′ = $s′)"))
    R = 2 * s * s′ / (s + s′)
    k = -((s′ - s) / (s′ + s))^2
    return R, k
end

function _conic_y_range(R, k, x_off, diameter)
    # extrema of y_surf over the aperture disc; Z is monotone in r
    T = promote_type(typeof(R), typeof(k), typeof(x_off), typeof(diameter))
    r_hi = abs(x_off) + diameter / 2
    r_lo = max(zero(T), abs(x_off) - diameter / 2)
    Z_off = _conic_sag(x_off, R, k)
    y_a = -(_conic_sag(r_lo, R, k) - Z_off)
    y_b = -(_conic_sag(r_hi, R, k) - Z_off)
    return y_a <= y_b ? (y_a, y_b) : (y_b, y_a)
end

function _conic_auto_thickness(R, k, x_off, diameter, ::Type{T}) where {T}
    y_min, y_max = _conic_y_range(R, k, x_off, diameter)
    return y_max - y_min + T(10e-3)
end

"""
    OffAxisConicMirror(R, k, x_off, diameter; thickness=nothing, hole_diameter=nothing)

Constructs an off-axis segment of a general conic-of-revolution [`Mirror`](@ref) (sphere,
paraboloid, ellipsoid or hyperboloid), offset by `x_off` from the parent vertex. See
[`ConicSDF`](@ref) for the frame convention (origin, parent axis, opening direction, vertex
location).

# Inputs

- `R`:          Radius of curvature at the parent vertex \\[m\\]; `R > 0` is concave (opens
                towards `-y`), `R < 0` is convex (opens towards `+y`). Must be non-zero.
- `k`:          Conic constant. `k = -1` is a paraboloid, `k = 0` a sphere, `-1 < k <= 0` a
                prolate ellipsoid, `k > 0` an oblate ellipsoid, `k < -1` a hyperboloid.
- `x_off`:      Off-axis distance from the parent vertex to the aperture center \\[m\\]
- `diameter`:   Mirror aperture diameter \\[m\\]
- `thickness`:  Substrate thickness \\[m\\], calculated automatically to ensure solid backing
                if `nothing` (default)
- `hole_diameter`:  Diameter of the central through-hole \\[m\\], no hole if `nothing` (default). Must satisfy `0 < hole_diameter < diameter`. The bore is parallel to the local +y axis through the aperture centre.

For `k > -1` the aperture must stay within the domain of the parent conic
(`abs(x_off) + diameter/2 < abs(R)/sqrt(1+k)`), otherwise an `ArgumentError` is thrown; for
`k <= -1` there is no such limit. See also [`ConicMirror`](@ref) for the on-axis case.

Note that `hole_axis` (collimated/focused bore) exists only for [`OffAxisParabolicMirror`](@ref).
"""
function OffAxisConicMirror(
        R::Real,
        k::Real,
        x_off::Real,
        diameter::Real;
        thickness::Union{Real, Nothing} = nothing,
        hole_diameter::Union{Real, Nothing} = nothing
    )
    T = float(promote_type(typeof(R), typeof(k), typeof(x_off), typeof(diameter),
        typeof(thickness === nothing ? 0.0 : thickness)))
    Rt, kt, x_offt, dt = T(R), T(k), T(x_off), T(diameter)
    if thickness === nothing
        # Validate (R, k, x_off, diameter) via the SDF constructor itself first: computing
        # the auto thickness below evaluates the sag function outside its domain if the
        # aperture is invalid, which would raise a DomainError instead of an ArgumentError.
        ConicSDF(Rt, kt, x_offt, dt, one(T))
        t = _conic_auto_thickness(Rt, kt, x_offt, dt, T)
    else
        t = T(thickness)
    end
    conic_sdf = ConicSDF(Rt, kt, x_offt, dt, t)
    if hole_diameter === nothing
        return Mirror(conic_sdf)
    end
    (y_min, _) = _conic_y_range(Rt, kt, x_offt, dt)
    return Mirror(_pierce_substrate(conic_sdf, dt, min(y_min, zero(T)), t, hole_diameter))
end

"""
    ConicMirror(R, k, diameter; thickness=nothing, hole_diameter=nothing)

Constructs an on-axis segment of a general conic-of-revolution [`Mirror`](@ref) (sphere,
paraboloid, ellipsoid or hyperboloid). The vertex lies at the origin and the mirror opens
towards the negative y-axis for `R > 0`. See [`OffAxisConicMirror`](@ref) for the off-axis
case and the full sign/domain conventions.

If `hole_diameter` is given, a cylindrical bore centred on the optical axis (the local
+y-axis) is subtracted from the substrate, piercing it completely (e.g. Cassegrain,
Ritchey-Chrétien, or Dall-Kirkham primary).

# Inputs

- `R`:              Radius of curvature at the vertex \\[m\\]; `R > 0` concave, `R < 0` convex
- `k`:              Conic constant; `k = -1` is a paraboloid (see [`ParabolicMirror`](@ref)),
                    `k = 0` a sphere (see [`SphericalMirror`](@ref))
- `diameter`:       Mirror aperture diameter \\[m\\]
- `thickness`:      Substrate thickness \\[m\\], calculated automatically to ensure solid backing
                    if `nothing` (default)
- `hole_diameter`:  Diameter of the central through-hole \\[m\\], no hole if `nothing` (default). Must satisfy `0 < hole_diameter < diameter`.
"""
function ConicMirror(
        R::Real,
        k::Real,
        diameter::Real;
        thickness::Union{Real, Nothing} = nothing,
        hole_diameter::Union{Real, Nothing} = nothing
    )
    return OffAxisConicMirror(R, k, 0, diameter; thickness, hole_diameter)
end
