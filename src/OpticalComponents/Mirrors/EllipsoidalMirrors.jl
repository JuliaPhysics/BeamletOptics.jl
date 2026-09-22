function _ellipsoid_from_conjugates(s, s′)
    R, k = _conic_from_conjugates(s, s′)
    if !(-1 < k <= 0)
        throw(ArgumentError(
            "s and s′ must have the same sign for an ellipsoid (got s = $s, s′ = $s′ ⇒ k = $k); " *
            "use (OffAxis)HyperbolicMirror for a virtual focus or ParabolicMirror for s′ → ∞"))
    end
    return R, k
end

"""
    OffAxisEllipsoidalMirror(s, s′, x_off, diameter; thickness=nothing, hole_diameter=nothing)

Constructs an off-axis segment of an ellipsoidal [`Mirror`](@ref) whose two real conjugate
foci lie at object/image distances `s`, `s′` from the parent vertex, offset by `x_off`. Both
foci lie on the same side (the `-y`, i.e. reflecting, side) of the parent vertex for positive
`s`, `s′`.

# Inputs

- `s`, `s′`:    Conjugate object/image distances from the parent vertex \\[m\\], measured
                positive towards `-y` (in front of the mirror). Must have the same sign.
- `x_off`:      Off-axis distance from the parent vertex to the aperture center \\[m\\]
- `diameter`:   Mirror aperture diameter \\[m\\]
- `thickness`:  Substrate thickness \\[m\\], calculated automatically to ensure solid backing
                if `nothing` (default)
- `hole_diameter`:  Diameter of the central through-hole \\[m\\], no hole if `nothing` (default). Must satisfy `0 < hole_diameter < diameter`. The bore is parallel to the local +y axis through the aperture centre.

The vertex radius of curvature and conic constant are derived via
`R = 2ss′/(s+s′)`, `k = -((s′-s)/(s′+s))^2`. See also [`EllipsoidalMirror`](@ref) for the
on-axis case.

Note that `hole_axis` (collimated/focused bore) exists only for [`OffAxisParabolicMirror`](@ref).
"""
function OffAxisEllipsoidalMirror(
        s::Real,
        s′::Real,
        x_off::Real,
        diameter::Real;
        thickness::Union{Real, Nothing} = nothing,
        hole_diameter::Union{Real, Nothing} = nothing
    )
    R, k = _ellipsoid_from_conjugates(s, s′)
    return OffAxisConicMirror(R, k, x_off, diameter; thickness, hole_diameter)
end

"""
    EllipsoidalMirror(s, s′, diameter; thickness=nothing, hole_diameter=nothing)

Constructs an on-axis segment of an ellipsoidal [`Mirror`](@ref) whose two real conjugate
foci lie at `(0, -s, 0)` and `(0, -s′, 0)`; the vertex lies at the origin. See
[`OffAxisEllipsoidalMirror`](@ref) for the off-axis case and the sign convention for `s`, `s′`.

If `hole_diameter` is given, a cylindrical bore centred on the optical axis (the local
+y-axis) is subtracted from the substrate, piercing it completely (e.g. Dall-Kirkham primary).

# Inputs

- `s`, `s′`:        Conjugate object/image distances from the vertex \\[m\\], same sign
- `diameter`:       Mirror aperture diameter \\[m\\]
- `thickness`:      Substrate thickness \\[m\\], calculated automatically to ensure solid backing
                    if `nothing` (default)
- `hole_diameter`:  Diameter of the central through-hole \\[m\\], no hole if `nothing` (default). Must satisfy `0 < hole_diameter < diameter`.
"""
function EllipsoidalMirror(
        s::Real,
        s′::Real,
        diameter::Real;
        thickness::Union{Real, Nothing} = nothing,
        hole_diameter::Union{Real, Nothing} = nothing
    )
    return OffAxisEllipsoidalMirror(s, s′, 0, diameter; thickness, hole_diameter)
end
