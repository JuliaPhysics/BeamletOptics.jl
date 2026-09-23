function _hyperboloid_from_conjugates(s, s′)
    R, k = _conic_from_conjugates(s, s′)
    if !(k < -1)
        throw(ArgumentError(
            "s and s′ must have opposite signs for a hyperboloid (got s = $s, s′ = $s′ ⇒ k = $k); " *
            "use (OffAxis)EllipsoidalMirror"))
    end
    return R, k
end

"""
    OffAxisHyperbolicMirror(s, s′, x_off, diameter; thickness=nothing, hole_diameter=nothing)

Constructs an off-axis segment of a hyperboloidal [`Mirror`](@ref) whose conjugate foci lie at
object/image distances `s`, `s′` from the parent vertex, offset by `x_off`. Exactly one focus
is virtual, i.e. `s` and `s′` have opposite signs.

# Inputs

- `s`, `s′`:    Conjugate object/image distances from the parent vertex \\[m\\], measured
                positive towards `-y` (real focus, in front of the mirror) and negative
                towards `+y` (virtual focus, behind the mirror). Must have opposite signs.
- `x_off`:      Off-axis distance from the parent vertex to the aperture center \\[m\\]
- `diameter`:   Mirror aperture diameter \\[m\\]
- `thickness`:  Substrate thickness \\[m\\], calculated automatically to ensure solid backing
                if `nothing` (default)
- `hole_diameter`:  Diameter of the central through-hole \\[m\\], no hole if `nothing` (default). Must satisfy `0 < hole_diameter < diameter`. The bore is parallel to the local +y axis through the aperture centre.

!!! note "Cassegrain secondary"
    The secondary of a Cassegrain/Gregory telescope sees the prime focus behind itself, so
    pass it as a negative `s`.

The vertex radius of curvature and conic constant are derived via
`R = 2ss′/(s+s′)`, `k = -((s′-s)/(s′+s))^2`. See also [`HyperbolicMirror`](@ref) for the
on-axis case.

Note that `hole_axis` (collimated/focused bore) exists only for [`OffAxisParabolicMirror`](@ref).
"""
function OffAxisHyperbolicMirror(
        s::Real,
        s′::Real,
        x_off::Real,
        diameter::Real;
        thickness::Union{Real, Nothing} = nothing,
        hole_diameter::Union{Real, Nothing} = nothing
    )
    R, k = _hyperboloid_from_conjugates(s, s′)
    return OffAxisConicMirror(R, k, x_off, diameter; thickness, hole_diameter)
end

"""
    HyperbolicMirror(s, s′, diameter; thickness=nothing, hole_diameter=nothing)

Constructs an on-axis segment of a hyperboloidal [`Mirror`](@ref) (e.g. a Cassegrain/Gregory
secondary, or a Ritchey-Chrétien primary) whose conjugate foci lie at `(0, -s, 0)` and `(0, -s′, 0)`; the vertex lies at the
origin. See [`OffAxisHyperbolicMirror`](@ref) for the off-axis case and the sign convention.

If `hole_diameter` is given, a cylindrical bore centred on the optical axis (the local
+y-axis) is subtracted from the substrate, piercing it completely (e.g. Ritchey-Chrétien primary).

!!! note "Cassegrain secondary"
    The secondary of a Cassegrain/Gregory telescope sees the prime focus behind itself, so
    pass it as a negative `s`.

# Inputs

- `s`, `s′`:        Conjugate object/image distances from the vertex \\[m\\], opposite signs
- `diameter`:       Mirror aperture diameter \\[m\\]
- `thickness`:      Substrate thickness \\[m\\], calculated automatically to ensure solid backing
                    if `nothing` (default)
- `hole_diameter`:  Diameter of the central through-hole \\[m\\], no hole if `nothing` (default). Must satisfy `0 < hole_diameter < diameter`.
"""
function HyperbolicMirror(
        s::Real,
        s′::Real,
        diameter::Real;
        thickness::Union{Real, Nothing} = nothing,
        hole_diameter::Union{Real, Nothing} = nothing
    )
    return OffAxisHyperbolicMirror(s, s′, 0, diameter; thickness, hole_diameter)
end
