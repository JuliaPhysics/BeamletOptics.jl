"""
    AbstractReflectiveOptic <: AbstractObject

A generic type to represent an [`AbstractObject`] which reflects incoming rays.

# Implementation reqs.

Subtypes of `AbstractReflectiveOptic` should implement all supertype reqs. as well as:

## Fields

- no specific fields required

## Getters/setters

- none required

## Functions

- `interact3d`:  the interaction logic should be akin to [`reflection3d`](@ref) for each surface crossing

# Additional information

The information provided below applies to the standard functional implementation of this type and may be overwritten
by specialized subtypes.

!!! info "Polarization ray tracing"
    Fresnel coefficients during reflection are set such that no reflection losses occur (i.e. `|rₚ| = |rₛ| = 1`).
"""
abstract type AbstractReflectiveOptic{T} <: AbstractObject{T} end

# FIXME Require reflectivity field/function for interaction with PolarizedRay

"""
    interact3d(AbstractReflectiveOptic, Ray)

Implements the reflection of a [`Ray`](@ref) via the normal at the intersection point on an optical surface.
"""
function interact3d(::AbstractSystem,
        ::AbstractReflectiveOptic,
        ::Beam{T, R},
        ray::R) where {T <: Real, R <: Ray{T}}
    normal = normal3d(intersection(ray))
    npos = position(ray) + length(ray) * direction(ray)
    ndir = reflection3d(direction(ray), normal)
    return BeamInteraction{T, R}(nothing,
        Ray{T}(npos, ndir, nothing, wavelength(ray), refractive_index(ray)))
end

"""
    interact3d(AbstractReflectiveOptic, PolarizedRay)

Implements the ideal reflection of a [`PolarizedRay`](@ref) via the normal at the intersection point on an optical surface.
A Jones matrix of [-1 0 0; 0 1 0] is assumed as per Peatross (2015, 2023 Ed. p. 154) and Yun et al. (see [`PolarizedRay`](@ref) for more information).
"""
function interact3d(::AbstractSystem,
        obj::AbstractReflectiveOptic,
        ::Beam{T, R},
        ray::R) where {T <: Real, R <: PolarizedRay{T}}
    normal = normal3d(intersection(ray))
    npos = position(ray) + length(ray) * direction(ray)
    ndir = reflection3d(direction(ray), normal)
    # Jones reflection matrix
    J = SPBasis(-1, 0, 0, 1)
    E0 = _calculate_global_E0(obj, ray, ndir, J)
    return BeamInteraction{T, R}(nothing,
        PolarizedRay{T}(
            npos, ndir, nothing, wavelength(ray), refractive_index(ray), E0))
end

"""
    Mirror{S <: AbstractShape} <: AbstractReflectiveOptic

Concrete implementation of a perfect mirror (R = 1) with arbitrary shape.

!!! warning "Reflecting surfaces"
    It is important to consider that **all** surfaces of this mirror type are reflecting!
"""
struct Mirror{T, S <: AbstractShape{T}} <: AbstractReflectiveOptic{T}
    shape::S
end

"""
    SquarePlanoMirror2D(edge_length)

Constructs a 2D square plano [`Mirror`](@ref) with a given `edge_length`.
The reflecting surface is normal to the y-axis.

# Inputs

- `edge_length`: the edge length of the square mirror in [m]
"""
function SquarePlanoMirror2D(size::T) where {T <: Real}
    shape = QuadraticFlatMesh(size)
    return Mirror(shape)
end

"""
    RectangularPlanoMirror(width, height, thickness)

Constructs a rectangular plano [`Mirror`](@ref) based on the input dimensions.
The front reflecting surface is normal to the y-axis and lies at the origin.

# Inputs

- `width`:      of the mirror in x-direction [m] 
- `height`:     of the mirror in z-direction [m] 
- `thickness`:  of the mirror in y-direction [m] 
"""
function RectangularPlanoMirror(width::W, height::H, thickness::T) where {W<:Real,H<:Real,T<:Real}
    shape = CuboidMesh(width, thickness, height)
    translate3d!(shape, [
        -width/2,       # x
        0,              # y
        -height/2,      # z
    ])
    set_new_origin3d!(shape)
    return Mirror(shape)
end

"""
    SquarePlanoMirror(width, thickness)

Constructs a square plano [`Mirror`](@ref) with equal width and height.
The front reflecting surface is normal to the y-axis and lies at the origin.
See also [`RectangularPlanoMirror`](@ref).

# Inputs

- `width`: the side length of the square mirror in x- and y-direction [m]
- `thickness`: of the mirror in [m]
"""
function SquarePlanoMirror(width::W, thickness::T) where {W<:Real,T<:Real}
    return RectangularPlanoMirror(width, width, thickness)
end

function _pierce_substrate(substrate, diameter, t::T, sag_max::T, hole_diameter) where {T}
    hd = T(hole_diameter)
    if !(0 < hd < T(diameter))
        throw(ArgumentError("hole_diameter must satisfy 0 < hole_diameter < diameter (got $hole_diameter, diameter $diameter)"))
    end
    margin = T(10e-3)
    half_height = (t + sag_max) / 2 + margin
    y_center = (t - sag_max) / 2
    bore = CylinderSDF(hd / 2, half_height)
    translate3d!(bore, [zero(T), y_center, zero(T)])
    return Mirror(substrate - bore)
end

"""
    RoundPlanoMirror <: AbstractReflectiveOptic

An ideal cylindrical mirror with planar reflecting surface, e.g. R = 1.
See also [`Mirror`](@ref).

# Fields

- `shape`: a [`PlanoSurfaceSDF`](@ref) that represents the substrate
"""
struct RoundPlanoMirror{T} <: AbstractReflectiveOptic{T}
    shape::PlanoSurfaceSDF{T}
end

"""
    RoundPlanoMirror(diameter, thickness; hole_diameter=nothing)

Returns a cylindrical, flat [`RoundPlanoMirror`](@ref) (or pierced [`Mirror`](@ref)) with perfect reflectivity based on:

# Inputs

- `diameter`: mirror diameter in [m]
- `thickness`: mirror substrate thickness in [m]
- `hole_diameter`: diameter of the central through-hole in [m], no hole if `nothing` (default)
"""
function RoundPlanoMirror(diameter::D, thickness::T; hole_diameter::Union{Real, Nothing} = nothing) where {D<:Real,T<:Real}
    shape = PlanoSurfaceSDF(thickness, diameter)
    if hole_diameter === nothing
        return RoundPlanoMirror(shape)
    end
    T_res = typeof(float(thickness))
    return _pierce_substrate(shape, diameter, T_res(thickness), zero(T_res), hole_diameter)
end

"""[`SphericalMirror`](@ref) shape type based on a [`UnionSDF`](@ref)"""
const SphericalMirrorShape{T} = UnionSDF{T, Tuple{ConcaveSphericalSurfaceSDF{T}, PlanoSurfaceSDF{T}}}

"""
    SphericalMirror <: AbstractReflectiveOptic

An ideal concave mirror with spherical reflecting surface, e.g. R = 1.
See also [`RoundPlanoMirror`](@ref).

# Fields

- `shape`: a [`SphericalMirrorShape`](@ref) that represents the substrate
"""
struct SphericalMirror{T} <: AbstractReflectiveOptic{T}
    shape::SphericalMirrorShape{T}
end

"""
    SphericalMirror(radius, thickness, diameter; hole_diameter=nothing)

Constructor for a spherical mirror with a concave reflecting surface. The component is aligned with the positive y-axis.
See also [`SphericalMirror`](@ref).

# Inputs

- `radius`: the spherical surface radius of curvature in [m]
- `thickness`: substrate thickness in [m]
- `diameter`: mirror outer diameter in [m]
- `hole_diameter`: diameter of the central through-hole in [m], no hole if `nothing` (default)
"""
function SphericalMirror(radius::Real, thickness::Real, diameter::Real; hole_diameter::Union{Real, Nothing} = nothing)
    cylinder = PlanoSurfaceSDF(thickness, diameter)
    concave = ConcaveSphericalSurfaceSDF(abs(radius), diameter)
    shape = concave + cylinder
    if hole_diameter === nothing
        return SphericalMirror(shape)
    end
    T = float(promote_type(typeof(radius), typeof(thickness), typeof(diameter)))
    r_sph = T(abs(radius))
    d_half = T(diameter) / 2
    sag_max = r_sph >= d_half ? r_sph - sqrt(r_sph^2 - d_half^2) : d_half
    return _pierce_substrate(shape, diameter, T(thickness), sag_max, hole_diameter)
end

# Former name, kept for backwards compatibility
Base.@deprecate_binding ConcaveSphericalMirror SphericalMirror

"""
    RightAnglePrismMirror <: AbstractReflectiveOptic

An ideal right angle prism mirror with planar reflecting surface, i.e. R = 1.
See also [`Mirror`](@ref).

# Fields

- `shape`: a [`RightAnglePrismSDF`](@ref) that represents the substrate
"""
struct RightAnglePrismMirror{T} <: AbstractReflectiveOptic{T}
    shape::RightAnglePrismSDF{T}
end

"""
    RightAnglePrismMirror(leg_length, height)

Constructs a right angle prism mirror. The primary surface is aligned with the pos. y-axis.

# Inputs

- `leg_length`: edge length in x and y in [m] 
- `height`: in z-axis in [m]
"""
function RightAnglePrismMirror(leg_length::Real, height::Real)
    shape = RightAnglePrismSDF(leg_length, height)
    zrotate3d!(shape, deg2rad(45+180))
    return RightAnglePrismMirror(shape)
end

"""
    OffAxisParabolicMirror(rfl, diameter; angle=90, thickness=nothing, hole_diameter=nothing, hole_axis=:collimated)

Constructs an Off-Axis Parabolic (OAP) [`Mirror`](@ref) from:

# Inputs

- `rfl`:            Reflected Focal Length (distance from aperture center to focus) [m]
- `diameter`:       Mirror aperture diameter [m]
- `angle`:          Deflection angle in degrees (default: 90°)
- `thickness`:      Substrate thickness [m], calculated automatically to ensure solid backing if `nothing` (default)
- `hole_diameter`:  Diameter of the through-hole [m], no hole if `nothing` (default). Must satisfy `0 < hole_diameter < diameter`.
- `hole_axis`:      Orientation of the through-hole. Options:
                    - `:collimated` (default): parallel to the collimated beam (local y-axis / substrate normal), centered at the aperture center `(0, 0, 0)`.
                    - `:focused`: angled towards the parent paraboloid focus `(-x_off, -f, 0)`, passing through the aperture center `(0, 0, 0)` (e.g. for collinear pump-probe beams).
"""
function OffAxisParabolicMirror(
        rfl::Real,
        diameter::Real;
        angle::Real = 90,
        thickness::Union{Real, Nothing} = nothing,
        hole_diameter::Union{Real, Nothing} = nothing,
        hole_axis::Symbol = :collimated
    )
    T = float(promote_type(typeof(rfl), typeof(diameter), typeof(angle), typeof(thickness === nothing ? 0.0 : thickness)))
    angle_rad = deg2rad(angle)

    f = T(rfl * (cos(angle_rad / 2)^2))
    x_off = T(rfl * sin(angle_rad))

    r_max = T(diameter / 2)
    sag_max = abs(-(((r_max + x_off)^2 - x_off^2) / (4 * f)))

    t = thickness === nothing ? max(T(diameter / 2), sag_max + T(10e-3)) : T(thickness)

    oap_sdf = ConicSDF(2f, -one(T), x_off, T(diameter), t)
    if hole_diameter === nothing
        return Mirror(oap_sdf)
    end

    hd = T(hole_diameter)
    if !(0 < hd < T(diameter))
        throw(ArgumentError("hole_diameter must satisfy 0 < hole_diameter < diameter (got $hole_diameter, diameter $diameter)"))
    end

    if hole_axis === :collimated
        return _pierce_substrate(oap_sdf, diameter, t, sag_max, hole_diameter)
    elseif hole_axis === :focused
        focus_vec = Point3(-x_off, -T(rfl * cos(angle_rad)), zero(T))
        u_dir = normalize(focus_vec)
        margin = T(10e-3)
        span = max(T(diameter), t + sag_max) + 2margin
        bore = CylinderSDF(hd / 2, span)
        align3d!(bore, u_dir)
        return Mirror(oap_sdf - bore)
    else
        throw(ArgumentError("hole_axis must be :collimated or :focused, got :$hole_axis"))
    end
end

"""
    ParabolicMirror(f, diameter; thickness=nothing, hole_diameter=nothing)

Constructs an on-axis parabolic [`Mirror`](@ref) with focal length `f`.
The vertex of the concave reflecting surface lies at the origin, the mirror opens towards the negative y-axis
and its focus lies at `(0, -f, 0)`. The shape is an [`OffAxisParaboloidSDF`](@ref) without off-axis offset.

If `hole_diameter` is given, a cylindrical bore centred on the optical axis (the local
+y-axis) is subtracted from the substrate, piercing it completely (a Cassegrain primary).

!!! warning "Reflective bore wall"
    Rays that graze into the hole will reflect off its wall rather than being absorbed.

# Inputs

- `f`:              Focal length \\[m\\]
- `diameter`:       Mirror aperture diameter \\[m\\]
- `thickness`:      Substrate thickness \\[m\\], rim sag + 10 mm if `nothing` (default)
- `hole_diameter`:  Diameter of the central through-hole \\[m\\], no hole if `nothing` (default). Must satisfy `0 < hole_diameter < diameter`.
"""
function ParabolicMirror(
        f::Real,
        diameter::Real;
        thickness::Union{Real, Nothing} = nothing,
        hole_diameter::Union{Real, Nothing} = nothing
    )
    T = float(promote_type(typeof(f), typeof(diameter), typeof(thickness === nothing ? 0.0 : thickness)))
    sag_max = T(diameter / 2)^2 / (4 * T(f))
    t = thickness === nothing ? sag_max + T(10e-3) : T(thickness)
    substrate = ConicSDF(2 * T(f), -one(T), zero(T), T(diameter), t)
    if hole_diameter === nothing
        return Mirror(substrate)
    end
    return _pierce_substrate(substrate, diameter, t, sag_max, hole_diameter)
end

function _conic_from_conjugates(s, s′)
    iszero(s + s′) && throw(ArgumentError("s + s′ must be non-zero (s = $s, s′ = $s′)"))
    R = 2 * s * s′ / (s + s′)
    k = -((s′ - s) / (s′ + s))^2
    return R, k
end

function _conic_auto_thickness(R, k, x_off, diameter, ::Type{T}) where {T}
    # extrema of y_surf over the aperture disc; Z is monotone in r
    r_hi = abs(x_off) + diameter / 2
    r_lo = max(zero(T), abs(x_off) - diameter / 2)
    Z_off = _conic_sag(x_off, R, k)
    y_a = -(_conic_sag(r_lo, R, k) - Z_off)
    y_b = -(_conic_sag(r_hi, R, k) - Z_off)
    return abs(y_a - y_b) + T(10e-3)
end

"""
    OffAxisConicMirror(R, k, x_off, diameter; thickness=nothing)

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

For `k > -1` the aperture must stay within the domain of the parent conic
(`abs(x_off) + diameter/2 < abs(R)/sqrt(1+k)`), otherwise an `ArgumentError` is thrown; for
`k <= -1` there is no such limit. See also [`ConicMirror`](@ref) for the on-axis case.
"""
function OffAxisConicMirror(
        R::Real,
        k::Real,
        x_off::Real,
        diameter::Real;
        thickness::Union{Real, Nothing} = nothing
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
    return Mirror(ConicSDF(Rt, kt, x_offt, dt, t))
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
    m = OffAxisConicMirror(R, k, 0, diameter; thickness)
    if hole_diameter === nothing
        return m
    end
    t = shape(m).thickness
    T = typeof(t)
    sag_max = abs(_conic_sag(T(diameter / 2), T(R), T(k)))
    return _pierce_substrate(shape(m), diameter, t, sag_max, hole_diameter)
end

"""
    OffAxisEllipsoidalMirror(s, s′, x_off, diameter; thickness=nothing)

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

The vertex radius of curvature and conic constant are derived via
`R = 2ss′/(s+s′)`, `k = -((s′-s)/(s′+s))^2`. See also [`EllipsoidalMirror`](@ref) for the
on-axis case.
"""
function OffAxisEllipsoidalMirror(
        s::Real,
        s′::Real,
        x_off::Real,
        diameter::Real;
        thickness::Union{Real, Nothing} = nothing
    )
    R, k = _conic_from_conjugates(s, s′)
    if !(-1 < k <= 0)
        throw(ArgumentError(
            "s and s′ must have the same sign for an ellipsoid (got s = $s, s′ = $s′ ⇒ k = $k); " *
            "use OffAxisHyperbolicMirror for a virtual focus or ParabolicMirror for s′ → ∞"))
    end
    return OffAxisConicMirror(R, k, x_off, diameter; thickness)
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
    R, k = _conic_from_conjugates(s, s′)
    if !(-1 < k <= 0)
        throw(ArgumentError(
            "s and s′ must have the same sign for an ellipsoid (got s = $s, s′ = $s′ ⇒ k = $k); " *
            "use HyperbolicMirror for a virtual focus or ParabolicMirror for s′ → ∞"))
    end
    return ConicMirror(R, k, diameter; thickness, hole_diameter)
end

"""
    OffAxisHyperbolicMirror(s, s′, x_off, diameter; thickness=nothing)

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

!!! note "Cassegrain secondary"
    The secondary of a Cassegrain/Gregory telescope sees the prime focus behind itself, so
    pass it as a negative `s`.

The vertex radius of curvature and conic constant are derived via
`R = 2ss′/(s+s′)`, `k = -((s′-s)/(s′+s))^2`. See also [`HyperbolicMirror`](@ref) for the
on-axis case.
"""
function OffAxisHyperbolicMirror(
        s::Real,
        s′::Real,
        x_off::Real,
        diameter::Real;
        thickness::Union{Real, Nothing} = nothing
    )
    R, k = _conic_from_conjugates(s, s′)
    if !(k < -1)
        throw(ArgumentError(
            "s and s′ must have opposite signs for a hyperboloid (got s = $s, s′ = $s′ ⇒ k = $k); " *
            "use OffAxisEllipsoidalMirror"))
    end
    return OffAxisConicMirror(R, k, x_off, diameter; thickness)
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
    R, k = _conic_from_conjugates(s, s′)
    if !(k < -1)
        throw(ArgumentError(
            "s and s′ must have opposite signs for a hyperboloid (got s = $s, s′ = $s′ ⇒ k = $k); " *
            "use EllipsoidalMirror"))
    end
    return ConicMirror(R, k, diameter; thickness, hole_diameter)
end
