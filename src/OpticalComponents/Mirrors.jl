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

by specialized subtypes.

!!! info "Coating"
    Reflective optics now support an `AbstractCoating`. The default is `Uncoated`, which assumes perfect reflection for mirrors unless specified otherwise in the interaction logic.
    However, for `Mirror`, it is recommended to use `SimpleCoating` or `MultilayerCoating` if realistic losses are needed.


!!! info "Polarization ray tracing"
    Fresnel coefficients during reflection are set such that no reflection losses occur (i.e. `|rₚ| = |rₛ| = 1`).
"""
abstract type AbstractReflectiveOptic{T, C} <: AbstractObject{T} end

coating(object::AbstractReflectiveOptic) = object.coating

# FIXME Require reflectivity field/function for interaction with PolarizedRay

"""
    interact3d(AbstractReflectiveOptic, Ray)

Implements the reflection of a [`Ray`](@ref) via the normal at the intersection point on an optical surface.
"""
function interact3d(::AbstractSystem,
        obj::AbstractReflectiveOptic,
        ::Beam{T, R},
        ray::R) where {T <: Real, R <: Ray{T}}
    normal = normal3d(intersection(ray))
    npos = position(ray) + length(ray) * direction(ray)
    ndir = reflection3d(direction(ray), normal)

    # TODO: Handle coating wavelength dependence if we want to track intensity/power in Ray{T} in the future.
    # Currently Ray{T} is geometric only.

    return BeamInteraction{T, R}(nothing,
        Ray{T}(npos, ndir, nothing, wavelength(ray), refractive_index(ray)))
end

"""
    interact3d(AbstractReflectiveOptic, PolarizedRay)

Implements the reflection of a [`PolarizedRay`](@ref).
Calculates Fresnel/Coating coefficients.
"""
function interact3d(system::AbstractSystem,
        obj::AbstractReflectiveOptic,
        ::Beam{T, R},
        ray::R) where {T <: Real, R <: PolarizedRay{T}}
    normal = normal3d(intersection(ray))
    npos = position(ray) + length(ray) * direction(ray)

    # Calculate angle of incidence
    θi = angle3d(direction(ray), -normal)

    # Get refractive indices
    # For a mirror, "transmission" goes into the substrate or backing.
    # We assume the ray is in n1.
    n1 = refractive_index(ray)
    # What is n2? For a mirror in air, n1=1.
    # If it's a solid substrate mirror, n2 might be the glass.
    # But often for a mirror we just care about the coating on top.
    # Let's assume n2 = n1 if not specified, or we need a way to get substrate index.
    # But AbstractReflectiveOptic doesn't mandate a refractive_index field.
    n2 = n1 # simplified assumption if unknown

    # Get coefficients from coating
    # Note: For efficient mirrors, we usually just want the reflection part.
    c_point = get_coating_at(coating(obj), obj, intersection(ray), normal)
    rs, rp, ts, tp = coating_coefficients(c_point, n1, n2, wavelength(ray), θi)

    ndir = reflection3d(direction(ray), normal)

    # Jones reflection matrix construction
    # We are in the s-p basis.
    # Standard Fresenl reflection vector J * Ein
    # J = [rs 0; 0 rp]  (if we follow the same basis conventions as interact3d for Lenses)
    # Note: Mirrors often defined with rp = -1 for perfect conductor at normal incidence?
    # Our coating_coefficients should handle valid r values.
    # The previous implementation assumed diagonal [-1, 1].

    J = SPBasis(rs, 0, 0, rp)

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
struct Mirror{T, S <: AbstractShape{T}, C <: AbstractCoating} <:
       AbstractReflectiveOptic{T, C}
    shape::S
    coating::C
end

Mirror(shape::S) where {S <: AbstractShape} = Mirror(shape, SimpleCoating(1.0))

"""
    SquarePlanoMirror2D(edge_length)

Constructs a 2D square plano [`Mirror`](@ref) with a given `edge_length`.
The reflecting surface is normal to the y-axis.

# Inputs

- `edge_length`: the edge length of the square mirror in [m]
"""
function SquarePlanoMirror2D(
        size::T, coating::C = SimpleCoating(1.0)) where {T <: Real, C <: AbstractCoating}
    shape = QuadraticFlatMesh(size)
    return Mirror(shape, coating)
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
function RectangularPlanoMirror(width::W,
        height::H,
        thickness::T,
        coating::C = SimpleCoating(1.0)) where {
        W <: Real, H <: Real, T <: Real, C <: AbstractCoating}
    shape = CuboidMesh(width, thickness, height)
    translate3d!(shape, [
        -width / 2,       # x
        0,              # y
        -height / 2      # z
    ])
    set_new_origin3d!(shape)
    return Mirror(shape, coating)
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
function SquarePlanoMirror(width::W, thickness::T,
        coating::C = SimpleCoating(1.0)) where {W <: Real, T <: Real, C <: AbstractCoating}
    return RectangularPlanoMirror(width, width, thickness, coating)
end

"""
    RoundPlanoMirror <: AbstractReflectiveOptic

An ideal cylindrical mirror with planar reflecting surface, e.g. R = 1.
See also [`Mirror`](@ref).

# Fields

- `shape`: a [`PlanoSurfaceSDF`](@ref) that represents the substrate
"""
struct RoundPlanoMirror{T, C <: AbstractCoating} <: AbstractReflectiveOptic{T, C}
    shape::PlanoSurfaceSDF{T}
    coating::C
end

"""
    RoundPlanoMirror(diameter, thickness)

Returns a cylindrical, flat [`RoundPlanoMirror`](@ref) with perfect reflectivity based on:

# Inputs

- `diameter`: mirror diameter in [m]
- `thickness`: mirror substrate thickness in [m]
"""
function RoundPlanoMirror(diameter::D, thickness::T,
        coating::C = SimpleCoating(1.0)) where {D <: Real, T <: Real, C <: AbstractCoating}
    shape = PlanoSurfaceSDF(thickness, diameter)
    return RoundPlanoMirror(shape, coating)
end

"""[`ConcaveSphericalMirror`](@ref) shape type based on a [`UnionSDF`](@ref)"""
const ConcaveSphericalMirrorShape{T} = UnionSDF{
    T, Tuple{ConcaveSphericalSurfaceSDF{T}, PlanoSurfaceSDF{T}}}

"""
    ConcaveSphericalMirror <: AbstractReflectiveOptic

An ideal concave mirror with spherical reflecting surface, e.g. R = 1.
See also [`RoundPlanoMirror`](@ref).

# Fields

- `shape`: a [`ConcaveSphericalMirrorShape`](@ref) that represents the substrate
"""
struct ConcaveSphericalMirror{T, C <: AbstractCoating} <: AbstractReflectiveOptic{T, C}
    shape::ConcaveSphericalMirrorShape{T}
    coating::C
end

"""
    ConcaveSphericalMirror(radius, thickness, diameter)

Constructor for a spherical mirror with a concave reflecting surface. The component is aligned with the positive y-axis.
See also [`ConcaveSphericalMirror`](@ref).

# Inputs

- `radius`: the spherical surface radius of curvature in [m]
- `thickness`: substrate thickness in [m]
- `diameter`: mirror outer diameter in [m]
"""
function ConcaveSphericalMirror(radius::Real, thickness::Real, diameter::Real,
        coating::C = SimpleCoating(1.0)) where {C <: AbstractCoating}
    cylinder = PlanoSurfaceSDF(thickness, diameter)
    concave = ConcaveSphericalSurfaceSDF(abs(radius), diameter)
    shape = concave + cylinder
    return ConcaveSphericalMirror(shape, coating)
end

"""
    RightAnglePrismMirror <: AbstractReflectiveOptic

An ideal right angle prism mirror with planar reflecting surface, i.e. R = 1.
See also [`Mirror`](@ref).

# Fields

- `shape`: a [`RightAnglePrismSDF`](@ref) that represents the substrate
"""
struct RightAnglePrismMirror{T, C <: AbstractCoating} <: AbstractReflectiveOptic{T, C}
    shape::RightAnglePrismSDF{T}
    coating::C
end

"""
    RightAnglePrismMirror(leg_length, height)

Constructs a right angle prism mirror. The primary surface is aligned with the pos. y-axis.

# Inputs

- `leg_length`: edge length in x and y in [m]
- `height`: in z-axis in [m]
"""
function RightAnglePrismMirror(leg_length::Real, height::Real,
        coating::C = SimpleCoating(1.0)) where {C <: AbstractCoating}
    shape = RightAnglePrismSDF(leg_length, height)
    zrotate3d!(shape, deg2rad(45 + 180))
    return RightAnglePrismMirror(shape, coating)
end
