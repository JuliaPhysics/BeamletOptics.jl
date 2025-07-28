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
abstract type AbstractReflectiveOptic{T, S <: AbstractShape{T}} <: AbstractObject{T, S} end

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
    J = SPBasis([-1 0 0; 0 1 0; 0 0 1])
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
struct Mirror{T, S <: AbstractShape{T}} <: AbstractReflectiveOptic{T, S}
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

"""
    RoundPlanoMirror <: AbstractReflectiveOptic

An ideal cylindrical mirror with planar reflecting surface, e.g. R = 1.
See also [`Mirror`](@ref).

# Fields

- `shape`: a [`PlanoSurfaceSDF`](@ref) that represents the substrate
"""
struct RoundPlanoMirror{T} <: AbstractReflectiveOptic{T, PlanoSurfaceSDF{T}}
    shape::PlanoSurfaceSDF{T}
end

"""
    RoundPlanoMirror(diameter, thickness)

Returns a cylindrical, flat [`RoundPlanoMirror`](@ref) with perfect reflectivity based on:

# Inputs

- `diameter`: mirror diameter in [m]
- `thickness`: mirror substrate thickness in [m]
"""
function RoundPlanoMirror(diameter::D, thickness::T) where {D<:Real,T<:Real}
    shape = PlanoSurfaceSDF(thickness, diameter)
    return RoundPlanoMirror(shape)
end

"""[`ConcaveSphericalMirror`](@ref) shape type based on a [`UnionSDF`](@ref)"""
const ConcaveSphericalMirrorShape{T} = UnionSDF{T, Tuple{ConcaveSphericalSurfaceSDF{T}, PlanoSurfaceSDF{T}}}

"""
    ConcaveSphericalMirror <: AbstractReflectiveOptic

An ideal concave mirror with spherical reflecting surface, e.g. R = 1.
See also [`RoundPlanoMirror`](@ref).

# Fields

- `shape`: a [`ConcaveSphericalMirrorShape`](@ref) that represents the substrate
"""
struct ConcaveSphericalMirror{T} <: AbstractReflectiveOptic{T, ConcaveSphericalMirrorShape{T}}
    shape::ConcaveSphericalMirrorShape{T}
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
function ConcaveSphericalMirror(radius::Real, thickness::Real, diameter::Real)
    cylinder = PlanoSurfaceSDF(thickness, diameter)
    concave = ConcaveSphericalSurfaceSDF(abs(radius), diameter)
    shape = concave + cylinder
    return ConcaveSphericalMirror(shape)
end

"""
    RightAnglePrismMirror <: AbstractReflectiveOptic

An ideal right angle prism mirror with planar reflecting surface, i.e. R = 1.
See also [`Mirror`](@ref).

# Fields

- `shape`: a [`RightAnglePrismSDF`](@ref) that represents the substrate
"""
struct RightAnglePrismMirror{T} <: AbstractReflectiveOptic{T, RightAnglePrismSDF{T}}
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