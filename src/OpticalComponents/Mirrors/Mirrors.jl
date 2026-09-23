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

# order of inclusion matters
include("MirrorUtils.jl")
include("PlanoMirrors.jl")
include("SphericalMirrors.jl")
include("ParabolicMirrors.jl")
include("ConicMirrors.jl")
include("EllipsoidalMirrors.jl")
include("HyperbolicMirrors.jl")
