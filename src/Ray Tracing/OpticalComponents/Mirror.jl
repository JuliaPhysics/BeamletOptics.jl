"""
    AbstractReflectiveOptic <: AbstractObject

A generic type to represent `AbstractObject`s which reflect incoming rays. The main function of `interact3d` should be akin to [`reflection3d`](@ref).
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
        ::AbstractReflectiveOptic,
        ::Beam{T, R},
        ray::R) where {T <: Real, R <: PolarizedRay{T}}
    normal = normal3d(intersection(ray))
    npos = position(ray) + length(ray) * direction(ray)
    ndir = reflection3d(direction(ray), normal)
    # Jones reflection matrix
    J = @SArray [-1 0 0; 0 1 0; 0 0 1]
    E0 = _calculate_global_E0(direction(ray), ndir, J, polarization(ray))
    return BeamInteraction{T, R}(nothing,
        PolarizedRay{T}(
            npos, ndir, nothing, wavelength(ray), refractive_index(ray), E0))
end

"""
    Mirror{S <: AbstractShape} <: AbstractReflectiveOptic

Concrete implementation of a perfect mirror with arbitrary shape.
"""
struct Mirror{T, S <: AbstractShape{T}} <: AbstractReflectiveOptic{T, S}
    shape::S
end

"""
    RoundPlanoMirror

An ideal cylindrical mirror with planar reflecting surface, e.g. R = 1.

# Fields

- `shape`: a [`PlanoSurfaceSDF`](@ref) that represents the substrate

!!! warning "Reflecting surfaces"
    It is important to consider that **all** surfaces of this mirror type are reflecting!
"""
struct RoundPlanoMirror{T} <: AbstractReflectiveOptic{T, PlanoSurfaceSDF{T}}
    shape::PlanoSurfaceSDF{T}
end

"""
    RoundPlanoMirror(diameter, thickness)

Returns a cylindrical, flat [`Mirror`](@ref) with perfect reflectivity based on:

- `diameter`: mirror diameter in [m]
- `thickness`: mirror substrate thickness in [m]
"""
function RoundPlanoMirror(diameter::D, thickness::T) where {D<:Real,T<:Real}
    shape = PlanoSurfaceSDF(thickness, diameter)
    return RoundPlanoMirror(shape)
end