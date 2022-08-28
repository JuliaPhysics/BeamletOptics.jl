"""
    AbstractSurface{T}

Abstract supertype for all surfaces which might be used for a lens or even a mirror.

The standard interface for a surface contains:
- sag(::AbstractSurface{T}, r::Real)
- mechanical_semi_diameter(::AbstractSurface{T})
- chip_zone(::AbstractSurface{T})
- clear_semi_diameter(s::AbstractSurface)
"""
abstract type AbstractSurface{T} end

"""
    sag(s::AbstractSurface, r::Real)

Calculates the sag(gita) of the given surface `s` at position `r`.
"""
sag(::AbstractSurface{T}, r::Real) where T = zero(T)

"""
    mechanical_semi_diameter(::AbstractSurface{T})

Returns the mechanical semi-diameter of the surface. This means the full mechanical diameter.
"""
mechanical_semi_diameter(::AbstractSurface{T}) where T = zero(T)

"""
    chip_zone(::AbstractSurface{T})

Returns the chip zone of the surface (a flat bevel at the edge of the surface without curvature)
"""
chip_zone(::AbstractSurface{T}) where T = zero(T)

"""
    clear_semi_diameter(s::AbstractSurface)

Returns the clear aperture of the surface which is the difference between the mechanical
semi-diameter and the chip zone.
"""
clear_semi_diameter(s::AbstractSurface) = mechanical_semi_diameter(s) - chip_zone(s)

"""
    PlanarSurface{T}

Represents a planar surface without any curvature but with a physical size (mechanical semi-diameter).
"""
struct PlanarSurface{T} <: AbstractSurface{T}
    r::T
end

mechanical_semi_diameter(s::PlanarSurface) = s.r

"""
    ConicSurface{T} <: AbstractSurface{T}

Abstract super-type for all conic surfaces. Extends the abstract surface interface by the methods:
- conic_constant(::ConicSurface)
- radius_of_curvature(::ConicSurface{T})
"""
abstract type ConicSurface{T} <: AbstractSurface{T} end

conic_constant(::ConicSurface) = error("Not implemented")
radius_of_curvature(::ConicSurface{T}) where T = typemax(T)

"""
    SphericSurface{T} <: ConicSurface{T}

Concrete type representing a spherical surface with a given curvature.
"""
mutable struct SphericSurface{T} <: ConicSurface{T}
    R::T
    mechanical_semi_diameter::T
    chip_zone::T
end
conic_constant(::SphericSurface{T}) where T = zero(T)
radius_of_curvature(s::SphericSurface) = s.R
mechanical_semi_diameter(s::SphericSurface) = s.mechanical_semi_diameter
chip_zone(s::SphericSurface) = s.chip_zone

sag(s::SphericSurface, r::Real) = (radius_of_curvature(s) - √(radius_of_curvature(s)^2 - r^2))

"""
    AsphericSurface{T} <: ConicSurface{T}

Concrete type representing an aspherical surface with a given curvature, conic constant
and aspheric coefficients. This type adresses "standard" aspheres according to ISO10110.
"""
mutable struct AsphericSurface{T} <: AbstractSurface{T}
    sphere::SphericSurface{T}
    K::T
    A_even::Vector{T}
    A_odd::Vector{T}
end
conic_constant(s::AsphericSurface) = s.K
# pass to aspheric surface
radius_of_curvature(s::AsphericSurface) = radius_of_curvature(s.sphere)
mechanical_semi_diameter(s::AsphericSurface) = mechanical_semi_diameter(s.sphere)
chip_zone(s::AsphericSurface) = chip_zone(s.sphere)

function sag(s::AsphericSurface, r::Real)
    return (r^2 / (s.R*(1+√(1-(1+s.K)*r^2/s.R^2)))
            + map(i->s.A_even[i-1]*s.r^(2i), 2:length(s.A_even))
            + map(i->s.A_odd[i]*s.r^(2i+1), 1:length(s.A_odd)))
end
