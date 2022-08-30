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
mechanical_semi_diameter(s::AbstractSurface{T}) where T = s.mechanical_semi_diameter

"""
    chip_zone(::AbstractSurface{T})

Returns the chip zone of the surface (a flat bevel at the edge of the surface without curvature)
"""
chip_zone(s::AbstractSurface{T}) where T = s.chip_zone

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
    mechanical_semi_diameter::T
end
chip_zone(s::PlanarSurface{T}) where T = zero(T)

"""
    ConicSurface{T} <: AbstractSurface{T}

Abstract super-type for all conic surfaces. Extends the abstract surface interface by the methods:
- conic_constant(::ConicSurface)
- radius_of_curvature(::ConicSurface{T})
"""
abstract type ConicSurface{T} <: AbstractSurface{T} end

conic_constant(s::ConicSurface) = s.K
radius_of_curvature(s::ConicSurface)= s.R

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

sag(s::SphericSurface, r::Real) = (radius_of_curvature(s) - √(radius_of_curvature(s)^2 - r^2))

"""
    AsphericSurface{T} <: ConicSurface{T}

Concrete type representing an aspherical surface with a given curvature, conic constant
and aspheric coefficients. This type adresses "standard" aspheres according to ISO10110.
"""
mutable struct AsphericSurface{T} <: ConicSurface{T}
    sphere::SphericSurface{T}
    K::T
    A_even::Vector{T}
    A_odd::Vector{T}
end

# pass to aspheric surface
radius_of_curvature(s::AsphericSurface) = radius_of_curvature(s.sphere)
mechanical_semi_diameter(s::AsphericSurface) = mechanical_semi_diameter(s.sphere)
chip_zone(s::AsphericSurface) = chip_zone(s.sphere)

function sag(s::AsphericSurface{T}, r::Real) where T
    return (r^2 / (radius_of_curvature(s)*(1+√(1-(1+conic_constant(s))*r^2/radius_of_curvature(s)^2)))
            + mapreduce(i->s.A_even[i-1]*r^(2i), +, 2:length(s.A_even), init=zero(T))
            + mapreduce(i->s.A_odd[i]*r^(2i+1), +, 1:length(s.A_odd), init=zero(T)))
end

mutable struct CylinderSurface{T} <: AbstractSurface{T}
    R::T
    mechanical_semi_diameter::T
    chip_zone::T
    d::Vector{T}
end
radius_of_curvature(s::CylinderSurface) = s.R
# Identical to sphere
sag(s::CylinderSurface, r::Real) = (radius_of_curvature(s) - √(radius_of_curvature(s)^2 - r^2))
