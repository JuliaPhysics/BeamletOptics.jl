"""
    AbstractSurface{T}

A generic type for a surface which is basically an information storage type in order to build
shapes (volumes) from a combination of surfaces.

# Implementation reqs.

ToDo

## Getters/setters

ToDo

## Functions:

ToDo
"""
abstract type AbstractSurface{T} end

"""
    AbstractRotationallySymmetricSurface{T} <: AbstractSurface{T}

A surface type which is rotationally symmetric around one axis.
"""
abstract type AbstractRotationallySymmetricSurface{T} <: AbstractSurface{T} end

"""
    radius(s::AbstractRotationallySymmetricSurface)

Returns the radius of curvature of the surface. This might return `Inf` for planar surfaces
or surfaces which cannot be described by just one curvature radius.
"""
radius(s::AbstractRotationallySymmetricSurface) = s.radius

"""
    diameter(s::AbstractRotationallySymmetricSurface)

Returns the clear optical diameter of the surface.
"""
diameter(s::AbstractRotationallySymmetricSurface) = s.diameter

"""
    mechanical_diameter(s::AbstractRotationallySymmetricSurface)

Returns the mechanical diameter of the surface.

!!! note
It is assumed that mechanical_diameter(s) >= diameter(s) always holds.
"""
mechanical_diameter(s::AbstractRotationallySymmetricSurface) = diameter(s)

"""
Returns the sagitta of the surface at it edge, i.e. at `diameter(s)`

"""
edge_sag(s::AbstractRotationallySymmetricSurface, ::AbstractSDF) = s.edge_sag

abstract type AbstractOrientationType end

struct ForwardOrientation <: AbstractOrientationType end
struct ForwardLeftMeniscusOrientation <: AbstractOrientationType end
struct ForwardRightMeniscusOrientation <: AbstractOrientationType end

struct BackwardOrientation <: AbstractOrientationType end
struct BackwardLeftMeniscusOrientation <: AbstractOrientationType end
struct BackwardRightMeniscusOrientation <: AbstractOrientationType end

sdf(s::AbstractRotationallySymmetricSurface, ::Union{Nothing, AbstractOrientationType}) = throw(ArgumentError(lazy"sdf of $(typeof(s)) not implemented"))
sdf(s::AbstractRotationallySymmetricSurface) = sdf(s, nothing)
