"""
    AbstractSurface{T}

A generic type for a surface which is basically an information storage type in order to build
shapes (volumes) from a combination of surfaces.
"""
abstract type AbstractSurface{T} end

"""
    AbstractRotationallySymmetricSurface{T} <: AbstractSurface{T}

A surface type which is rotationally symmetric around one axis.

# Implementation reqs.

Subtypes of `AbstractShape` should implement the following:

## Getters/setters

- [`radius`](@ref) : Returns the radius of curvature of the `AbstractRotationallySymmetricSurface`
- [`diameter`](@ref) : Returns the clear optical diameter of the `AbstractRotationallySymmetricSurface`
- [`mechanical_diameter`](@ref) : Returns the mechanical diameter of the `AbstractRotationallySymmetricSurface`
- [`edge_sag`](@ref) : Returns the edge sagitta of the `AbstractRotationallySymmetricSurface`

## Functions:
- [`sdf(::AbstractRotationallySymmetricSurface, ::Union{Nothing, AbstractOrientationType})`](@ref) :
    Converts the surface specification of `AbstractRotationallySymmetricSurface` into an `AbstractSDF`
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

"""
    sdf(::AbstractRotationallySymmetricSurface, ::Union{Nothing, AbstractOrientationType})

Takes the surface specification and an optional `AbstractOrientationType` as
trait parameter and returns a corresponding `AbstractSDF` type.

# Surface vs. volume based tracing

This function is a mere convenience provider for users coming from other optic simulations frameworks which
are surface oriented. The goal of this function is to return the best matching closed volume SDF
which posesses a surface with the given specs on one side and most often a boundary and planar surface on the other side.

!!! info
    Always keep in mind that this package performs closed-volume baced ray tracing using either SDFs or meshes.
"""
sdf(s::AbstractRotationallySymmetricSurface, ::Union{Nothing, AbstractOrientationType}) = throw(ArgumentError(lazy"sdf of $(typeof(s)) not implemented"))
sdf(s::AbstractRotationallySymmetricSurface) = sdf(s, nothing)


"""
    CircularFlatSurface{T} <: AbstractRotationallySurface{T}

A type representing a planar circular surface, which is only parametrized by its `diameter`.

# Fields
- `diameter::T`: The diameter of the planar surface

"""
struct CircularFlatSurface{T} <: AbstractRotationallySymmetricSurface{T}
    diameter::T
end

sdf(::CircularFlatSurface, ::AbstractOrientationType) = nothing
