"""
    AbstractShape{T<:Real}

A generic type for a shape that exists in 3D-space. Must have a `pos`ition and orientation.
Types used to describe the geometry of a shape should be subtypes of `Real`.

# Implementation reqs.

Subtypes of `AbstractShape` should implement the following:

## Fields:

- `pos`: a 3D-vector that stores the current `position` of the object-specific coordinate system
- `dir`: a 3x3-matrix that represents the orthonormal basis of the object and therefore, the `orientation`

## Getters/setters

- [`position`](@ref) / `position!`: gets or sets the `pos`ition vector of the `AbstractShape`
- [`orientation`](@ref) / `orientation!`: gets or sets the orientation matrix of the `AbstractShape`

## Kinematic:

`AbstractShape`s are [`BeamletOptics.Movable`](@ref) with an [`BeamletOptics.Oriented`](@ref) frame,
see [`BeamletOptics.AbstractKinematicTrait`](@ref). The default primitives
`translate3d!(::Movable, shape, offset)` and `rotate3d!(::Movable, shape, R::AbstractMatrix)` act on
`position`/`orientation`; subtypes with additional geometry data (e.g. mesh vertices) dispatch their own.

## Ray Tracing:

- [`intersect3d`](@ref): returns the intersection between an `AbstractShape` and `AbstractRay`, or lack thereof. See also [`Intersection`](@ref)
- [`bounding_sphere`](@ref) (optional): a sphere that encloses the shape. The solver then skips the
  shape for every ray that misses the sphere, without a call of `intersect3d`

## Rendering (with Makie):

Refer to the [`render!`](@ref) documentation.
"""
abstract type AbstractShape{T <: Real} end

kinematic_trait_of(::AbstractShape) = Movable(Oriented())

"Enforces that `shape` has to have the field `pos` or implement `position()`."
Base.position(shape::AbstractShape) = shape.pos
position!(shape::AbstractShape, pos) = (shape.pos = pos)

"Enforces that `shape` has to have the field `dir` or implement `orientation()`."
orientation(shape::AbstractShape) = shape.dir
orientation!(shape::AbstractShape, dir) = (shape.dir = dir)

"""
    bounding_sphere(shape::AbstractShape)

Returns `nothing`, or the `center` and the `radius` in [m] of a sphere that encloses the `shape`,
as a tuple `(center, radius)`. The `center` is given in the local frame of the `shape`, i.e.
relative to its [`position`](@ref) and along the axes of its [`orientation`](@ref), such that the
sphere does not change when the shape is moved. See [`world_bounding_sphere`](@ref) for the sphere in
world coordinates.

The method is an optional part of the [`AbstractShape`](@ref) interface. With the default `nothing`
every ray is tested via [`intersect3d`](@ref). With a sphere, a ray that misses the sphere is not
tested against the shape at all, which pays off for shapes with a costly `intersect3d`.

The sphere must enclose every point of the shape, for all its parameters: a hit outside of the
sphere is lost without a warning. It should also be tight, since a ray that hits the sphere is
tested as before. A safety margin is not needed, the solver adds its own tolerance. Use
[`render_bounding_sphere!`](@ref) to look at the result.

```julia
# a cylinder of the radius `r` and the height `h`, with the origin at the center of its base
bounding_sphere(c::MyCylinder) = (Point3(0, c.h / 2, 0), sqrt(c.r^2 + (c.h / 2)^2))
```
"""
bounding_sphere(::AbstractShape) = nothing

"""
    world_bounding_sphere(shape::AbstractShape)

Returns `nothing`, or the [`bounding_sphere`](@ref) of the `shape` as a tuple `(center, radius)`
with the `center` in world coordinates, for the current position and orientation of the `shape`.
"""
world_bounding_sphere(shape::AbstractShape) = world_bounding_sphere(bounding_sphere(shape), shape)
world_bounding_sphere(::Nothing, ::AbstractShape) = nothing
function world_bounding_sphere(sphere::Tuple, shape::AbstractShape)
    center, radius = sphere
    return (position(shape) + orientation(shape) * center, radius)
end

"""
    translate3d!(::Movable, shape::AbstractShape, offset)

Translates the `pos`ition of `shape` by the `offset`-vector.
"""
function translate3d!(::Movable, shape::AbstractShape, offset)
    position!(shape, position(shape) + offset)
    return nothing
end

"""
    rotate3d!(::Movable, shape::AbstractShape, R::AbstractMatrix)

Rotates the `dir`-matrix of `shape` by the rotation matrix `R`.
"""
function rotate3d!(::Movable, shape::AbstractShape, R::AbstractMatrix)
    orientation!(shape, R * orientation(shape))
    return nothing
end
