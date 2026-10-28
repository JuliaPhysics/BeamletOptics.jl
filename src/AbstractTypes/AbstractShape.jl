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
- [`bounding_sphere_of`](@ref) (optional): computes a sphere that encloses the shape. A solve then skips
  the shape for every ray that misses the sphere, without a call of `intersect3d`

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

"""
    translate3d!(x, offset, carried::IdSet)
    rotate3d!(x, R::AbstractMatrix, pivot::AbstractVector, carried::IdSet)

Carried motion: moves `x` because its parent container moves, which is not a command to `x`
itself. Containers ([`ObjectGroup`](@ref)s, [`MultiShape`](@ref) objects, composite SDFs and
beam groups) move their members with these methods.

`carried` holds every part that the current rigid motion has already moved, compared by identity
(`===`). A part that can be reached along several paths, e.g. a shape listed twice by a
`MultiShape` object or an object in two subgroups, therefore moves exactly once. A container
starts a motion in its own primitive with an empty `IdSet{Any}()` and passes it on to its parts.

For a shape or a beam, the member types without parts, the carried method calls the public
`translate3d!(x, offset)` or `rotate3d!(x, R, pivot)` if `x` is not in `carried` yet. A container
implements its own carried methods: it adds itself to `carried`, rotates its own frame about
`pivot` (reading its `position` before its parts move) and carries its parts with the same
`carried`. An object type that overrides the primitives, e.g. an [`AbstractDetector`](@ref) that
empties itself, also overrides its carried methods, see [`AbstractObject`](@ref).
"""
function translate3d!(x::AbstractShape, offset, carried::IdSet)
    if _mark_carried!(carried, x)
        translate3d!(x, offset)
    end
    return nothing
end

function rotate3d!(x::AbstractShape, R::AbstractMatrix, pivot::AbstractVector, carried::IdSet)
    if _mark_carried!(carried, x)
        rotate3d!(x, R, pivot)
    end
    return nothing
end
