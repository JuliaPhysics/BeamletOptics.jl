"""
    AbstractEntity

A generic type for something that exists independently.

# Implementation reqs.
Subtypes of `AbstractObject` should implement the following:

# Fields
- `id`: a unique identifier (UUID 4)
"""
abstract type AbstractEntity end

id(entity::AbstractEntity) = entity.id

"Defines the intersection between an `entity` and something, defaults to no intersection."
function intersect3d(entity::AbstractEntity, ::Any)
    @warn lazy"No intersect3d method defined for:" typeof(entity)
    return nothing
end

"Defines the interaction between an `entity` and something, defaults to nothing which ends the ray tracer."
function interact3d(entity::AbstractEntity, ::Any)
    @warn lazy"No interact3d method defined for:" typeof(entity)
    return nothing
end

"""
    AbstractSystem <: AbstractEntity

A generic type for a container type which holds objects, rays, etc. 
"""
abstract type AbstractSystem <: AbstractEntity end

"""
    AbstractRay{T<:Real} <: AbstractEntity

A generic type for a ray/line in 3D-space. Must have a `pos`ition, `dir`ection and `len`gth, defined as 3D-vectors and a scalar.
"""
abstract type AbstractRay{T<:Real} <: AbstractEntity end

position(ray::AbstractRay) = ray.pos
position!(ray::AbstractRay, pos) = (ray.pos .= pos)

direction(ray::AbstractRay) = ray.dir
direction!(ray::AbstractRay, dir) = (ray.dir .= dir)

Base.length(ray::AbstractRay) = ray.len
length!(ray::AbstractRay, len) = (ray.len = len)

"""
    AbstractObject <: AbstractEntity

A generic type for 2D/3D objects. Optical elements are supposed to fall under this type. Must have a `shape` field which defines the
volume or surface with which a ray can interact.
"""
abstract type AbstractObject <: AbstractEntity end

shape(o::AbstractObject) = o.shape

"""
    intersect3d(object::AbstractObject, ray::AbstractRay)

In general, the intersection between an `object` and a `ray` is defined as the intersection with `shape(object)`.
"""
function intersect3d(object::AbstractObject, ray::AbstractRay)
    intersection = intersect3d(shape(object), ray)
    # Ensure that the intersection knows about the object id
    if !isnothing(intersection)
        intersection.id = id(object)
    end
    return intersection
end

"""
    AbstractShape{T<:Real} <: AbstractEntity

A generic type for a shape that exists in 3D-space. Must have a `pos`ition and `dir`ection.\\
Types used to describe an `object` should be subtypes of `Real`.\\

# Implementation reqs.
Subtypes of `AbstractShape` should implement the following:

## Fields:
- `pos`: a 3D-vector that stores the current `position` of the object-specific coordinate system
- `dir`: a 3x3-matrix that represents the orthonormal basis of the object and therefore, the `orientation`

## Kinematic:
- `translate3d!`: the object is moved by a translation vector relative to its current position
- `rotate3d!`: the object is rotated by an angle around a reference vector
- `xrotate3d!`: rotation around the x-axis
- `yrotate3d!`: rotation around the y-axis
- `zrotate3d!`: rotation around the z-axis
- `reset_translation3d!`: return the `object` to the global origin
- `reset_rotation3d!`: rotate the `object` back into its original state

## Ray Tracing:
- `intersect3d`: returns the intersection between an `AbstractShape` and `AbstractRay`, or lack thereof. See also `Intersection{T}`

## Rendering (with GLMakie):
- `render_shape!`: plot the `shape` into an `Axis3` environment
- `render_shape_normals!`: plot the `shape` surface normals into an `Axis3` environment (optional)
"""
abstract type AbstractShape{T<:Real} <: AbstractEntity end

"Enforces that `object` has to have the field `pos` or implement `position()`."
position(shape::AbstractShape) = shape.pos
position!(shape::AbstractShape, pos) = (shape.pos .= pos)

"Enforces that `object` has to have the field `dir` or implement `orientation()`."
orientation(shape::AbstractShape) = shape.dir
orientation!(shape::AbstractShape, dir) = (shape.dir .= dir)

"""
    translate3d!(shape::AbstractShape, offset)

Translates the `pos`ition of `shape` by the `offset`-vector.
"""
function translate3d!(shape::AbstractShape, offset)
    shape.pos .+= offset
    return nothing
end

"""
    rotate3d!(shape::AbstractShape, axis, θ)

Rotates the `dir`-matrix of `shape` around the reference-`axis` by an angle of `θ`.
"""
function rotate3d!(shape::AbstractShape, axis, θ)
    R = rotate3d(axis, θ)
    orientation!(shape, orientation(shape) * R)
    return nothing
end

xrotate3d!(shape::AbstractShape{T}, θ) where T = rotate3d!(shape, @SVector(T[one(T), zero(T), zero(T)]), θ)
yrotate3d!(shape::AbstractShape{T}, θ) where T = rotate3d!(shape, @SVector(T[zero(T), one(T), zero(T)]), θ)
zrotate3d!(shape::AbstractShape{T}, θ) where T = rotate3d!(shape, @SVector(T[zero(T), zero(T), one(T)]), θ)

function reset_translation3d!(shape::AbstractShape{T}) where T
    shape.pos .= zeros(T, 3)
    return nothing
end

reset_rotation3d!(shape::AbstractShape{T}) where T = (shape.dir .= Matrix{T}(I, 3, 3))

render_shape!(::Any, ::AbstractShape) = nothing
render_shape_normals!(::Any, ::AbstractShape) = nothing

"""
    isinfrontof(shape::AbstractShape, ray::AbstractRay)

A simple test to check if a `shape` lies "in front of" a `ray`.
The forward direction is here defined as the ray `orientation`.
Only works well if `ray` is **outside** of the volume of `shape`.
Can be dispatched to return more accurate results for subtypes of `AbstractShape`.
"""
isinfrontof(shape::AbstractShape, ray::AbstractRay) = isinfrontof(position(shape), position(ray), direction(ray))
