"""
    AbstractShape{T<:Real}

A generic type for a shape that exists in 3D-space. Must have a `pos`ition and orientation.\\
Types used to describe the geometry of a shape should be subtypes of `Real`.\\

# Implementation reqs.

Subtypes of `AbstractShape` should implement the following:

## Fields:

- `pos`: a 3D-vector that stores the current `position` of the object-specific coordinate system
- `dir`: a 3x3-matrix that represents the orthonormal basis of the object and therefore, the `orientation`

## Getters/setters

- [`position`](@ref) / `position!`: gets or sets the `pos`ition vector of the `AbstractShape`
- [`orientation`](@ref) / `orientation!`: gets or sets the orientation matrix of the `AbstractShape`

## Kinematic:

- [`translate3d!`](@ref): the object is moved by a translation vector relative to its current position
- [`translate_to3d!`](@ref): the object is moved towards the target position
- [`rotate3d!`](@ref): the object is rotated by an angle around a reference vector
- [`xrotate3d!`](@ref): rotation around the x-axis
- [`yrotate3d!`](@ref): rotation around the y-axis
- [`zrotate3d!`](@ref): rotation around the z-axis
- [`align3d!`](@ref): align local shape y-axis with target vector
- [`reset_translation3d!`](@ref): return the `object` to the global origin
- [`reset_rotation3d!`](@ref): rotate the `object` back into its original state

## Ray Tracing:

- [`intersect3d`](@ref): returns the intersection between an `AbstractShape` and `AbstractRay`, or lack thereof. See also [`Intersection`](@ref)

## Rendering (with Makie):

Refer to the [`render!`](@ref) documentation.
"""
abstract type AbstractShape{T <: Real} end

"Enforces that `shape` has to have the field `pos` or implement `position()`."
Base.position(shape::AbstractShape) = shape.pos
position!(shape::AbstractShape, pos) = (shape.pos = pos)

"Enforces that `shape` has to have the field `dir` or implement `orientation()`."
orientation(shape::AbstractShape) = shape.dir
orientation!(shape::AbstractShape, dir) = (shape.dir = dir)

"""
    translate3d!(shape::AbstractShape, offset)

Translates the `pos`ition of `shape` by the `offset`-vector.
"""
function translate3d!(shape::AbstractShape, offset)
    position!(shape, position(shape) + offset)
    return nothing
end

"""
    translate_to3d!(shape::AbstractShape, target)

Translates the `shape` to the `target` position. 
"""
function translate_to3d!(shape::AbstractShape, target)
    current = position(shape)
    translate3d!(shape, target - current)
    return nothing
end

"""
    rotate3d!(shape::AbstractShape, axis, θ)

Rotates the `dir`-matrix of `shape` around the reference-`axis` by an angle of `θ`.
"""
function rotate3d!(shape::AbstractShape, axis, θ)
    R = rotate3d(axis, θ)
    orientation!(shape, R * orientation(shape))
    return nothing
end

"""Rotates the `dir`-matrix of `shape` around the global x-axis by an angle of `θ`."""
function xrotate3d!(shape::AbstractShape{T}, θ) where {T}
    rotate3d!(shape, Point3(one(T), zero(T), zero(T)), θ)
end
"""Rotates the `dir`-matrix of `shape` around the global y-axis by an angle of `θ`."""
function yrotate3d!(shape::AbstractShape{T}, θ) where {T}
    rotate3d!(shape, Point3(zero(T), one(T), zero(T)), θ)
end
"""Rotates the `dir`-matrix of `shape` around the global z-axis by an angle of `θ`."""
function zrotate3d!(shape::AbstractShape{T}, θ) where {T}
    rotate3d!(shape, Point3(zero(T), zero(T), one(T)), θ)
end

"""Returns the `shape` to the global origin."""
function reset_translation3d!(shape::AbstractShape{T}) where {T}
    position!(shape, Point3(zero(T)))
    return nothing
end

"""
    align3d!(shape, target_axis)

Rotates the `shape` such that its local y-axis aligns with the `target_axis`.
"""
function align3d!(shape::AbstractShape, target_axis)
    R = align3d(orientation(shape)[:,2], target_axis)
    orientation!(shape, R * orientation(shape))
    return nothing
end

"""Resets the `shape` rotation angles to zero."""
reset_rotation3d!(shape::AbstractShape{T}) where {T} = orientation!(shape, Matrix{T}(I, 3, 3))