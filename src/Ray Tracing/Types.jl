"""
    AbstractEntity

A generic type for something that exists independently.
"""
abstract type AbstractEntity end

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

"Defines the interaction between a `ray` and an `entity`."
interact3d(::AbstractEntity, ::AbstractRay{T}) where T = nothing

"""
    AbstractObject{T<:Real} <: AbstractEntity

A generic type for something that exists in 3D-space. Must have a `pos`ition and `dir`ection.\\
Types used to describe an `object` should be subtypes of `Real`.\\

# Implementation reqs.
Subtypes of `AbstractObject` should implement the following:

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
- `intersect3d`: returns the intersection between an `AbstractObject` and `AbstractRay`, or lack thereof. See also `Intersection{T}`
- `interact3d`: returns an `Interaction{T}`. See also `Interaction{T}`.

## Rendering (with GLMakie):
- `render_object!`: plot the `object` into an `Axis3` environment
- `render_object_normals!`: plot the `object` surface normals into an `Axis3` environment (optional)
"""
abstract type AbstractObject{T<:Real} <: AbstractEntity end

"Enforces that `object` has to have the field `pos` or implement `position()`."
position(object::AbstractObject) = object.pos
position!(object::AbstractObject, pos) = (object.pos .= pos)

"Enforces that `object` has to have the field `dir` or implement `orientation()`."
orientation(object::AbstractObject) = object.dir
orientation!(object::AbstractObject, dir) = (object.dir .= dir)

"""
    translate3d!(object::AbstractObject, offset)

Translates the `pos`ition of `object` by the `offset`-vector.
"""
function translate3d!(object::AbstractObject, offset)
    object.pos .+= offset
    return nothing
end

"""
    rotate3d!(object::AbstractObject, axis, θ)

Rotates the `dir`-matrix of `object` around the reference-`axis` by an angle of `θ`.
"""
function rotate3d!(object::AbstractObject, axis, θ)
    R = rotate3d(axis, θ)
    orientation!(object, orientation(object) * R)
    return nothing
end

xrotate3d!(object::AbstractObject{T}, θ) where T = rotate3d!(object, T[one(T), zero(T), zero(T)], θ)
yrotate3d!(object::AbstractObject{T}, θ) where T = rotate3d!(object, T[zero(T), one(T), zero(T)], θ)
zrotate3d!(object::AbstractObject{T}, θ) where T = rotate3d!(object, T[zero(T), zero(T), one(T)], θ)

function reset_translation3d!(object::AbstractObject{T}) where T
    object.pos .= zeros(T, 3)
    return nothing
end

reset_rotation3d!(object::AbstractObject{T}) where T = (object.dir .= Matrix{T}(I, 3, 3))

function intersect3d(object::AbstractObject{O}, ::AbstractRay{R}) where {O, R}
    @warn lazy"No intersect3d method defined for:" typeof(object)
    return nothing
end

function interact3d(object::AbstractObject{O}, ::AbstractRay{R}) where {O, R}
    @warn lazy"No interact3d method defined for:" typeof(object)
    return nothing
end

render_object!(::Any, ::AbstractObject) = nothing
render_object_normals!(::Any, ::AbstractObject) = nothing

"""
    isinfrontof(object::AbstractObject, ray::AbstractRay)

A simple test to check if an `object` lies "in front of" a `ray`.
The forward direction is here defined as the ray `orientation`.
Only works well if `ray` is **outside** of the volume of `object`.
Can be dispatched to return more accurate results for subtypes of `AbstractObject`
"""
function isinfrontof(object::AbstractObject, ray::AbstractRay)
    los = position(object) - position(ray)
    if fast_dot3d(direction(ray), los) <= 0
        return false
    else
        return true
    end
end
