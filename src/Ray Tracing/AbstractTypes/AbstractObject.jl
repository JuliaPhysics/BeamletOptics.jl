"""
    AbstractObject

A generic type for 2D/3D objects that can be used to model optical elements. The geometry of the object is represented via a [`AbstractShape`](@ref).
The optical effect that occurs between the object and an incoming ray/beam of light is modeled via its [`interact3d`](@ref) method.

# Implementation reqs.

Subtypes of `AbstractObject` must implement the following:

## Fields:

- `shape`: stores an [`AbstractShape`](@ref) that represents the object geometry, called via `shape(object)`

## Getters/setters

- [`position`](@ref) / `position!`: gets or sets the `pos`ition vector of the `AbstractObject`
- [`orientation`](@ref) / `orientation!`: gets or sets the orientation matrix of the `AbstractObject`

## Functions:

- [`intersect3d`](@ref): returns the intersection between the object and an incoming ray, defaults to the intersection with `shape(object)`
- [`interact3d`](@ref): defines the optical interaction, should return `nothing` or an [`AbstractInteraction`](@ref)
- for the kinematic API, all corresponding functions should be forwarded to the underlying [`AbstractShape`](@ref) i.e. `rotate3d!(shape(object))`
"""
abstract type AbstractObject{T, S <: AbstractShape{T}} end

shape(object::AbstractObject) = object.shape

"Enforces that `object` has to have the field `pos` or implement `position()`."
position(object::AbstractObject) = position(shape(object))
position!(object::AbstractObject) = position!(shape(object))

"Enforces that `object` has to have the field `dir` or implement `orientation()`."
orientation(object::AbstractObject) = orientation(shape(object))
orientation!(object::AbstractObject) = orientation!(shape(object))

translate3d!(object::AbstractObject, offset) = translate3d!(shape(object), offset)

translate_to3d!(object::AbstractObject, target) = translate_to3d!(shape(object), target)

rotate3d!(object::AbstractObject, axis, θ) = rotate3d!(shape(object), axis, θ)

xrotate3d!(object::AbstractObject, θ) = rotate3d!(object, Point3(1, 0, 0), θ)
yrotate3d!(object::AbstractObject, θ) = rotate3d!(object, Point3(0, 1, 0), θ)
zrotate3d!(object::AbstractObject, θ) = rotate3d!(object, Point3(0, 0, 1), θ)

reset_translation3d!(object::AbstractObject) = reset_translation3d!(shape(object))

reset_rotation3d!(object::AbstractObject) = reset_rotation3d!(shape(object))

"""
    AbstractObjectGroup <: AbstractObject

Container type for groups of optical elements, based on a tree-like data structure. Intended for easier kinematic handling of connected elements.

# Implementation reqs.

Subtypes of `AbstractObjectGroup` must implement the following:

## Fields:

- `objects`: stores objects or additional subgroups of objects, allows for hierarchical structures

## Functions:

- for the kinematic API, all corresponding functions of [`AbstractObject`](@ref) must be implemented
"""
abstract type AbstractObjectGroup{T} <: AbstractObject{T, AbstractShape{T}} end

AbstractTrees.children(group::AbstractObjectGroup) = group.objects
