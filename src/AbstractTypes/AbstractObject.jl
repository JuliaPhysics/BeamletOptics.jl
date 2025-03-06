"""
    AbstractObject

A generic type for 2D/3D objects that can be used to model optical elements. The geometry of the object is represented via a [`AbstractShape`](@ref).
The optical effect that occurs between the object and an incoming ray/beam of light is modeled via its [`interact3d`](@ref) method.

# Implementation reqs.

Subtypes of `AbstractObject` must implement the following:

## Shape trait

An `AbstractObject` can consist of a single [`AbstractShape`](@ref), e.g. a lens element, or a collection of functionally dependant shapes, e.g. a cube beamsplitter.
In order to model this, the API implementation of an `AbstractObject` requires the definition of the **shape trait**. This trait allows the dispatch onto specialized methods to handle
the kinematic interface and tracing methods for objects consisting of one or more shapes.

- `shape_trait_of`: defines the shape type of the `AbstractObject`, refer to [`AbstractShapeTrait`](@ref) for more information

!!! info "Default shape trait"
    Unless specified otherwise, the `shape_trait_of` an `AbstractObject` is defined as [`SingleShape`](@ref). This requires `object.shape` as a dedicated field.
    For [`MultiShape`](@ref)s the getter function `shape(object)` must return a tuple of all shapes that make up the object.

## Getters/setters

All kinematic functions defined for the [`AbstractShape`](@ref) can also be called for a `AbstractObject`. In this case, the shape trait will define how the specific movement function
is dispatched.

## Functions:

- [`interact3d`](@ref): defines the optical interaction, the return type must be `Nothing` or an [`AbstractInteraction`](@ref)
"""
abstract type AbstractObject{T <: Real, S <: AbstractShape{T}} end

"Default trait"
shape_trait_of(::AbstractObject) = SingleShape()

"Dispatch the shape function based on the [`AbstractShapeTrait`](@ref) of the [`AbstractObject`](@ref)"
shape(object::AbstractObject) = shape(shape_trait_of(object), object)

"Enforces that `object` has to have the field `pos` or implement `position()`."
position(object::AbstractObject) = position(shape_trait_of(object), object)
position!(object::AbstractObject, pos) = position!(shape_trait_of(object), object, pos)

"Enforces that `object` has to have the field `dir` or implement `orientation()`."
orientation(object::AbstractObject) = orientation(shape_trait_of(object), object)
orientation!(object::AbstractObject, dir) = orientation!(shape_trait_of(object), object, dir)

translate3d!(object::AbstractObject, offset) = translate3d!(shape_trait_of(object), object, offset)

translate_to3d!(object::AbstractObject, target) = translate_to3d!(shape_trait_of(object), object, target)

rotate3d!(object::AbstractObject, axis, θ) = rotate3d!(shape_trait_of(object), object, axis, θ)

xrotate3d!(object::AbstractObject, θ) = rotate3d!(object, Point3(1, 0, 0), θ)
yrotate3d!(object::AbstractObject, θ) = rotate3d!(object, Point3(0, 1, 0), θ)
zrotate3d!(object::AbstractObject, θ) = rotate3d!(object, Point3(0, 0, 1), θ)

align3d!(object::AbstractObject, axis) = align3d!(shape_trait_of(object), object, axis)

reset_translation3d!(object::AbstractObject) = reset_translation3d!(shape_trait_of(object), object)

reset_rotation3d!(object::AbstractObject) = reset_rotation3d!(shape_trait_of(object), object)

"""
    AbstractObjectGroup <: AbstractObject

Container type for groups of optical elements, based on a tree-like data structure. Intended for easier kinematic handling of connected elements.
See also [`ObjectGroup`](@ref) for a concrete implementation.

# Implementation reqs.

Subtypes of `AbstractObjectGroup` must implement the following:

## Fields:

- `objects`: stores objects or additional subgroups of objects, allows for hierarchical structures

## Functions:

- for the kinematic API, all corresponding functions of [`AbstractObject`](@ref) must be implemented
"""
abstract type AbstractObjectGroup{T} <: AbstractObject{T, AbstractShape{T}} end

AbstractTrees.children(group::AbstractObjectGroup) = group.objects
