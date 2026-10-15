"""
    ObjectGroup <: AbstractObjectGroup

A tree-like storage container for groups of objects. Can store individual objects and subgroups.
Main purpose is handling of, i.e., groups of lenses.

## Fields

- `center`: a point in 3D space which is regarded as the reference origin of the group
- `dir`: a 3x3 matrix that describes the common `orientation` of the group
- `objects`: stores [`AbstractObject`](@ref), can also store subgroups of type [`AbstractObjectGroup`](@ref)

Each of its `objects` is traced and interacts on its own, read them via [`objects`](@ref)`(group)`.

## Kinematic

An `ObjectGroup` implements the kinematic API, see [`BeamletOptics.AbstractKinematicTrait`](@ref). The following logic is applied to

- [`translate3d!`](@ref): all objects in the group are translated by the offset vector
- [`translate_to3d!`](@ref): all objects are moved in parallel such that the group `center` is equal to the target position
- [`rotate3d!`](@ref): all objects are rotated around the `center` point with respect to their relative position
- [`set_pivot3d!`](@ref): moves the group `center` (the pivot used above) without moving any of its `objects`

The `objects` must be either all static or all movable, otherwise the constructor throws an `ArgumentError`.
The group takes that kinematic class, see [`BeamletOptics.AbstractKinematicTrait`](@ref).
"""
mutable struct ObjectGroup{T, O <: Tuple{Vararg{ObjectOrGroup}}} <: AbstractObjectGroup{T}
    dir::SMatrix{3, 3, T, 9}
    center::Point3{T}
    const objects::O
end

ObjectGroup(v::AbstractArray, T = Float64) = ObjectGroup(tuple(v...), T)
function ObjectGroup(v::V, T = Float64) where {V <: Tuple}
    _check_kinematic_members(v)
    ObjectGroup{T, V}(SMatrix{3,3}(one(T)*I), Point3{T}(0), v)
end

Base.show(::IO, ::MIME"text/plain", group::ObjectGroup) = print_tree(group)
