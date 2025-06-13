"""
    ObjectGroup <: AbstractObjectGroup

A tree-like storage container for groups of objects. Can store individual objects and subgroups.
Main purpose is handling of, i.e., groups of lenses.

## Fields

- `center`: a point in 3D space which is regarded as the reference origin of the group
- `dir`: a 3x3 matrix that describes the common `orientation` of the group
- `objects`: stores [`AbstractObject`](@ref), can also store subgroups of type [`AbstractObjectGroup`](@ref)

## Kinematic

A `ObjectGroup` implements the kinematic functions of [`AbstractObject`](@ref). The following logic is applied to

- [`translate3d!`](@ref): all objects in the group are translated by the offset vector
- [`translate_to3d!`](@ref): all objects are moved in parallel such that the group `center` is equal to the target position
- [`rotate3d!`](@ref): all objects are rotated around the `center` point with respect to their relative position
"""
mutable struct ObjectGroup{T, O <: Tuple{Vararg{AbstractObject}}} <: AbstractObjectGroup{T}
    dir::SMatrix{3, 3, T, 9}
    center::Point3{T}
    const objects::O
end

shape_trait_of(::ObjectGroup) = MultiShape()

shape(o::ObjectGroup) = o.objects

Base.position(group::ObjectGroup) = group.center
position!(group::ObjectGroup, pos) = (group.center = pos)

orientation(group::ObjectGroup) = group.dir
orientation!(group::ObjectGroup, dir) = (group.dir = dir)

ObjectGroup(v::AbstractArray, T = Float64) = ObjectGroup(tuple(v...), T)
function ObjectGroup(v::V, T = Float64) where {V <: Tuple}
    ObjectGroup{T, V}(SMatrix{3,3}(one(T)*I), Point3{T}(0), v)
end

"""
    objects(group::ObjectGroup)

Exposes all objects/subgroups stored within the group.
"""
objects(group::ObjectGroup) = group.objects

Base.show(::IO, ::MIME"text/plain", group::ObjectGroup) = print_tree(group)