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
- [`set_pivot3d!`](@ref): moves the group `center` (the pivot used above) without moving any of its `objects`
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

"""
    set_pivot3d!(group::ObjectGroup, pivot)

Moves the kinematic pivot (`center`) of `group` to `pivot`, without moving any of its
`objects`. The pivot is the reference point used by [`rotate3d!`](@ref),
[`translate_to3d!`](@ref) and [`reset_translation3d!`](@ref); the group [`orientation`](@ref)
is unchanged.

# Example

Rotate a lens group about its first surface instead of the group's default center:

```julia
group = ObjectGroup([lens1, lens2])
set_pivot3d!(group, position(lens1))
rotate3d!(group, [0, 0, 1], deg2rad(5))   # rotates about lens1's position, not the old center
```
"""
function set_pivot3d!(group::ObjectGroup{T}, pivot) where {T}
    position!(group, Point3{T}(pivot))
    return nothing
end

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
