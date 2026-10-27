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

## Changing the objects

Objects and subgroups are added to and removed from the top level of the group via `push!`,
`delete!`, `pop!` and `popat!`, like those of a [`System`](@ref). The constructor copies the vector or tuple it is
given, such that these do not change it.
"""
mutable struct ObjectGroup{T} <: AbstractObjectGroup{T}
    dir::SMatrix{3, 3, T, 9}
    center::Point3{T}
    const objects::Vector{ObjectOrGroup}
end

function ObjectGroup(v::Union{AbstractArray, Tuple}, T = Float64)
    for x in v
        x isa ObjectOrGroup || throw(TypeError(:ObjectGroup, "objects", ObjectOrGroup, x))
    end
    _check_kinematic_members(v)
    ObjectGroup{T}(SMatrix{3,3}(one(T)*I), Point3{T}(0), collect(ObjectOrGroup, v))
end

"""
    push!(group::ObjectGroup, objects::ObjectOrGroup...) -> group

Adds the `objects` (or object groups) at the top level of the `group`. The pose of the `group`
(its pivot and [`orientation`](@ref)) and the poses of the `objects` stay as they are, i.e. the
`objects` are added where they are.

Throws an `ArgumentError` and adds nothing if an object is already part of the `group`, directly
or within a subgroup, or if the `objects` and those of the `group` would not be all static or all
movable, see [`BeamletOptics.AbstractKinematicTrait`](@ref). Objects are compared by identity.

The `group` does not know the systems that hold it: an object that is also part of such a system
is then traced twice, so remove it from the system first, e.g. via `delete!(system, object)`.
Beams that were solved before must be solved again, see `push!(system, objects...)`.
"""
function Base.push!(group::ObjectGroup, objs::ObjectOrGroup...)
    # objects and subgroups, also the `group` itself, such that a group can not hold itself
    nodes = ObjectOrGroup[o for o in PreOrderDFS(group)]
    for obj in objs, node in PreOrderDFS(obj)
        _is_leaf_of(node, nodes) && throw(ArgumentError(
            "the $(nameof(typeof(node))) is already part of the group"))
        push!(nodes, node)
    end
    _check_kinematic_members((group.objects..., objs...))
    append!(group.objects, objs)
    return group
end

"""
    delete!(group::ObjectGroup, object::ObjectOrGroup) -> group

Removes the top-level `object` (or subgroup) from the `group`, compared by identity. The pose of
the `group` and the pose of the `object` stay as they are. Nothing happens for an `object` that is
not part of the `group`. An object within a subgroup throws an `ArgumentError`: remove it from
its subgroup instead. A group without objects remains a valid, empty group.

Beams that were solved before must be solved again, see `delete!(system, object)`.
"""
function Base.delete!(group::ObjectGroup, obj::ObjectOrGroup)
    i = findfirst(o -> o === obj, group.objects)
    if isnothing(i)
        for sub in group.objects
            sub isa AbstractObjectGroup || continue
            any(o -> o === obj, PreOrderDFS(sub)) && throw(ArgumentError(
                "the $(nameof(typeof(obj))) is part of a subgroup of the group, remove it from its subgroup instead"))
        end
        return group
    end
    deleteat!(group.objects, i)
    return group
end

"""
    pop!(group::ObjectGroup) -> ObjectOrGroup

Removes the last top-level object (or subgroup) of the `group` and returns it. The pose of the
`group` and the pose of the object stay as they are, see `delete!(group, object)`. An empty
`group` throws an `ArgumentError`.
"""
Base.pop!(group::ObjectGroup) = pop!(group.objects)

"""
    popat!(group::ObjectGroup, i::Integer) -> ObjectOrGroup

Removes the `i`-th top-level object (or subgroup) of the `group` and returns it. The index counts
the objects and subgroups as they were added, a subgroup counts as one. An index out of bounds
throws a `BoundsError`. The pose of the `group` and the pose of the object stay as they are, see
`delete!(group, object)`.
"""
Base.popat!(group::ObjectGroup, i::Integer) = popat!(group.objects, i)

Base.show(::IO, ::MIME"text/plain", group::ObjectGroup) = print_tree(group)
