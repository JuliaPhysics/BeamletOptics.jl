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
    # Heterogeneous: objects and subgroups
    const objects::Vector{ObjectOrGroup}
end

function ObjectGroup(v::Union{AbstractArray, Tuple}, T = Float64)
    for x in v
        x isa ObjectOrGroup || throw(TypeError(:ObjectGroup, "objects", ObjectOrGroup, x))
    end
    check_kinematic_members(v)
    ObjectGroup{T}(SMatrix{3,3}(one(T)*I), Point3{T}(0), collect(ObjectOrGroup, v))
end

"""
    push!(group::ObjectGroup, objects::ObjectOrGroup...) -> group

Adds the `objects` (or object groups) at the top level of the `group`. The pose of the `group`
(its pivot and [`orientation`](@ref)) and the poses of the `objects` stay as they are, i.e. the
`objects` are added where they are.

The `objects` join the current frame of the `group`, like parts put onto a mount that was moved
before: [`reset_rotation3d!`](@ref) and [`reset_translation3d!`](@ref) of the `group` undo its
rotation and translation so far and move the added `objects` by the same amount, although these
did not take part in that motion.

Throws an `ArgumentError` and adds nothing if an object is already part of the `group`, directly
or within a subgroup, or if the `objects` and those of the `group` would not be all static or all
movable, see [`BeamletOptics.AbstractKinematicTrait`](@ref). Objects are compared by identity.

The `group` does not know the systems that hold it: an object that is also part of such a system
is then traced twice, so remove it from the system first, e.g. via `delete!(system, object)`.
Beams that were solved before must be solved again, see `push!(system, objects...)`. A
[`StaticSystem`](@ref) keeps the objects that its groups had when it was built and does not trace
the new `objects`: trace the `group` in a `System`, or build a new `StaticSystem`. The `group`
must not be changed while a system that holds it is being solved.
"""
function Base.push!(group::ObjectGroup, objects::ObjectOrGroup...)
    # Collect all nodes, including the group
    nodes = ObjectOrGroup[o for o in PreOrderDFS(group)]
    for object in objects, node in PreOrderDFS(object)
        _is_leaf_of(node, nodes) && throw(ArgumentError(
            "the $(nameof(typeof(node))) is already part of the group"))
        push!(nodes, node)
    end
    check_kinematic_members((group.objects..., objects...))
    append!(group.objects, objects)
    return group
end

"""
    delete!(group::ObjectGroup, object::ObjectOrGroup) -> group

Removes the top-level `object` (or subgroup) from the `group`, compared by identity. The pose of
the `group` and the pose of the `object` stay as they are. Nothing happens for an `object` that is
not part of the `group`. An object within a subgroup throws an `ArgumentError`: remove it from
its subgroup instead. A group without objects remains a valid, empty group.

Beams that were solved before must be solved again, see `delete!(system, object)`. A
[`StaticSystem`](@ref) keeps the objects that its groups had when it was built and still traces the
removed `object`: trace the `group` in a `System`, or build a new `StaticSystem`. The `group`
must not be changed while a system that holds it is being solved.
"""
function Base.delete!(group::ObjectGroup, object::ObjectOrGroup)
    i = findfirst(o -> o === object, group.objects)
    if isnothing(i)
        # Search the subgroups
        for sub in group.objects
            any(o -> o === object, PreOrderDFS(sub)) && throw(ArgumentError(
                "the $(nameof(typeof(object))) is part of a subgroup of the group, remove it from its subgroup instead"))
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
