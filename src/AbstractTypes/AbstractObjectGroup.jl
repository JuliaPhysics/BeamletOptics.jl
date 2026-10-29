"""
    AbstractObjectGroup

Container type for groups of optical elements, based on a tree-like data structure. Intended for easier kinematic handling of connected elements.
See also [`ObjectGroup`](@ref) for a concrete implementation.

A ray or beam interacts with the object of the group that it hits. A component whose parts interact
as one, e.g. a cemented doublet, is an [`AbstractObject`](@ref) with the [`MultiShape`](@ref) trait.

!!! tip "Own assembly types"
    An assembly of separate parts, e.g. a telescope from its mirrors, can be a subtype of
    `AbstractObjectGroup` with named parts and its own constructor: it needs the fields below and
    no `interact3d` method, and the kinematic API, [`set_pivot3d!`](@ref), systems and rendering
    work as for an [`ObjectGroup`](@ref). See `templates/08_custom_group.jl` of the agent skill.

# Implementation reqs.

Subtypes of `AbstractObjectGroup` must implement the following:

## Fields:

- `objects`: a tuple or a vector of [`AbstractObject`](@ref)s or additional subgroups, allows for hierarchical structures.
  Only [`ObjectGroup`](@ref) can be changed after its construction, via `push!` and `delete!`.
- `center`: a `Point3{T}` which is regarded as the reference origin (pivot) of the group
- `dir`: a `SMatrix{3,3,T,9}` that describes the local coordinate system of the group

Since the kinematic API modifies `center` and `dir`, the subtype must be a `mutable struct`.

## Functions:

If the fields above do not exist, the following getters/setters must be dispatched:

- [`objects`](@ref): getter for the `objects` field or equivalent return type
- `position` / `position!`: gets or sets the reference origin (pivot)
- [`orientation`](@ref) / `orientation!`: gets or sets the orientation matrix

## Kinematic

An `AbstractObjectGroup` is a container: it takes the kinematic class of its `objects`, i.e.
[`BeamletOptics.Movable`](@ref) with an [`BeamletOptics.Oriented`](@ref) frame for movable objects, see
[`BeamletOptics.AbstractKinematicTrait`](@ref). The constructors must check that the `objects` are either
all static or all movable, e.g. with `allequal(BeamletOptics.is_static, objects)`. The primitives `translate3d!(::Movable, group, offset)` and
`rotate3d!(::Movable, group, R::AbstractMatrix)` move all `objects` and the pose of the group. A group that is
nested in another group moves as a frame (`position!`/`orientation!`), its own primitives are not called then.
"""
abstract type AbstractObjectGroup{T} end

"""
    ObjectOrGroup

Union of [`AbstractObject`](@ref) and [`AbstractObjectGroup`](@ref): what a [`System`](@ref) or a
group stores and what [`render!`](@ref) draws.
"""
const ObjectOrGroup = Union{AbstractObject, AbstractObjectGroup}

"""
    objects(group::AbstractObjectGroup)

Exposes all objects/subgroups stored within the group.
"""
objects(group::AbstractObjectGroup) = group.objects

AbstractTrees.children(group::AbstractObjectGroup) = objects(group)

"""
    _leaf_objects(tree)

Iterates the objects at all depths of `tree`, an object, an object group or a collection of them.
`Leaves(tree)` without the nodes that have no children but are no objects, i.e. empty groups and an
empty collection.
"""
_leaf_objects(tree) = Iterators.filter(_is_object, Leaves(tree))

_is_object(::AbstractObject) = true
_is_object(::Any) = false

# The type of `_leaf_objects(tree)`, e.g. of the objects of a `System`, see `intersect3d(parts, ray)`
const _LeafObjects = Iterators.Filter{typeof(_is_object), <:Leaves}

kinematic_trait_of(group::AbstractObjectGroup) = _container_trait(objects(group))

Base.position(group::AbstractObjectGroup) = group.center
position!(group::AbstractObjectGroup{T}, pos) where {T} = (group.center = Point3{T}(pos))

orientation(group::AbstractObjectGroup) = group.dir
orientation!(group::AbstractObjectGroup{T}, dir) where {T} = (group.dir = SMatrix{3, 3, T, 9}(dir))

function initialize!(group::AbstractObjectGroup)
    foreach(initialize!, objects(group))
    return nothing
end

"""
    set_pivot3d!(group::AbstractObjectGroup, pivot)

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
function set_pivot3d!(group::AbstractObjectGroup, pivot)
    position!(group, pivot)
    return nothing
end

"""
    _unwind(group) -> (subgroups, leaves)

Unwinds the tree of `group` for one rigid motion, each node once, compared by identity (`===`):
the `subgroups`, including `group`, move their `center` and `dir`; the `leaves` are its objects,
which are moved with their public verbs, see `_movable_parts`.
"""
function _unwind(group::AbstractObjectGroup)
    subgroups = unique(objectid, Iterators.filter(!_is_object, PreOrderDFS(group)))
    leaves = _movable_parts(collect(_leaf_objects(group)))
    return subgroups, leaves
end

"""
    translate3d!(::Movable, group::AbstractObjectGroup, offset)

Moves all objects of the `group` and its `center` by `offset`. Every object and subgroup moves
once, also if the `group` holds it along several paths, see the containers in
[`BeamletOptics.AbstractKinematicTrait`](@ref).
"""
function translate3d!(::Movable, group::AbstractObjectGroup, offset)
    subgroups, leaves = _unwind(group)
    foreach(subgroup -> position!(subgroup, position(subgroup) .+ offset), subgroups)
    foreach(object -> translate3d!(object, offset), leaves)
    return nothing
end

"""
    rotate3d!(::Movable, group::AbstractObjectGroup, R::AbstractMatrix)

Rotates all objects of the `group` by `R` about the group `center` and updates the group
[`orientation`](@ref) to `R * orientation(group)`. Every object and subgroup moves once, see
`translate3d!(::Movable, group::AbstractObjectGroup, offset)`.
"""
function rotate3d!(::Movable, group::AbstractObjectGroup, R::AbstractMatrix)
    pivot = position(group)
    subgroups, leaves = _unwind(group)
    for subgroup in subgroups
        # Rotate center about pivot
        p = position(subgroup)
        orientation!(subgroup, R * orientation(subgroup))
        position!(subgroup, pivot + R * (p - pivot))
    end
    foreach(object -> rotate3d!(object, R, pivot), leaves)
    return nothing
end
