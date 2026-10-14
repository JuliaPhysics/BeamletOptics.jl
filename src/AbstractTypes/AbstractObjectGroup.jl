"""
    AbstractObjectGroup

Container type for groups of optical elements, based on a tree-like data structure. Intended for easier kinematic handling of connected elements.
See also [`ObjectGroup`](@ref) for a concrete implementation.

A ray or beam interacts with the object of the group that it hits. A component whose parts interact
as one, e.g. a cemented doublet, is an [`AbstractObject`](@ref) with the [`MultiShape`](@ref) trait.

# Implementation reqs.

Subtypes of `AbstractObjectGroup` must implement the following:

## Fields:

- `objects`: a tuple of [`AbstractObject`](@ref)s or additional subgroups, allows for hierarchical structures
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
all static or all movable. The primitives `translate3d!(::Movable, group, offset)` and
`rotate3d!(::Movable, group, R::AbstractMatrix)` move all `objects` and the pose of the group.
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
    translate3d!(::Movable, group::AbstractObjectGroup, offset)

Moves all objects of the `group` and its `center` by `offset`.
"""
function translate3d!(::Movable, group::AbstractObjectGroup, offset)
    position!(group, position(group) .+ offset)
    foreach(member -> translate3d!(member, offset), objects(group))
    return nothing
end

"""
    rotate3d!(::Movable, group::AbstractObjectGroup, R::AbstractMatrix)

Rotates all objects of the `group` by `R` about the group `center` and updates the group
[`orientation`](@ref) to `R * orientation(group)`.
"""
function rotate3d!(::Movable, group::AbstractObjectGroup, R::AbstractMatrix)
    orientation!(group, R * orientation(group))
    pivot = position(group)
    foreach(member -> rotate3d!(member, R, pivot), objects(group))
    return nothing
end
