"""
    AbstractObjectGroup

Container type for groups of optical elements, based on a tree-like data structure.
Intended purely as a **kinematic** bundling convenience (e.g. moving a set of otherwise
unrelated components together as one rigid assembly) — it is deliberately *not* an
[`AbstractObject`](@ref): unlike a [`MultiShape`](@ref) `AbstractObject` (e.g. a
[`CubeBeamsplitter`](@ref) or [`DoubletLens`](@ref)), a group has no optical identity of
its own and is never a target of `intersect3d`/`interact3d` in a [`System`](@ref) —
`System`/`StaticSystem` transparently expand any `AbstractObjectGroup` (and any
`MultiShape` composite) down to its constituent `AbstractObject`s for tracing, see
[`objects(::System)`](@ref).
See also [`ObjectGroup`](@ref) for a concrete implementation.

# Implementation reqs.

Subtypes of `AbstractObjectGroup` must implement the following:

## Fields:

- `objects`: stores objects or additional subgroups of objects, allows for hierarchical structures

## Functions:

- the kinematic API (`translate3d!`, `translate_to3d!`, `rotate3d!`, `xrotate3d!`/`yrotate3d!`/`zrotate3d!`,
  `align3d!`, `reset_translation3d!`, `reset_rotation3d!`, `position`/`position!`, `orientation`/`orientation!`)
  is provided below by reusing [`MultiShape`](@ref)'s implementation directly (`shape(group) = objects(group)`
  makes an `AbstractObjectGroup` a valid "owner" for it) — subtypes only need to implement `objects(group)`,
  `position`/`position!` and `orientation`/`orientation!` themselves.
"""
abstract type AbstractObjectGroup{T} end

AbstractTrees.children(group::AbstractObjectGroup) = objects(group)

shape(group::AbstractObjectGroup) = objects(group)

translate3d!(group::AbstractObjectGroup, offset) = translate3d!(MultiShape(), group, offset)
translate_to3d!(group::AbstractObjectGroup, target) = translate_to3d!(MultiShape(), group, target)
rotate3d!(group::AbstractObjectGroup, axis, θ) = rotate3d!(MultiShape(), group, axis, θ)

xrotate3d!(group::AbstractObjectGroup, θ) = rotate3d!(group, Point3(1, 0, 0), θ)
yrotate3d!(group::AbstractObjectGroup, θ) = rotate3d!(group, Point3(0, 1, 0), θ)
zrotate3d!(group::AbstractObjectGroup, θ) = rotate3d!(group, Point3(0, 0, 1), θ)

align3d!(group::AbstractObjectGroup, target_vec) = align3d!(MultiShape(), group, target_vec)

reset_translation3d!(group::AbstractObjectGroup) = reset_translation3d!(MultiShape(), group)
reset_rotation3d!(group::AbstractObjectGroup) = reset_rotation3d!(MultiShape(), group)
