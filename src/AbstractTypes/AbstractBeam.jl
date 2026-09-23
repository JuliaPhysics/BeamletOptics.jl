"""
    AbstractBeam{T <: Real, R <: AbstractRay{T}}

A generic type for a container type which holds rays, beams etc.

# Parametrization:

A subtype of `AbstractBeam` is parameterized by its main data type `T <: Real`, as well as the underlying ray representation `R <: AbstractRay{T}`.
If a beam is to be compatible with different [`AbstractRay`](@ref) implementations, it must be parameterized by `T` and `R`.
However, it can also be set to a fixed type for `T` and `R`, i.e. `MyBeam <: AbstractBeam{Float32, MyRay}`.

# Implementation reqs.

Subtypes of `AbstractBeam` must implement the following:

## Fields:

- `parent`: a [`Nullable`](@ref) field that holds the same type as the subtype, used for tree navigation
- `children`: a vector that holds the same type as the subtype, used for sub-beam tracking, i.e. beamsplitting

## Functions:

- `_modify_beam_head!`: modifies the beam path for retracing purposes
- `_last_beam_intersection`: returns the last `Beam` intersection
- `empty!`: resets the beam to its unsolved state
- `first_ray`: returns the start ray on the optical axis of the beam; for beamlets the first chief ray.
  Defines the generic `position`/`direction` of the beam (the pivot for rotations)

The tree functions `parent`, `children` and `isroot` (from `AbstractTrees`) work via the `parent`/`children` fields.

## Kinematic API:

`AbstractBeam`s are [`BeamletOptics.Movable`](@ref) with a [`BeamletOptics.Directed`](@ref) frame, see
[`BeamletOptics.AbstractKinematicTrait`](@ref). Only root beams can be moved.
Every move resets the beam to its untraced start state via `empty!`; moving a child beam throws an `ArgumentError`.
Rotations are applied about the `position` of the beam, i.e. the start of `first_ray`.
[`reset_rotation3d!`](@ref) throws an `ArgumentError`, since a beam has no orientation.

To support the kinematic API, a subtype additionally implements one of the following:

- `_component_beams`: for beams made up of component [`Beam`](@ref)s (e.g. [`GaussianBeamlet`](@ref)),
  returns a tuple of these beams. The generic `translate3d!`/`rotate3d!` then delegate to them.
- `translate3d!(::Movable, beam, offset)` and `rotate3d!(::Movable, beam, R::AbstractMatrix)`: for beams that
  store rays directly (e.g. [`Beam`](@ref)), move the start ray(s) via the ray verbs.
"""
abstract type AbstractBeam{T <: Real, R <: AbstractRay{T}} end

kinematic_trait_of(::AbstractBeam) = Movable(Directed())

AbstractTrees.NodeType(::Type{T}) where {T <: AbstractBeam} = HasNodeType()
AbstractTrees.nodetype(::Type{T}) where {T <: AbstractBeam} = T

AbstractTrees.ParentLinks(::Type{<:AbstractBeam}) = AbstractTrees.StoredParents()
AbstractTrees.parent(beam::AbstractBeam) = beam.parent
parent!(beam::B, parent::B) where {B <: AbstractBeam} = (beam.parent = parent)
# `isroot(beam)` is provided by AbstractTrees via `parent`, i.e. `isnothing(parent(beam))`

AbstractTrees.children(b::AbstractBeam) = b.children

AbstractTrees.printnode(io::IO, node::B; kw...) where {B <: AbstractBeam} = show(io, B)

"""
    children!(beam::B, child::B) where {B<:AbstractBeam}

Handles the inclusion of adding a single `child` to an existing `beam`. The function behaves as follows:

1. If no previous children exist, add child
2. If `beam` already has a single child, modify child beam starting ray (retracing)
3. Else throw error
"""
function children!(beam::B, child::B) where {B <: AbstractBeam}
    if isempty(children(beam))
        # Link parent and add child to tree
        parent!(child, beam)
        push!(children(beam), child)
        return nothing
    end
    if length(children(beam)) == 1
        _modify_beam_head!(first(children(beam)), child)
        return nothing
    end
    return error("Adding child to beam failed")
end

function children!(beam::B, _children::AbstractVector{B}) where {B <: AbstractBeam}
    if isempty(children(beam))
        # Link parent and add children to tree
        parent!.(_children, Ref(beam))
        append!(children(beam), _children)
        return nothing
    end
    if length(children(beam)) == length(_children)
        for (i, child) in enumerate(children(beam))
            _modify_beam_head!(child, _children[i])
        end
        return nothing
    end
    return error("Adding children to beam failed")
end

_drop_beams!(b::B) where {B <: AbstractBeam} = (b.children = Vector{B}())

function _modify_beam_head!(::B, ::B) where {B <: AbstractBeam}
    throw(ArgumentError(lazy"_modify_beam_head not implemented for $B"))
end

function _last_beam_intersection(::B) where {B <: AbstractBeam}
    throw(ArgumentError(lazy"_last_beam_intersection not implemented for $B"))
end

function Base.empty!(::B) where {B <: AbstractBeam}
    throw(ArgumentError(lazy"empty! not implemented for $B"))
end

"""
    first_ray(beam::AbstractBeam)

Returns the start ray on the optical axis of the `beam`; for beamlets the first chief ray.
Defines the generic `position` and `direction` of the `beam`, which are used as the pivot for rotations.
"""
function first_ray(::B) where {B <: AbstractBeam}
    throw(ArgumentError(lazy"first_ray not implemented for $B"))
end

Base.position(b::AbstractBeam) = position(first_ray(b))
direction(b::AbstractBeam) = direction(first_ray(b))

"""
    _component_beams(beam::AbstractBeam)

Returns a tuple of the component [`Beam`](@ref)s of a composite `beam` (e.g. a [`GaussianBeamlet`](@ref)),
chief beam first. Required for the kinematic API of composite beams.
"""
function _component_beams(::B) where {B <: AbstractBeam}
    throw(ArgumentError(lazy"_component_beams not implemented for $B"))
end

"""
    translate3d!(::Movable, beam::AbstractBeam, offset)

Resets the root `beam` to its untraced start state and moves it by `offset`.
Composite beams delegate to their component beams. Throws an `ArgumentError` for child beams.
"""
function translate3d!(::Movable, b::AbstractBeam, offset)
    isroot(b) || throw(ArgumentError("cannot move a child beam; move its root beam instead"))
    empty!(b)
    foreach(c -> translate3d!(c, offset), _component_beams(b))
    return nothing
end

"""
    rotate3d!(::Movable, beam::AbstractBeam, R::AbstractMatrix)

Resets the root `beam` to its untraced start state and rotates it by `R` about its `position`
(the chief ray start for beamlets). Composite beams delegate to their component beams.
Throws an `ArgumentError` for child beams.
"""
function rotate3d!(::Movable, b::AbstractBeam, R::AbstractMatrix)
    isroot(b) || throw(ArgumentError("cannot move a child beam; move its root beam instead"))
    empty!(b)
    p = position(b)
    foreach(c -> rotate3d!(c, R, p), _component_beams(b))
    return nothing
end
