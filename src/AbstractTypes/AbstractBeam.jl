"""
    AbstractBeam{T <: Real, R <: AbstractRay{T}}

A generic type for a container type which holds rays, beams etc.

# Parametrization

A subtype of `AbstractBeam` is parameterized by its main data type `T <: Real`, as well as the underlying ray representation `R <: AbstractRay{T}`.
If a beam is to be compatible with different [`AbstractRay`](@ref) implementations, it must be parameterized by `T` and `R`.
However, it can also be set to a fixed type for `T` and `R`, i.e. `MyBeam <: AbstractBeam{Float32, MyRay}`.

# Implementation reqs.

Subtypes of `AbstractBeam` must implement the following:

## Fields

- `parent`: a [`Nullable`](@ref) field that holds the same type as the subtype, used for tree navigation
- `children`: a vector that holds the same type as the subtype, used for sub-beam tracking, i.e. beamsplitting

## Getters/setters

The tree functions `parent`, `children` and `isroot` (from `AbstractTrees`) work via the `parent`/`children` fields.

## Functions

- `trace_system!(system, beam; r_max, kwargs...)`: traces the beam from its last ray through the `system`
  until no object is hit or an interaction returns `nothing`
- `push!(beam, interaction)`: appends the result of an interaction to the beam
- `_component_beams(beam)`: returns a tuple of the beams that store the rays of `beam`, chief beam first.
  A beam that stores its rays itself (e.g. [`Beam`](@ref)) returns `(beam,)`

A beam that stores its rays itself must also implement `rays` and `empty!`, otherwise the generic
methods recurse. Everything else (`rays`, `first_ray`, `empty!`, `position`, `direction`, the lengths
and the kinematic API) has generic methods built on `_component_beams`.

## Kinematic API

`AbstractBeam`s are [`BeamletOptics.Movable`](@ref) with a [`BeamletOptics.Directed`](@ref) frame, see
[`BeamletOptics.AbstractKinematicTrait`](@ref). A subtype needs no move methods of its own: the generic
`translate3d!`/`rotate3d!` move the start ray of every component beam. Rotations are applied about the
`position` of the beam, i.e. the start of `first_ray`.
[`reset_rotation3d!`](@ref) throws an `ArgumentError`, since a beam has no orientation.

!!! warning
    Only root beams can be moved. Every move resets the beam to its untraced start state via `empty!`;
    moving a child beam throws an `ArgumentError`.
"""
abstract type AbstractBeam{T <: Real, R <: AbstractRay{T}} end

kinematic_trait_of(::AbstractBeam) = Movable(Directed())

AbstractTrees.NodeType(::Type{T}) where {T <: AbstractBeam} = HasNodeType()
AbstractTrees.nodetype(::Type{T}) where {T <: AbstractBeam} = T

AbstractTrees.ParentLinks(::Type{<:AbstractBeam}) = AbstractTrees.StoredParents()
AbstractTrees.parent(beam::AbstractBeam) = beam.parent
parent!(beam::B, parent::B) where {B <: AbstractBeam} = (beam.parent = parent)

AbstractTrees.children(beam::AbstractBeam) = beam.children

AbstractTrees.printnode(io::IO, node::B; kw...) where {B <: AbstractBeam} = show(io, B)

"""
    children!(beam::B, child::B) where {B<:AbstractBeam}

    children!(beam::B, children::AbstractVector{B}) where {B<:AbstractBeam}

Attaches the `child` or the `children` to `beam` and links `beam` as their parent, e.g. the
transmitted and the reflected beam at a beamsplitter. The children continue `beam`: its length,
optical path and reference plane are counted on, see [`relaunch!`](@ref BeamletOptics.relaunch!)
for beams that start their own reference.

A beam is split once, at its end: a `beam` that already has children throws an `ArgumentError`.
"""
children!(beam::B, child::B) where {B <: AbstractBeam} = children!(beam, [child])

function children!(beam::B, _children::AbstractVector{B}) where {B <: AbstractBeam}
    isempty(children(beam)) ||
        throw(ArgumentError("the beam already has children; a beam is split once, at its end"))
    # Link parent and add children to tree
    parent!.(_children, Ref(beam))
    append!(children(beam), _children)
    return nothing
end

_drop_beams!(beam::B) where {B <: AbstractBeam} = (beam.children = Vector{B}())

"""
    relaunch!(beam::B, new_beams::AbstractVector{B}) where {B <: AbstractBeam}

Attaches `new_beams` as children of `beam` that start their own optical reference: they are
traced and removed together with `beam`, but their length, optical path and phase
are counted from their own start instead of being continued from `beam`. The new beams carry
the full phase in their amplitude.

Use it in the `interact3d` of a component that re-emits light elsewhere, e.g. the output end
of a fiber or a solver coupled through a field: the light leaves at another place than where
`beam` ended, so the geometric continuation of [`children!`](@ref BeamletOptics.children!)
(parent length, reference plane of a beamlet) does not apply.

The `new_beams` replace any children that `beam` already has. A beam that the component keeps
and attaches again is solved from its first ray, like every beam.
"""
function relaunch!(beam::B, new_beams::AbstractVector{B}) where {B <: AbstractBeam}
    _drop_beams!(beam)
    append!(children(beam), new_beams)
    return nothing
end

"""
    _component_beams(beam::AbstractBeam)

Returns a tuple of the beams that store the rays of `beam`, chief beam first: `(beam,)` for a beam
that stores its rays itself (e.g. [`Beam`](@ref)), the component beams of a composite beam (e.g. a
[`GaussianBeamlet`](@ref)). The generic beam methods, including the kinematic API, are built on it.
"""
_component_beams(beam::AbstractBeam) = error(lazy"_component_beams not defined for $(typeof(beam))")

"""
    rays(beam::AbstractBeam)

Returns the rays of the `beam`; for beamlets the rays of the chief beam.
"""
rays(beam::AbstractBeam) = rays(first(_component_beams(beam)))

"""
    first_ray(beam::AbstractBeam)

Returns the start ray on the optical axis of the `beam`; for beamlets the first chief ray.
Defines the generic `position` and `direction` of the `beam`, which are used as the pivot for rotations.
"""
first_ray(beam::AbstractBeam) = first(rays(beam))

Base.position(beam::AbstractBeam) = position(first_ray(beam))
direction(beam::AbstractBeam) = direction(first_ray(beam))

"""
    empty!(beam::AbstractBeam)

Resets the `beam` to its untraced start state: resets every component beam (e.g. the `chief`,
`waist` and `divergence` beams of a [`GaussianBeamlet`](@ref)) and drops all child beams.
"""
function Base.empty!(beam::AbstractBeam)
    foreach(empty!, _component_beams(beam))
    _drop_beams!(beam)
    return beam
end

refractive_index(beam::AbstractBeam, id::Int) = refractive_index(rays(beam)[id])

function refractive_index!(beam::AbstractBeam, id::Int, n_new::Real)
    for b in _component_beams(beam)
        refractive_index!(rays(b)[id], n_new)
    end
    return nothing
end

isentering(beam::AbstractBeam, id::Int) = isentering(rays(beam)[id])

"""
    isparentbeam(beam, ray)

Tests if the given `beam` contains the `ray` as a part of its solution.
"""
function isparentbeam(beam::AbstractBeam, ray::AbstractRay)
    return any(b -> any(r -> r === ray, rays(b)), _component_beams(beam))
end

"""
    Base.length(beam::AbstractBeam)

Calculate the length of a beam up to the point of the last intersection.

!!! tip
    Use [`optical_path_length`](@ref) to get the optical path length instead.
"""
function Base.length(beam::AbstractBeam{T}) where {T}
    # Recursively get length of beam parents
    l0 = length_parent(beam)
    l = length_rays(beam)
    return l + l0
end

"""
    optical_path_length(beam::AbstractBeam)

Calculate the optical path length of the `beam`, i.e. ``\\mathrm{OPL} = n \\cdot l``.
"""
function optical_path_length(beam::AbstractBeam{T}) where {T}
    p = AbstractTrees.parent(beam)
    l0 = isnothing(p) ? zero(T) : optical_path_length(p)
    for ray in rays(beam)
        if isnothing(intersection(ray))
            break
        end

        l0 += optical_path_length(ray)
    end

    return l0
end

function length_parent(beam::AbstractBeam{T}) where {T}
    p = AbstractTrees.parent(beam)
    if isnothing(p)
        return zero(T)
    else
        return length(p)::T
    end
end

"""
    length_before(beam::AbstractBeam, id::Integer)

Returns the geometric length of the `beam`, its parents included, up to the start of ray `id`.
"""
function length_before(beam::AbstractBeam{T}, id::Integer) where {T}
    l = length_parent(beam)
    for j in 1:(id - 1)
        l += length(rays(beam)[j])
    end
    return l
end

function length_rays(beam::AbstractBeam{T}) where {T}
    l = zero(T)
    for ray in rays(beam)
        if isnothing(intersection(ray))
            break
        end
        l += length(ray)::T
    end
    return l
end

"""
    point_on_beam(beam::AbstractBeam, t::Real)

Function to find a point given a specific distance `t` along the beam. Return the ray `index` aswell.
For negative distances, assume first ray backwards.
"""
function point_on_beam(beam::AbstractBeam{T}, t::Real)::Tuple{Point3{T}, Int} where {T}
    # Initialize counter to track cumulative length
    p = AbstractTrees.parent(beam)
    if isnothing(p)
        temp = zero(T)
    else
        temp = length(p)
    end
    numEl = length(rays(beam))
    for (index, ray) in enumerate(rays(beam))
        # Catch final ray
        if index == numEl
            break
        end
        temp += length(ray)
        # If the specified distance `t` is less than the cumulative length,
        # calculate the local ray length `b` and find the point along the ray
        if t < temp
            b = temp - t
            point = position(ray) + (length(ray) - b) * direction(ray)
            return point, index
        end
    end
    # If no solution at this point assume final ray with infinite length
    ray = last(rays(beam))
    b = t - temp
    point = position(ray) + b * direction(ray)
    return point, numEl
end

"""
    reset_root!(beam::AbstractBeam)

Resets the root `beam` to its untraced start state before a move; throws an `ArgumentError` for a child beam.
"""
function reset_root!(beam::AbstractBeam)
    isroot(beam) || throw(ArgumentError("cannot move a child beam; move its root beam instead"))
    empty!(beam)
    return nothing
end

# The hit at which a beam interacts `at` its ray, see interact3d: the ray of a beam, or the segment
# `ray_id` of the chief beam of a beamlet
intersection(::AbstractBeam, ray::AbstractRay) = intersection(ray)
intersection(beam::AbstractBeam, ray_id::Int) = intersection(rays(beam)[ray_id])

"""
    translate3d!(::Movable, beam::AbstractBeam, offset)

Resets the root `beam` to its untraced start state and moves the start ray of every component beam
by `offset`. Throws an `ArgumentError` for child beams.
"""
function translate3d!(::Movable, beam::AbstractBeam, offset)
    reset_root!(beam)
    for b in _component_beams(beam)
        translate3d!(first_ray(b), offset)
    end
    return nothing
end

"""
    rotate3d!(::Movable, beam::AbstractBeam, R::AbstractMatrix)

Resets the root `beam` to its untraced start state and rotates the start ray of every component beam
by `R` about the `position` of the `beam` (the chief ray start for beamlets). Throws an `ArgumentError`
for child beams.
"""
function rotate3d!(::Movable, beam::AbstractBeam, R::AbstractMatrix)
    reset_root!(beam)
    pivot = position(beam)
    for b in _component_beams(beam)
        rotate3d!(first_ray(b), R, pivot)
    end
    return nothing
end
