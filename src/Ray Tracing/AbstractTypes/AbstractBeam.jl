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
- `children`: a vector that holds the same type as the subtype, used for sub-beam tracking, i.e. beam splitting

## Functions:

- `_modify_beam_head!`: modifies the beam path for retracing purposes
- `_last_beam_intersection`: returns the last `Beam` intersection
"""
abstract type AbstractBeam{T <: Real, R <: AbstractRay{T}} end

AbstractTrees.NodeType(::Type{T}) where {T <: AbstractBeam} = HasNodeType()
AbstractTrees.nodetype(::Type{T}) where {T <: AbstractBeam} = T

AbstractTrees.ParentLinks(::Type{<:AbstractBeam}) = AbstractTrees.StoredParents()
AbstractTrees.parent(beam::AbstractBeam) = beam.parent
parent!(beam::B, parent::B) where {B <: AbstractBeam} = (beam.parent = parent)

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
