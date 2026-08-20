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

- `_reset_beam!`: resets the beam tree back to just its head ray(s), ready for a fresh trace
"""
abstract type AbstractBeam{T <: Real, R <: AbstractRay{T}} end

AbstractTrees.NodeType(::Type{T}) where {T <: AbstractBeam} = HasNodeType()
AbstractTrees.nodetype(::Type{T}) where {T <: AbstractBeam} = T

AbstractTrees.ParentLinks(::Type{<:AbstractBeam}) = AbstractTrees.StoredParents()
AbstractTrees.parent(beam::AbstractBeam) = beam.parent
function parent!(beam::B, parent::B) where {B <: AbstractBeam}
    if beam !== parent
        beam.parent = parent
    end
    return nothing
end

AbstractTrees.children(b::AbstractBeam) = b.children

AbstractTrees.printnode(io::IO, node::B; kw...) where {B <: AbstractBeam} = show(io, B)

"""
    children!(beam::B, child::B) where {B<:AbstractBeam}

Handles the inclusion of adding a single `child` to an existing `beam`. The function behaves as follows:

1. If no previous children exist, add child
2. Else throw error
"""
function children!(beam::B, child::B) where {B <: AbstractBeam}
    if isempty(children(beam))
        # Link parent and add child to tree
        parent!(child, beam)
        push!(children(beam), child)
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
    return error("Adding children to beam failed")
end

function children!(beam::B, _children::Tuple{Vararg{B}}) where {B <: AbstractBeam}
    if isempty(children(beam))
        # Link parent and add children to tree
        for child in _children
            parent!(child, beam)
            push!(children(beam), child)
        end
        return nothing
    end
    return error("Adding children to beam failed")
end

_drop_beams!(b::B) where {B <: AbstractBeam} = (b.children = Vector{B}())

"""
    _reset_beam!(beam::AbstractBeam)

Resets `beam` back to just its head ray(s) — drops all children, truncates ray vector(s)
to the first ray, and clears that ray's intersection — so that [`solve_system!`](@ref)
can resolve it from scratch.
"""
function _reset_beam!(::B) where {B <: AbstractBeam}
    throw(ArgumentError(lazy"_reset_beam! not implemented for $B"))
end

"""
    optical_power(beam::AbstractBeam)

Returns the optical power of the beam. Default fallback returns `1.0` (representing 100% normalized power).
"""
optical_power(::AbstractBeam{T, R}) where {T <: Real, R} = one(T)

"""
    AbstractBeamGroup

Provides a generic container type interface for bundles of [`Beam`](@ref)s.
This interface assumes that there exists a central beam around which the bundle propagates,
e.g. akin to an optical axis.

# Implementation reqs.

Subtypes of `AbstractBeamGroup` must implement the following:

## Fields:

- `beams`: a vector or tuple of [`Beam`](@ref)s

## Functions:

If the `beams` field does not exist, the following getters must be dispatched:

- `beams`: getter for the `beams` field or equivalent return type
- `position`: getter for the starting position of the central [`Beam`](@ref)
- `direction`: getter for the starting direction of the central [`Beam`](@ref)
- `wavelength`: getter for the common wavelength of the beam bundle
"""
abstract type AbstractBeamGroup{T <: Real, R <: AbstractRay{T}} end

beams(bg::AbstractBeamGroup) = bg.beams

Base.length(bg::AbstractBeamGroup) = length(beams(bg))
Base.iterate(bg::AbstractBeamGroup, state...) = iterate(beams(bg), state...)
Base.getindex(bg::AbstractBeamGroup, i::Int) = getindex(beams(bg), i)

Base.position(bg::AbstractBeamGroup) = position(first(rays(first(beams(bg)))))
direction(bg::AbstractBeamGroup) = direction(first(rays(first(beams(bg)))))

wavelength(bg::AbstractBeamGroup) = wavelength(first(rays(first(beams(bg)))))

function Base.show(io::IO, ::MIME"text/plain", bg::AbstractBeamGroup)
    println(io, "Subtype of AbstractBeamGroup")
    println(io, "   # of beams: $(length(beams(bg)))")
    return nothing
end
