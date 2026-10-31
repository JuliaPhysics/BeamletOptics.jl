"""
    AbstractMedium

A medium in which light propagates between two boundaries, e.g. the glass of a lens or the air of
a laboratory. The bodies of a system provide their medium via [`medium`](@ref), and a
[`Transition`](@ref) holds the media on both sides of a boundary.

Further media, e.g. absorbing, birefringent or graded-index ones, will extend this interface with
the beam models that support them.

# Implementation reqs.

## Fields

None required.

## Getters/setters

None required.

## Functions

- `refractive_index(medium, λ)`: the refractive index at the wavelength `λ` in [m]
"""
abstract type AbstractMedium end

refractive_index(medium::AbstractMedium, ::Real) =
    error(lazy"refractive_index not defined for $(typeof(medium))")

"""
    IsotropicMedium(n)

A homogeneous, isotropic and transparent medium of the refractive index `n`, either a number or a
[`RefractiveIndex`](@ref), i.e. a callable that returns the index for a wavelength `λ` in [m].
"""
struct IsotropicMedium{N} <: AbstractMedium
    n::N
end

refractive_index(medium::IsotropicMedium{<:Real}, ::Real) = medium.n
refractive_index(medium::IsotropicMedium, λ::Real) = medium.n(λ)

"""
    SystemMedium(system)

The medium that surrounds the objects of the `system`, of the refractive index
`refractive_index(system, λ)`, see [`AbstractSystem`](@ref).
"""
struct SystemMedium{S <: AbstractSystem} <: AbstractMedium
    system::S
end

refractive_index(medium::SystemMedium, λ::Real) = refractive_index(medium.system, λ)

"""
    medium(x)

Returns the [`AbstractMedium`](@ref) of the body `x`, or `nothing` if `x` is no body, e.g. a mirror,
a detector or a thin beamsplitter. A ray that enters a body propagates in its medium, see
[`transition`](@ref). The medium of a system is a [`SystemMedium`](@ref).

A refracting component declares its medium via this method; every
[`AbstractRefractiveOptic`](@ref) is a body of an [`IsotropicMedium`](@ref) of its refractive index.
"""
medium(::AbstractObject) = nothing
medium(system::AbstractSystem) = SystemMedium(system)

"""
    Enclosure

The bodies that enclose a ray, innermost first: a persistent list, which the rays of a beam and the
children of a beamsplitter share. A ray outside of all bodies, i.e. in the medium of the system, has
the enclosure `nothing`. The solver updates it at every boundary, see [`transition`](@ref): the ray
enters a body ([`enter`](@ref)) or leaves one ([`leave`](@ref)). Hence the medium behind a boundary
is the body around the one that the ray leaves, e.g. the water of a cuvette around a lens.

The bodies must be closed volumes that are nested or touch, not partially overlapping ones.

# Fields

- `body`: the innermost body, see [`medium`](@ref)
- `outer`: the enclosure of `body`, or `nothing`
"""
struct Enclosure
    body::AbstractObject
    outer::Union{Nothing, Enclosure}
end

"""Returns the innermost body of the `enclosure`, or `nothing` outside of all bodies."""
innermost(enclosure::Enclosure) = enclosure.body
innermost(::Nothing) = nothing

"""
    enter(enclosure, body)

Returns the `enclosure` within the `body`, the new innermost body. A `body` that the `enclosure`
holds already, or no body (`nothing`), leaves it unchanged.
"""
enter(enclosure::Nullable{Enclosure}, body::AbstractObject) =
    encloses(enclosure, body) ? enclosure : Enclosure(body, enclosure)
enter(enclosure::Nullable{Enclosure}, ::Nothing) = enclosure

"""
    leave(enclosure, body)

Returns the `enclosure` without the `body` that the ray leaves, its innermost body for nested
bodies. A `body` deeper within it is removed with a warning once, since bodies overlap then; a
`body` that it does not hold leaves it unchanged.
"""
leave(::Nothing, body) = nothing
leave(enclosure::Enclosure, ::Nothing) = enclosure

function leave(enclosure::Enclosure, body::AbstractObject)
    enclosure.body === body && return enclosure.outer
    encloses(enclosure, body) || return enclosure
    @warn "A ray leaves a body that is not its innermost one: bodies overlap" maxlog = 1
    return enclosure_without(enclosure, body)
end

enclosure_without(enclosure::Enclosure, body) = enclosure.body === body ? enclosure.outer :
    Enclosure(enclosure.body, enclosure_without(enclosure.outer, body))

"""Returns `true` if the `enclosure` holds the `body`."""
encloses(::Nothing, body) = false
encloses(enclosure::Enclosure, body) = enclosure.body === body || encloses(enclosure.outer, body)
