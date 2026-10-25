"""
    AbstractMedium

A medium in which light propagates between two boundaries, e.g. the glass of a lens or the air of
a laboratory. The bodies of a system provide their medium via [`medium`](@ref), and a
[`Transition`](@ref) holds the media on both sides of a boundary.

# Interface

- `refractive_index(medium, λ)`: the refractive index at the wavelength `λ` in [m]

Further media, e.g. absorbing, birefringent or graded-index ones, will extend this interface with
the beam models that support them.
"""
abstract type AbstractMedium end

"""
    IsotropicMedium(n)

A homogeneous, isotropic and transparent medium of the refractive index `n`, either a number or a
[`RefractiveIndex`](@ref), i.e. a callable that returns the index for a wavelength `λ` in [m].
"""
struct IsotropicMedium{N} <: AbstractMedium
    n::N
end

refractive_index(m::IsotropicMedium{<:Real}, ::Real) = m.n
refractive_index(m::IsotropicMedium, λ::Real) = m.n(λ)

"""
    SystemMedium(system)

The medium that surrounds the objects of the `system`, of the refractive index
`refractive_index(system, λ)`, see [`AbstractSystem`](@ref).
"""
struct SystemMedium{S <: AbstractSystem} <: AbstractMedium
    system::S
end

refractive_index(m::SystemMedium, λ::Real) = refractive_index(m.system, λ)

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

"""Returns the innermost body of the enclosure `e`, or `nothing` outside of all bodies."""
innermost(e::Enclosure) = e.body
innermost(::Nothing) = nothing

"""Returns the enclosure `e` within the `body`, which is the new innermost body."""
enter(e::Union{Nothing, Enclosure}, body::AbstractObject) = Enclosure(body, e)

"""
    leave(e, body)

Returns the enclosure `e` without the `body`, which the ray leaves. It is the innermost body of `e`
for nested bodies. A `body` deeper within `e` is removed as well, with a warning once, since bodies
overlap then. A `body` that `e` does not hold leaves `e` unchanged, e.g. for a ray that a component
created without an enclosure.
"""
leave(::Nothing, body) = nothing
leave(e::Enclosure, ::Nothing) = e

function leave(e::Enclosure, body::AbstractObject)
    e.body === body && return e.outer
    encloses(e, body) || return e
    @warn "A ray leaves a body that is not its innermost one: bodies overlap" maxlog = 1
    return _without(e, body)
end

_without(e::Enclosure, body) = e.body === body ? e.outer : Enclosure(e.body, _without(e.outer, body))

"""Returns `true` if the enclosure `e` holds the `body`."""
encloses(::Nothing, body) = false
encloses(e::Enclosure, body) = e.body === body || encloses(e.outer, body)
