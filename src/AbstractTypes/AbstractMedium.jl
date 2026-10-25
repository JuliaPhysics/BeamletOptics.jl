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
