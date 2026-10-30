"""
    Beam{T, R <: AbstractRay{T}} <: AbstractBeam{T, R}

Stores the rays that are calculated from geometric optics when propagating through an optical system.
The `Beam` type is parametrically defined by the [`AbstractRay`](@ref) subtype that it stores.

# Fields

- `rays`: vector of `AbstractRay` objects, representing the rays that make up the beam
- `parent`: reference to the parent beam, if any ([`Nullable`](@ref) to account for the root beam which has no parent)
- `children`: vector of child beams, each child beam represents a branching or bifurcation of the original beam, i.e. beam-splitting
"""
mutable struct Beam{T, R <: AbstractRay{T}} <: AbstractBeam{T, R}
    rays::Vector{R}
    parent::Nullable{Beam{T, R}}
    children::Vector{Beam{T, R}}
end

"""
    rays(beam::Beam)

Returns the vector of rays that make up the `beam`.
"""
rays(b::Beam) = b.rays

wavelength(b::Beam) = wavelength(first(rays(b)))


Base.push!(b::Beam, ray::AbstractRay) = push!(b.rays, ray)

"""
    empty!(beam::Beam)

Resets the `beam` to its untraced start state: all rays but the first are removed, the
intersection of the start ray is cleared and all child beams are dropped.
"""
function Base.empty!(b::Beam)
    deleteat!(rays(b), 2:length(rays(b)))
    empty!(first(rays(b)))
    _drop_beams!(b)
    return b
end

_component_beams(b::Beam) = (b,)

function Beam(ray::R) where {T, R <: AbstractRay{T}}
    Beam{T, R}([ray], nothing, Vector{Beam{T, R}}())
end

"""
    Beam(pos, dir, λ=1e-6)

Spawns a [`Beam`](@ref) at the start `pos`ition in the specified `dir`ection
with the wavelength `λ = 1000 nm`.
"""
function Beam(pos::AbstractArray{P}, dir::AbstractArray{D}, λ::L=1e-6) where {P,D,L}
    T = promote_type(P,D,L)
    ray = Ray(pos, dir, λ)
    return Beam{T, Ray{T}}([ray], nothing, Vector{Beam{T, Ray{T}}}())
end

"""
    Beam(pos, dir, λ, E0)

Spawns a [`Beam`](@ref) of [`PolarizedRay`](@ref)s at the start `pos`ition in the specified
`dir`ection with the wavelength `λ` and the electric field vector `E0`, which may be complex,
e.g. `[1, 0, im] / √2` for circular polarization. The number type of the beam is that of the
[`PolarizedRay`](@ref), i.e. real also for a complex `E0`.
"""
Beam(pos::AbstractArray, dir::AbstractArray, λ, E0) = Beam(PolarizedRay(pos, dir, λ, E0))

"""
    BeamInteraction <: AbstractInteraction

This type is used to store the new [`AbstractRay`](@ref) resulting from on optical interaction
between a [`Beam`](@ref) and some [`AbstractObject`](@ref).

# Fields

- `hint`: optional [`Hint`](@ref) for the solver
- `ray`: new [`AbstractRay`](@ref) resulting from the interaction
"""
mutable struct BeamInteraction{T <: Real, R <: AbstractRay{T}} <: AbstractInteraction
    hint::Nullable{Hint}
    ray::R
end

Base.push!(b::Beam, interaction::BeamInteraction) = push!(b, interaction.ray)

"""
    isparaxial(system, beam, threshold=π/4)

Tests the angle between the `beam` direction and surface normal at each intersection.
Mainly intended as a check for [`GaussianBeamlet`](@ref).
"""
function isparaxial(::AbstractSystem, beam::Beam, threshold::Real = π / 4)
    # Test if refractive elements are hit with angle larger than threshold
    for ray in beam.rays
        if isnothing(intersection(ray))
            break
        end
        target = object((intersection(ray)))
        # Test if refractive element
        if !isa(target, AbstractRefractiveOptic)
            continue
        end
        # Test angle between ray and its intersection
        angle = angle3d(ray)
        if angle > π / 2 # flip sector
            angle = π - angle
        end
        if angle > threshold # rad
            return false
        end
    end
    return true
end

function Base.show(io::IO, ::MIME"text/plain", beam::Beam)
    for (i, ray) in enumerate(rays(beam))
        println(io, "Ray $i:")
        if isnothing(intersection(ray))
            println(io, "    No intersection")
            println(io, "    Pos.: $(position(ray))")
            println(io, "    Dir.: $(direction(ray))")
        else
            println(io, "    Intersects with $(typeof(object(intersection(ray))))")
            println(io, "    Pos.: $(position(ray))")
            println(io, "    Dir.: $(direction(ray))")
            println(io, "    End.: $(hit_point(ray))")
        end
    end
    return nothing
end
