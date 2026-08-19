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


Base.push!(b::Beam, ray::AbstractRay) = push!(b.rays, ray)

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

Spawns a [`Beam`](@ref) at the start `pos`ition in the specified `dir`ection
with the wavelength `λ` and field vector `E0`.
"""
function Beam(pos::AbstractArray{P}, dir::AbstractArray{D}, λ::L, E0) where {P,D,L}
    T = promote_type(P,D,L)
    ray = PolarizedRay(pos, dir, λ, E0)
    return Beam{T, PolarizedRay{T}}([ray], nothing, Vector{Beam{T, PolarizedRay{T}}}())
end

"""
    Beam(pos, dir, λ, E0)

Spawns a [`Beam`](@ref) of [`PolarizedRay`](@ref)s at the start `pos`ition in the specified `dir`ection
with the wavelength `λ` and electric field vector `E0`
"""
function Beam(pos::AbstractArray{P}, dir::AbstractArray{D}, λ::L, E0::Vector{E}) where {P,D,L,E}
    T = promote_type(P,D,L,E)
    ray = PolarizedRay(pos, dir, λ, E0)
    return Beam{T, PolarizedRay{T}}([ray], nothing, Vector{Beam{T, PolarizedRay{T}}}())
end

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

function _reset_beam!(beam::Beam)
    deleteat!(rays(beam), 2:length(rays(beam)))
    intersection!(first(rays(beam)), nothing)
    _drop_beams!(beam)
    return nothing
end

"""
    Base.length(beam::Beam)

Calculate the length of a beam up to the point of the last intersection.

!!! tip
    Use [`optical_path_length`](@ref) to get the optical path length instead.
"""
function Base.length(beam::Beam{T}) where {T}
    # Recursively get length of beam parents
    l0 = length_parent(beam)
    l = length_rays(beam)
    return l + l0
end

"""
    optical_path_length(beam::Beam)

Calculate the optical path length of the `beam`, i.e. ``\\mathrm{OPL} = n \\cdot l``.
"""
function optical_path_length(beam::Beam{T}) where {T}
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

function length_parent(beam::Beam{T}) where {T}
    p = AbstractTrees.parent(beam)
    if isnothing(p)
        return zero(T)
    else
        return length(p)::T
    end
end

function length_rays(beam::Beam{T}) where {T}
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
    point_on_beam(beam::Beam, t::Real)

Function to find a point given a specific distance `t` along the beam. Return the ray `index` aswell.
For negative distances, assume first ray backwards.
"""
function point_on_beam(beam::B, t::Real)::Tuple{Point3{T}, Int} where {T, R <: AbstractRay{T}, B <: Beam{T, R}}
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

"""
    isparentbeam(beam, ray)

Tests if the given `beam` contains the `ray` as a part of its solution.
"""
function isparentbeam(beam::Beam, _ray::AbstractRay)
    for ray in rays(beam)
        if ray === _ray
            return true
        end
    end
    return false
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
            println(io, "    End.: $(position(ray) .+ length(ray) .* direction(ray))")
        end
    end
    return nothing
end
