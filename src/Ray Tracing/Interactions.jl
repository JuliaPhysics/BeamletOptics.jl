"""
    Interaction{T}

Describes how a `ray` and an `object` interact with each other.
Returns the new parameters for the construction of the next ray.
"""
struct Interaction{T}
    pos::Vector{T}
    dir::Vector{T}
    parameters::Parameters{T}
    hint::Nullable{UUID}
end

position(interaction::Interaction) = interaction.pos

direction(interaction::Interaction) = interaction.dir

parameters(interaction::Interaction) = interaction.parameters

hint(interaction::Interaction) = interaction.hint

abstract type AbstractMirror <: AbstractObject end

struct Mirror{S <: AbstractShape} <: AbstractMirror
    id::UUID
    shape::S
end

function interact3d(::AbstractMirror, ray::AbstractRay{R}) where {R<:Real}
    normal = intersection(ray).n
    npos = position(ray) + length(ray) * direction(ray)
    ndir = reflection3d(direction(ray), normal)
    return Interaction{R}(npos, ndir, ray.parameters, nothing)
end

abstract type AbstractPrism <: AbstractObject end

struct Prism{S <: AbstractShape} <: AbstractPrism
    id::UUID
    shape::S
    ref_index::Function
end

refractive_index(prism::Prism) = prism.ref_index

function interact3d(prism::Prism, ray::Ray{R}) where {R}
    # Check dir. of ray and surface normal
    normal = intersection(ray).n
    λ = wavelength(ray)
    if dot(direction(ray), normal) < 0
        # "Outside prism"
        n1 = 1.0
        n2 = refractive_index(prism)(λ)
    else
        # "Inside prism"
        n1 = refractive_index(prism)(λ)
        n2 = 1.0
        normal *= -1
    end
    # Calculate new dir. and pos.
    ndir = refraction3d(direction(ray), normal, n1, n2)
    npos = position(ray) + length(ray) * direction(ray)
    # Hint is the current object ID
    return Interaction{R}(npos, ndir, Parameters(λ, n2), id(prism))
end

abstract type AbstractDetector <: AbstractObject end

mutable struct Photodetector{S <: AbstractShape, T} <: AbstractDetector
    const id::UUID
    const shape::S
    grid::Matrix{T}
    beams::Nullable{Vector{UUID}}
end

function interact3d(pd::Photodetector, ray::Ray)
    if isnothing(pd.beams)
        pd.beams = [id(ray)]
    else
        append!(pd.beams, [id(ray)])
    end
    return nothing
end