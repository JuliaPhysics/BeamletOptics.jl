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

rays(b::Beam) = b.rays

Base.push!(b::Beam, ray::AbstractRay) = push!(b.rays, ray)

function Beam(ray::R) where {T, R <: AbstractRay{T}}
    Beam{T, R}([ray], nothing, Vector{Beam{T, R}}())
end

function Beam(pos::AbstractArray{P}, dir::AbstractArray{D}, λ::L) where {P,D,L}
    T = promote_type(P,D,L)
    ray = Ray(pos, dir, λ)
    return Beam{T, Ray{T}}([ray], nothing, Vector{Beam{T, Ray{T}}}())
end

function Beam(pos::AbstractArray{P}, dir::AbstractArray{D}, λ::L, E0::Vector{E}) where {P,D,L,E}
    T = promote_type(P,D,L,E)
    ray = PolarizedRay(pos, dir, λ, E0)
    return Beam{T, PolarizedRay{T}}([ray], nothing, Vector{Beam{T, PolarizedRay{T}}}())
end

struct BeamInteraction{T <: Real, R <: AbstractRay{T}} <: AbstractInteraction
    hint::Nullable{Hint}
    ray::R
end

Base.push!(b::Beam, interaction::BeamInteraction) = push!(b, interaction.ray)

function Base.replace!(beam::Beam{<:Real, <:Ray}, interaction::BeamInteraction{<:Real, <:Ray}, index::Int)
    position!(rays(beam)[index], position(interaction.ray))
    direction!(rays(beam)[index], direction(interaction.ray))
    wavelength!(rays(beam)[index], wavelength(interaction.ray))
    refractive_index!(rays(beam)[index], refractive_index(interaction.ray))
end

function Base.replace!(beam::Beam{<:Real, <:PolarizedRay}, interaction::BeamInteraction{<:Real, <:PolarizedRay}, index::Int)
    position!(rays(beam)[index], position(interaction.ray))
    direction!(rays(beam)[index], direction(interaction.ray))
    wavelength!(rays(beam)[index], wavelength(interaction.ray))
    refractive_index!(rays(beam)[index], refractive_index(interaction.ray))
    polarization!(rays(beam)[index], polarization(interaction.ray))
end

function _modify_beam_head!(old::Beam{T, R}, new::Beam{T, R}) where {T<:Real, R<:Ray{T}}
    position!(first(rays(old)), position(first(rays(new))))
    direction!(first(rays(old)), direction(first(rays(new))))
    wavelength!(first(rays(old)), wavelength(first(rays(new))))
    refractive_index!(first(rays(old)), refractive_index(first(rays(new))))
    return nothing
end

function _modify_beam_head!(old::Beam{T, R}, new::Beam{T, R}) where {T<:Real, R<:PolarizedRay{T}}
    position!(first(rays(old)), position(first(rays(new))))
    direction!(first(rays(old)), direction(first(rays(new))))
    wavelength!(first(rays(old)), wavelength(first(rays(new))))
    refractive_index!(first(rays(old)), refractive_index(first(rays(new))))
    polarization!(first(rays(old)), polarization(first(rays(new))))
    return nothing
end

_last_beam_intersection(beam::Beam) = intersection(last(rays(beam)))

"""
    length(beam::Beam; opl::Bool=false)

Calculate the length of a beam up to the point of the last intersection.
The `opl` keyword can be used to calculate the `o`ptical `p`ath `l`ength instead, i.e. ``OPL = \\sum_{i=0}^{k} n_i \\cdot l_i``.
Default is the geometrical length.
"""
function Base.length(beam::Beam{T}; opl::Bool=false) where {T}
    # Recursively get length of beam parents
    p = AbstractTrees.parent(beam)
    l0 = zero(T)
    if !isnothing(p)
        l0 = length(p; opl)
    end
    # Calculate ray lengths
    l = zero(T)
    for ray in rays(beam)
        if isnothing(intersection(ray))
            break
        end
        l += length(ray; opl)
    end
    return l + l0
end

"""
    point_on_beam(beam::Beam, t::Real)

Function to find a point given a specific distance `t` along the beam. Return the ray `index` aswell.
For negative distances, assume first ray backwards.
"""
function point_on_beam(beam::Beam{T}, t::Real)::Tuple{Point3{T}, Int} where {T}
    # Initialize counter to track cumulative length
    p = AbstractTrees.parent(beam)
    if isnothing(p)
        temp = zero(t)
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
