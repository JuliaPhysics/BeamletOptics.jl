"""
    Beam

Stores the rays that are calculated from geometric optics when propagating through an optical system.

# Fields
- `id`: beam ID (uuid4)
- `rays`: vector of `Ray` objects, representing the rays that make up the beam
- `parent`: reference to the parent beam, if any ([`Nullable`](@ref) to account for the root beam which has no parent)
- `children`: vector of child beams, each child beam represents a branching or bifurcation of the original beam, i.e. beam-splitting
"""
mutable struct Beam{T} <: AbstractBeam{T}
    id::UUID
    rays::Vector{Ray{T}}
    parent::Nullable{Beam{T}}
    children::Vector{Beam{T}}
end

rays(b::Beam) = b.rays

Base.push!(b::Beam, ray::AbstractRay) = push!(b.rays, ray)

Beam(ray::Ray{T}) where {T} = Beam{T}(uuid4(), [ray], nothing, Vector{Beam{T}}())

struct BeamInteraction{R <: Real} <: AbstractInteraction
    id::Nullable{UUID}
    ray::Ray{R}
end

Base.push!(b::Beam, interaction::BeamInteraction) = push!(b, interaction.ray)

function Base.replace!(beam::Beam{T}, interaction::BeamInteraction{T}, index::Int) where {T}
    position!(rays(beam)[index], position(interaction.ray))
    direction!(rays(beam)[index], direction(interaction.ray))
    parameters!(rays(beam)[index], parameters(interaction.ray))
end

function _modify_beam_head!(old::Beam{T}, new::Beam{T}) where {T <: Real}
    position!(first(rays(old)), position(first(rays(new))))
    direction!(first(rays(old)), direction(first(rays(new))))
    parameters!(first(rays(old)), parameters(first(rays(new))))
    return nothing
end

_last_beam_intersection(beam::Beam) = intersection(last(rays(beam)))

"""
    length(beam::Beam)

Calculate the length of a beam up to the point of the last intersection.
"""
function Base.length(beam::Beam{T}) where {T}
    # Recursively get length of beam parents
    p = AbstractTrees.parent(beam)
    l0 = zero(T)
    if !isnothing(p)
        l0 = length(p)
    end
    # Calculate ray lengths
    l = zero(T)
    for ray in rays(beam)
        if isnothing(intersection(ray))
            break
        end
        l += length(ray)
    end
    return l + l0
end

"""
    point_on_beam(beam::Beam, t::Real)

Function to find a point given a specific distance `t` along the beam. Return the ray `index` aswell.
For negative distances, assume first ray backwards.
"""
function point_on_beam(beam::Beam{T}, t::Real) where {T}
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
function isparaxial(system::AbstractSystem, beam::Beam, threshold::Real = π / 4)
    # Test if refractive elements are hit with angle larger than threshold
    for ray in beam.rays
        if isnothing(intersection(ray))
            break
        end
        target = object(system, id(intersection(ray)))
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

function isparentbeam(beam::Beam, ray_id::UUID)
    for ray in rays(beam)
        if id(ray) == ray_id
            return true
        end
    end
    return false
end

"""
    isparentbeam(beam, ray)

Tests if the given `beam` contains the `ray` as a part of its solution.
"""
isparentbeam(beam::Beam, ray::AbstractRay) = isparentbeam(beam, ray.id)

function Base.show(io::IO, ::MIME"text/plain", beam::Beam)
    for (i, ray) in enumerate(rays(beam))
        println(io, "Ray $i:")
        if isnothing(intersection(ray))
            println(io, "    No intersection")
            println(io, "    Pos.: $(position(ray))")
            println(io, "    Dir.: $(direction(ray))")
        else
            println(io, "    Intersects with object #$(intersection(ray).id)")
            println(io, "    Pos.: $(position(ray))")
            println(io, "    Dir.: $(direction(ray))")
            println(io, "    End.: $(position(ray) .+ length(ray) .* direction(ray))")
        end
    end
    return nothing
end
