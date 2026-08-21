"""
    Beam{T, R <: AbstractRay{T}} <: AbstractBeam{T, R}

Stores the ray propagation segments calculated from geometrical optics when propagating through an optical system.
The `Beam` type is parametrically defined by the [`AbstractRay`](@ref) subtype that it stores.

# Fields
- `head_ray`: the initial spawned ray at the origin of the beam
- `segments`: vector of [`RaySegment`](@ref)s, representing each propagation step from boundary to boundary
- `parent`: reference to the parent beam, if any ([`Nullable`](@ref) to account for the root beam)
- `children`: vector of child beams created during branching (e.g. at beamsplitters)
"""
mutable struct Beam{T, R <: AbstractRay{T}} <: AbstractBeam{T, R}
    head_ray::R
    segments::Vector{RaySegment{T, R}}
    parent::Nullable{Beam{T, R}}
    children::Vector{Beam{T, R}}
end

segments(b::Beam) = b.segments

@inline function Base.getproperty(b::Beam, s::Symbol)
    if s === :rays
        return rays(b)
    else
        return getfield(b, s)
    end
end


"""
    rays(beam::Beam)

Returns the sequence of rays that make up the `beam`. If the beam is not yet solved, returns `[beam.head_ray]`.
"""
function rays(b::Beam{T, R}) where {T, R}
    if isempty(b.segments)
        return R[b.head_ray]
    else
        return R[seg.ray for seg in b.segments]
    end
end

function Base.push!(b::Beam{T, R}, seg::RaySegment{T, R}) where {T, R}
    if isempty(b.segments) || seg.accumulated_opl > optical_path_length(seg)
        push!(b.segments, seg)
    else
        prev_opl = last(b.segments).accumulated_opl
        adjusted = RaySegment(seg.ray, seg.intersection, seg.hit_object, prev_opl + optical_path_length(seg))
        push!(b.segments, adjusted)
    end
end
Base.push!(b::Beam{T, R}, ray::R) where {T, R} = push!(b, RaySegment(ray))


function Beam(ray::R) where {T, R <: AbstractRay{T}}
    return Beam{T, R}(ray, Vector{RaySegment{T, R}}(), nothing, Vector{Beam{T, R}}())
end

"""
    Beam(pos, dir, λ=1e-6)

Spawns a [`Beam`](@ref) at the start `pos`ition in the specified `dir`ection
with wavelength `λ = 1000 nm`.
"""
function Beam(pos::AbstractArray{P}, dir::AbstractArray{D}, λ::L=1e-6) where {P,D,L}
    T = promote_type(P,D,L)
    ray = Ray(pos, dir, λ, 1.0)
    return Beam{T, Ray{T}}(ray, Vector{RaySegment{T, Ray{T}}}(), nothing, Vector{Beam{T, Ray{T}}}())
end

function Beam(pos::AbstractArray{P}, dir::AbstractArray{D}, λ::L, n::N) where {P,D,L,N<:Real}
    T = promote_type(P,D,L,N)
    ray = Ray(pos, dir, λ, n)
    return Beam{T, Ray{T}}(ray, Vector{RaySegment{T, Ray{T}}}(), nothing, Vector{Beam{T, Ray{T}}}())
end

"""
    Beam(pos, dir, λ, E0)

Spawns a [`Beam`](@ref) of [`PolarizedRay`](@ref)s at the start `pos`ition in the specified `dir`ection
with the wavelength `λ` and electric field vector `E0`.
"""
function Beam(pos::AbstractArray{P}, dir::AbstractArray{D}, λ::L, E0::AbstractArray) where {P,D,L}
    T = promote_type(P,D,L,eltype(E0))
    ray = PolarizedRay(pos, dir, λ, 1.0, E0)
    return Beam{T, PolarizedRay{T}}(ray, Vector{RaySegment{T, PolarizedRay{T}}}(), nothing, Vector{Beam{T, PolarizedRay{T}}}())
end


"""
    BeamInteraction <: AbstractInteraction

This type is used to store the new [`AbstractRay`](@ref) resulting from an optical interaction
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
    empty!(beam.segments)
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
    l0 = length_parent(beam)
    l = length_rays(beam)
    return l + l0
end

"""
    optical_path_length(beam::Beam)

Calculate the optical path length of the `beam`, i.e. ``\\mathrm{OPL} = \\sum n_i \\cdot l_i``.
"""
function optical_path_length(beam::Beam{T}) where {T}
    p = AbstractTrees.parent(beam)
    l0 = isnothing(p) ? zero(T) : optical_path_length(p)
    return isempty(beam.segments) ? l0 : (l0 + last(beam.segments).accumulated_opl)
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
    for seg in beam.segments
        if isnothing(seg.intersection)
            break
        end
        l += length(seg)
    end
    return l
end

"""
    point_on_beam(beam::Beam, t::Real)

Function to find a point given a specific distance `t` along the beam. Returns `(point, segment_index)`.
"""
function point_on_beam(beam::B, t::Real)::Tuple{Point3{T}, Int} where {T, R <: AbstractRay{T}, B <: Beam{T, R}}
    p = AbstractTrees.parent(beam)
    temp = isnothing(p) ? zero(T) : length(p)
    if isempty(beam.segments)
        return position(beam.head_ray) + T(t) * direction(beam.head_ray), 1
    end
    numEl = length(beam.segments)
    for (index, seg) in enumerate(beam.segments)
        if index == numEl || isnothing(seg.intersection)
            break
        end
        temp += length(seg)
        if t < temp
            b = temp - t
            point = position(seg) + (length(seg) - b) * direction(seg)
            return point, index
        end
    end
    seg = last(beam.segments)
    b = t - temp
    point = position(seg) + b * direction(seg)
    return point, numEl
end

"""
    isparaxial(system, beam, threshold=π/4)

Tests if the `beam` direction angle with respect to the surface normal exceeds `threshold` at refractive interfaces.
Reflections (e.g. at mirrors or beamsplitters) are excluded from the paraxial check.
"""
function isparaxial(::AbstractSystem, beam::Beam, threshold::Real = π / 4)
    segs = beam.segments
    n_segs = length(segs)
    for i in 1:n_segs
        seg = segs[i]
        int = seg.intersection
        isnothing(int) && break
        
        if i < n_segs
            next_seg = segs[i + 1]
            n_surf = normal3d(int)
            is_refracted = (dot(direction(seg), n_surf) * dot(direction(next_seg), n_surf) > 0)
            if !is_refracted
                continue
            end
        end

        angle = angle3d(direction(seg.ray), normal3d(int))
        if angle > π / 2
            angle = π - angle
        end
        if angle > threshold
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
    if beam.head_ray === _ray
        return true
    end
    for seg in beam.segments
        if seg.ray === _ray
            return true
        end
    end
    return false
end

function Base.show(io::IO, ::MIME"text/plain", beam::Beam)
    if isempty(beam.segments)
        println(io, "Beam (unsolved):")
        println(io, "    Origin: $(position(beam.head_ray))")
        println(io, "    Dir.:   $(direction(beam.head_ray))")
        return nothing
    end
    for (i, seg) in enumerate(beam.segments)
        println(io, "Segment $i:")
        if isnothing(seg.intersection)
            println(io, "    No intersection")
            println(io, "    Pos.: $(position(seg))")
            println(io, "    Dir.: $(direction(seg))")
        else
            println(io, "    Intersection: distance=$(length(seg)), normal=$(normal3d(seg.intersection))")
            println(io, "    Pos.: $(position(seg))")
            println(io, "    Dir.: $(direction(seg))")
            println(io, "    End.: $(propagate(seg.ray, length(seg)))")
        end
    end
    return nothing
end

