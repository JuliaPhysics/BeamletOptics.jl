"""
    Beam{T, R <: AbstractRay{T}} <: AbstractBeam{T, R}

Stores the ray propagation sequence calculated from geometrical optics when propagating through an optical system.
The `Beam` type is parametrically defined by the [`AbstractRay`](@ref) subtype that it stores.

# Fields
- `rays`: vector of [`AbstractRay`](@ref)s, representing each propagation step from boundary to boundary
- `parent`: reference to the parent beam, if any ([`Nullable`](@ref) to account for the root beam)
- `children`: vector of child beams created during branching (e.g. at beamsplitters)
"""
mutable struct Beam{T, R <: AbstractRay{T}} <: AbstractBeam{T, R}
    rays::Vector{R}
    parent::Nullable{Beam{T, R}}
    children::Vector{Beam{T, R}}
end

rays(b::Beam) = b.rays

function Base.push!(b::Beam{T, R}, ray::AbstractRay{T}) where {T, R}
    if isempty(b.rays) || accumulated_opl(ray) > optical_path_length(ray)
        push!(b.rays, ray)
    else
        prev_opl = accumulated_opl(last(b.rays))
        new_opl = prev_opl + optical_path_length(ray)
        adjusted_ray = with_accumulated_opl(ray, new_opl)
        push!(b.rays, adjusted_ray)
    end
end

function Beam(ray::Ray{T}) where {T}
    return Beam{T, Ray{T}}(Ray{T}[ray], nothing, Vector{Beam{T, Ray{T}}}())
end

function Beam(ray::PolarizedRay{T}) where {T}
    return Beam{T, PolarizedRay{T}}(PolarizedRay{T}[ray], nothing, Vector{Beam{T, PolarizedRay{T}}}())
end

function Beam{T, R}(ray::R) where {T, R <: AbstractRay{T}}
    return Beam{T, R}(R[ray], nothing, Vector{Beam{T, R}}())
end

"""
    Beam(pos, dir, λ=1e-6)

Spawns a [`Beam`](@ref) at the start `pos`ition in the specified `dir`ection
with wavelength `λ = 1000 nm`.
"""
function Beam(pos::AbstractArray{P}, dir::AbstractArray{D}, λ::L=1e-6) where {P,D,L}
    T = promote_type(P,D,L)
    ray = Ray(pos, dir, λ, 1.0)
    return Beam{T, Ray{T}}(Ray{T}[ray], nothing, Vector{Beam{T, Ray{T}}}())
end

function Beam(pos::AbstractArray{P}, dir::AbstractArray{D}, λ::L, n::N) where {P,D,L,N<:Real}
    T = promote_type(P,D,L,N)
    ray = Ray(pos, dir, λ, n)
    return Beam{T, Ray{T}}(Ray{T}[ray], nothing, Vector{Beam{T, Ray{T}}}())
end

"""
    Beam(pos, dir, λ, E0)

Spawns a [`Beam`](@ref) of [`PolarizedRay`](@ref)s at the start `pos`ition in the specified `dir`ection
with the wavelength `λ` and electric field vector `E0`.
"""
function Beam(pos::AbstractArray{P}, dir::AbstractArray{D}, λ::L, E0::AbstractArray) where {P,D,L}
    T = promote_type(P,D,L,eltype(E0))
    ray = PolarizedRay(pos, dir, λ, 1.0, E0)
    return Beam{T, PolarizedRay{T}}(PolarizedRay{T}[ray], nothing, Vector{Beam{T, PolarizedRay{T}}}())
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
    if !isempty(beam.rays)
        first_ray = first(beam.rays)
        reset_seg = OpenSegment(position(first_ray), direction(first_ray), 0.0)
        reset_ray = with_segment(first_ray, reset_seg)
        empty!(beam.rays)
        push!(beam.rays, reset_ray)
    end
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
    return isempty(beam.rays) ? l0 : (l0 + accumulated_opl(last(beam.rays)))
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
    for ray in beam.rays
        if isnothing(intersection(ray))
            break
        end
        l += length(ray)
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
    if isempty(beam.rays)
        error("Cannot query point on empty beam")
    end
    numEl = length(beam.rays)
    for (index, ray) in enumerate(beam.rays)
        if index == numEl || isnothing(intersection(ray))
            break
        end
        temp += length(ray)
        if t < temp
            b = temp - t
            point = position(ray) + (length(ray) - b) * direction(ray)
            return point, index
        end
    end
    ray = last(beam.rays)
    b = t - temp
    point = position(ray) + b * direction(ray)
    return point, numEl
end

"""
    isparaxial(system, beam, threshold=π/4)

Tests if the `beam` direction angle with respect to the surface normal exceeds `threshold` at refractive interfaces.
Reflections (e.g. at mirrors or beamsplitters) are excluded from the paraxial check.
"""
function isparaxial(::AbstractSystem, beam::Beam, threshold::Real = π / 4)
    rays_vec = beam.rays
    n_rays = length(rays_vec)
    for i in 1:n_rays
        ray = rays_vec[i]
        int = intersection(ray)
        isnothing(int) && break
        
        if i < n_rays
            next_ray = rays_vec[i + 1]
            n_surf = normal3d(int)
            is_refracted = (dot(direction(ray), n_surf) * dot(direction(next_ray), n_surf) > 0)
            if !is_refracted
                continue
            end
        end

        angle = angle3d(direction(ray), normal3d(int))
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
    for ray in beam.rays
        if ray === _ray
            return true
        end
    end
    return false
end

function Base.show(io::IO, ::MIME"text/plain", beam::Beam)
    if isempty(beam.rays)
        println(io, "Beam (empty)")
        return nothing
    end
    if length(beam.rays) == 1 && isnothing(intersection(first(beam.rays)))
        println(io, "Beam (unsolved):")
        println(io, "    Origin: $(position(first(beam.rays)))")
        println(io, "    Dir.:   $(direction(first(beam.rays)))")
        return nothing
    end
    for (i, ray) in enumerate(beam.rays)
        println(io, "Segment $i:")
        if isnothing(intersection(ray))
            println(io, "    No intersection")
            println(io, "    Pos.: $(position(ray))")
            println(io, "    Dir.: $(direction(ray))")
        else
            println(io, "    Intersection: distance=$(length(ray)), normal=$(normal3d(ray))")
            println(io, "    Pos.: $(position(ray))")
            println(io, "    Dir.: $(direction(ray))")
            println(io, "    End.: $(propagate(ray, length(ray)))")
        end
    end
    return nothing
end
