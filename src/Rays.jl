"""
    Ray{T, S} <: AbstractRay{T, S}

Immutable struct to store geometrical ray information and physical parameters in R³.

# Fields
- `segment::S`: the trajectory segment of the ray (subtyping [`AbstractSegment`](@ref))
- `λ::T`: wavelength in [m]
- `n::T`: refractive index of current medium
"""
struct Ray{T <: Real, S <: AbstractSegment{T}} <: AbstractRay{T, S}
    segment::S
    λ::T
    n::T
end

"""
    Ray(segment, λ=1000e-9, n=1.0)

Constructs a `Ray` from an existing `segment`.
"""
function Ray(seg::S, λ::L = 1000e-9, n::N = 1.0) where {T <: Real, S <: AbstractSegment{T}, L <: Real, N <: Real}
    F = promote_type(T, L, N)
    return Ray{F, S}(seg, F(λ), F(n))
end

"""
    Ray(ray::Ray, segment)

Resolves or replaces the trajectory segment of an existing `ray` with a new `segment`.
"""
Ray{T, S}(ray::Ray{T}, seg::S) where {T <: Real, S <: AbstractSegment{T}} = Ray{T, S}(seg, ray.λ, ray.n)
Ray(ray::Ray{T}, seg::S) where {T <: Real, S <: AbstractSegment{T}} = Ray{T, S}(seg, ray.λ, ray.n)
with_segment(ray::Ray, seg::AbstractSegment) = Ray(ray, seg)

"""
    Ray(pos, dir, λ=1000e-9, n=1.0)

Constructs an immutable `Ray` with an [`OpenSegment`](@ref) where:
- `pos`: is the `Ray` origin
- `dir`: is the `Ray` direction of propagation, normalized to unit length
- `λ`: wavelength in \\[m\\] (default 1000 nm)
- `n`: refractive index (default 1.0)
"""
function (::Type{Ray{T}})(pos::AbstractArray{P},
        dir::AbstractArray{D},
        λ::L = 1000e-9,
        n::N = 1.0) where {T <: Real, P <: Real, D <: Real, L <: Real, N <: Real}
    seg = OpenSegment(pos, dir, zero(T))
    return Ray{T, typeof(seg)}(seg, T(λ), T(n))
end

function Ray(pos::AbstractArray{P},
        dir::AbstractArray{D},
        λ::L = 1000e-9,
        n::N = 1.0) where {P <: Real, D <: Real, L <: Real, N <: Real}
    F = promote_type(P, D, L, N)
    return Ray{F}(pos, dir, F(λ), F(n))
end

function Ray(pos::Point3{T},
        dir::Point3{T},
        λ::Real = 1000e-9,
        n::Real = 1.0) where {T <: Real}
    F = promote_type(T, typeof(λ), typeof(n))
    return Ray{F}(pos, dir, F(λ), F(n))
end