"""
    Ray{T} <: AbstractRay{T}

Immutable struct to store geometrical ray information in R³.

# Fields
- `pos`: a point in R³ that describes the `Ray` origin
- `dir`: a normalized vector in R³ that describes the `Ray` propagation direction
- `λ`: wavelength in [m]
- `n`: refractive index of current medium
"""
struct Ray{T <: Real} <: AbstractRay{T}
    pos::Point3{T}
    dir::Point3{T}
    λ::T
    n::T
end

"""
    Ray(pos, dir, λ=1000e-9, n=1.0)

Constructs an immutable `Ray` where:
- `pos`: is the `Ray` origin
- `dir`: is the `Ray` direction of propagation, normalized to unit length
- `λ`: wavelength in [m] (default 1000 nm)
- `n`: refractive index (default 1.0)
"""
function Ray(pos::AbstractArray{P},
        dir::AbstractArray{D},
        λ::L = 1000e-9,
        n::N = 1.0) where {P <: Real, D <: Real, L <: Real, N <: Real}
    F = promote_type(P, D, L, N)
    if isapprox(norm(dir), 0, atol = 1e-14)
        throw(ErrorException("Direction vector to short for normalization."))
    end
    return Ray{F}(
        Point3{F}(pos),
        normalize(Point3{F}(dir)),
        F(λ),
        F(n))
end

function Ray(pos::Point3{T},
        dir::Point3{T},
        λ::Real = 1000e-9,
        n::Real = 1.0) where {T <: Real}
    F = promote_type(T, typeof(λ), typeof(n))
    return Ray{F}(
        Point3{F}(pos),
        normalize(Point3{F}(dir)),
        F(λ),
        F(n))
end

# 5-argument compatibility constructor for legacy call sites
function Ray(pos::AbstractArray,
        dir::AbstractArray,
        ::Nullable{<:AbstractIntersection},
        λ::Real = 1000e-9,
        n::Real = 1.0)
    return Ray(pos, dir, λ, n)
end

function Ray{T}(pos::AbstractArray,
        dir::AbstractArray,
        ::Nullable{<:AbstractIntersection},
        λ::Real = 1000e-9,
        n::Real = 1.0) where {T <: Real}
    return Ray{T}(Point3{T}(pos), normalize(Point3{T}(dir)), T(λ), T(n))
end