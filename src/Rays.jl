"""
    Ray{T} <: AbstractRay{T}

Mutable struct to store ray information.

# Fields

- `pos`: a point in R³ that describes the `Ray` origin
- `dir`: a normalized vector in R³ that describes the `Ray` direction
- `intersection`: refer to [`Intersection`](@ref)
- `λ`: wavelength in [m]
- `n`: refractive index along the beam path
"""
mutable struct Ray{T} <: AbstractRay{T}
    pos::Point3{T}
    dir::Point3{T}
    intersection::Nullable{Intersection{T}}
    λ::T
    n::T
end

"""
    Ray(pos, dir, λ=1000e-9)

Constructs a `Ray` where:

- `pos`: is the `Ray` origin
- `dir`: is the `Ray` direction of propagation, normalized to unit length

Optionally, a wavelength `λ` can be specified. The start refractive index is assumed to be in vacuum (n = 1).
"""
function Ray(pos::AbstractArray{P},
        dir::AbstractArray{D},
        λ::L = 1000e-9) where {P <: Real, D <: Real, L<:Real}
    F = promote_type(P, D, L)
    dir = normalize(dir)
    return Ray{F}(
        Point3{F}(pos),
        Point3{F}(dir),
        nothing,
        F(λ),
        F(1))
end