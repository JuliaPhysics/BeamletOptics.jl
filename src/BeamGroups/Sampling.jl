#=
Sampling of the rays of the sources: the reference vector of the sampling pattern, the sampling
types stored in a source and the beams they generate, see `source_beams`
=#

"""
    sampling_basis(dir, basis, T)

Returns the unit-length reference vector that seeds the azimuthal sampling of a beam group.
`dir` must already be normalized.

If `basis` is `nothing`, a vector normal to `dir` is picked deterministically via
[`normal3d`](@ref). Otherwise `basis` is projected into the plane normal to `dir`, which
lets a caller rotate the sampling pattern of a source about its own axis. For a tuple of grid
axes `(ex, ey)`, the first axis is used.

Throws an `ArgumentError` if `basis` has no significant component in that plane, i.e. if it is
zero or parallel to `dir`. The test is relative to `norm(basis)` because projecting an
unnormalized vector leaves a residual that scales with its length; an absolute
tolerance would let a parallel `basis` through and yield `NaN` sampling.
"""
sampling_basis(dir::AbstractVector, ::Nothing, ::Type{<:Real}) = normal3d(dir)

function sampling_basis(dir::AbstractVector, basis::AbstractVector, ::Type{T}) where {T <: Real}
    b1 = basis - dot(basis, dir) * dir
    if norm(b1) ≤ sqrt(eps(T)) * norm(basis)
        throw(ArgumentError("Source `basis` must not be zero or parallel to `dir`"))
    end
    return normalize(b1)
end

sampling_basis(dir::AbstractVector, basis::Tuple{AbstractVector, AbstractVector}, ::Type{T}) where {T <: Real} =
    sampling_basis(dir, first(basis), T)

"""
    rays_per_ring(weights, num_rays) -> Vector{Int}

The numbers of rays on the rings of a ring sampling, `num_rays` in total, in proportion to the
`weights` (circumferences) of the rings, half of it for the outermost ring (trapezoidal rule).
"""
function rays_per_ring(weights::AbstractVector{<:Real}, num_rays::Integer)
    w = collect(float(eltype(weights)), weights)
    w[end] /= 2
    total = sum(w)
    ds = total / num_rays
    counts = round.(Int, w / ds)
    # Assign the rounding remainder to the outermost ring
    counts[end] += num_rays - sum(counts)
    return counts
end

"""
    cone_ring_directions(dir, b1, θ, num_rings, num_rays) -> Vector

The `num_rays` directions on the cones around the unit vector `dir` of a ring sampling with the
half spread angle `θ` [rad], without the center: `num_rings - 1` cones, the outermost at `θ`,
each starting at the azimuth of `b1`. The rays per cone follow `rays_per_ring`. The angles are
computed in the precision of `dir` and `θ`.
"""
function cone_ring_directions(dir::AbstractVector, b1::AbstractVector, θ::Real, num_rings::Integer, num_rays::Integer)
    T = promote_type(float(eltype(dir)), float(typeof(θ)))
    b2 = normal3d(dir, b1)
    θ_NA = LinRange(0, θ, num_rings)
    ndirs = [rotate3d(b2, step(θ_NA) * i) * dir for i in eachindex(θ_NA[2:end])]
    counts = rays_per_ring(norm.(ndirs .- dot.(ndirs, Ref(dir)) .* Ref(dir)) .* Ref(2π), num_rays)
    dirs = Vector{eltype(ndirs)}()
    for (ndir, n) in zip(ndirs, counts)
        iszero(n) && continue
        R = rotate3d(dir, 2 * T(π) / n)
        d = ndir
        for _ in 1:n
            push!(dirs, d)
            d = R * d
        end
    end
    return dirs
end

"""
    golden_angle

The golden angle `π(3 - √5) ≈ 2.39996` rad, the azimuthal increment of the sunflower samplings.
"""
const golden_angle = π * (3 - √5)

"""
    AbstractSourceSampling

How the rays of a [`CollimatedSource`](@ref) or [`PointSource`](@ref) are sampled. The source
stores its sampling in the field `sampling`, such that [`num_rays!`](@ref) can regenerate its
rays with another number of rays. A source built from given beams stores
[`BeamletOptics.NoSampling`](@ref). The generic constructors
`CollimatedSource(pos, dir, diameter, sampling, λ; num_rays)` and
`PointSource(pos, dir, sampling, λ; num_rays)` build a source from any sampling; the catalog
constructors call them.

# Implementation reqs.

Parameters of the sampling other than the number of rays (e.g. the number of rings or the
half spread angle) are fields of the subtype.

## Functions

A subtype `S <: AbstractSourceSampling` implements:

- `source_beams(s::S, pos, dir, b1, args..., λ, num_rays::Int, T) -> Vector{<:Beam}`: the beams
  of the source, see [`BeamletOptics.source_beams`](@ref). `args` are the arguments of the source
  type, i.e. the `diameter` [m] for a `CollimatedSource` and none for a `PointSource`. The
  sampling pattern must depend on the pose only via `pos`, `dir` and `b1`, i.e. it must move
  rigidly with the source. The first beam must be the center ray, which starts at `pos` along
  `dir`. Throws an `ArgumentError` if `num_rays` is below `min_num_rays(s)`.
- `min_num_rays(s::S) -> Union{Nothing, Int}` (optional): the fewest rays of the sampling,
  default `1`, see [`BeamletOptics.min_num_rays`](@ref).
- `numerical_aperture(s::S)`: the numerical aperture of the half spread angle, required for a
  sampling of a `PointSource`. The fallback throws an `ErrorException`.

# Implementations

- [`BeamletOptics.NoSampling`](@ref): given beams, can not be regenerated
- [`BeamletOptics.DiscRings`](@ref), [`BeamletOptics.DiscSunflower`](@ref),
  [`BeamletOptics.DiscLine`](@ref): `CollimatedSource`
- [`BeamletOptics.ConeRings`](@ref), [`BeamletOptics.ConeSunflower`](@ref),
  [`BeamletOptics.ConeFan`](@ref): `PointSource`
"""
abstract type AbstractSourceSampling end

"""
    NoSampling()

The beams of the source were given to its constructor, e.g.
`CollimatedSource(beams, diameter, pos, dir)`. They can not be regenerated:
[`num_rays!`](@ref) throws an `ArgumentError` and `min_num_rays` returns `nothing`.
"""
struct NoSampling <: AbstractSourceSampling end

"""
    DiscRings(num_rings::Int)

Concentric rings of rays around a center ray on the disc of a [`CollimatedSource`](@ref), the center
counted as the innermost of the `num_rings` rings. The
outermost ring lies on the edge of the disc; the rays are distributed over the rings in
proportion to the area of the disc each ring represents, i.e. to its circumference, but half of
it for the outermost ring, which represents only the inner half of its annulus. Every ray then
represents about the same area (trapezoidal rule), e.g. for spot sizes or point spread
functions, which weight every ray equally. Needs at least `20 num_rings` rays.
"""
struct DiscRings <: AbstractSourceSampling
    num_rings::Int
end

"""
    DiscSunflower()

Sunflower (Fibonacci) sampling of the disc of a [`UniformDiscSource`](@ref): every ray
represents the same area of the disc. The spiral starts with the center ray. Needs at least one
ray.
"""
struct DiscSunflower <: AbstractSourceSampling end

"""
    DiscLine()

Equidistant rays on the line of a [`UniformLineSource`](@ref), i.e. on the diameter of the disc
along the sampling reference vector. The center ray comes first, followed by the other rays from
one end of the line to the other, both ends included. An even number of rays has a second ray at
the center, see [`UniformLineSource`](@ref). Needs at least one ray.
"""
struct DiscLine <: AbstractSourceSampling end

"""
    check_half_angle(θ) -> θ

Returns the half spread angle `θ` [rad] of a cone sampling, or throws an `ArgumentError` unless
`θ < π` in the precision of `θ`.
"""
function check_half_angle(θ::T) where {T <: Real}
    θ < float(T)(π) || throw(ArgumentError("Point source opening half-angle θ must be < π"))
    return θ
end

"""
    ConeRings(num_rings::Int, θ::Real)

Concentric cones of rays around a center ray of a [`PointSource`](@ref), the center counted as the
innermost of the `num_rings` cones, the outermost with the
half spread angle `θ` [rad]. The rays are distributed over the cones in proportion to the solid
angle each cone represents, i.e. to its circumference, but half of it for the outermost cone,
which represents only the inner half of its band. Every ray then represents about the same solid
angle (trapezoidal rule). Needs at least `20 num_rings` rays. `θ ≥ π` throws an `ArgumentError`.
"""
struct ConeRings{T <: Real} <: AbstractSourceSampling
    num_rings::Int
    θ::T
    ConeRings(num_rings::Integer, θ::T) where {T <: Real} = new{T}(num_rings, check_half_angle(θ))
end

"""
    ConeSunflower(θ::Real)

Sunflower (Fibonacci) sampling of the spherical cap with the half spread angle `θ` [rad] of a
[`UniformPointSource`](@ref): every ray represents the same solid angle. The spiral starts with
the center ray. Needs at least one ray. `θ ≥ π` throws an `ArgumentError`.
"""
struct ConeSunflower{T <: Real} <: AbstractSourceSampling
    θ::T
    ConeSunflower(θ::T) where {T <: Real} = new{T}(check_half_angle(θ))
end

"""
    ConeFan(θ::Real)

Rays with equidistant angles in the plane fan of a [`UniformFanSource`](@ref), i.e. in the plane
spanned by the source direction and the sampling reference vector, between `-θ` and `θ` [rad]
to the source direction. The center ray comes first, followed by the other rays from one edge of
the fan to the other, both edges included. An even number of rays has a second ray at the center,
see [`UniformFanSource`](@ref). Needs at least one ray. `θ ≥ π` throws an `ArgumentError`.
"""
struct ConeFan{T <: Real} <: AbstractSourceSampling
    θ::T
    ConeFan(θ::T) where {T <: Real} = new{T}(check_half_angle(θ))
end

"""
    sunflower_fraction(k, n, T)

The area (or solid angle) fraction of the `k`-th of `n` sunflower samples (`k = 0 … n-1`) in the
number type `T`: `0` for the center, `(k + ½)/n` for the others.
"""
sunflower_fraction(k::Int, n::Int, ::Type{T}) where {T <: Real} = k == 0 ? zero(T) : T((k + 0.5) / n)

"""
    fan_coordinates(n, T) -> Vector{T}

The `n` sample positions of a line or fan on `[-1, 1]`: the center first, then equidistant
samples from `-1` to `1`. An even `n` repeats the center.
"""
function fan_coordinates(n::Int, ::Type{T}) where {T <: Real}
    h = (n - 1) ÷ 2
    side = [T(j / h) for j in 1:h]
    x = [zero(T); -reverse(side); side]
    return iseven(n) ? push!(x, zero(T)) : x
end

"""
    check_num_rays(sampling, num_rays)

Throws an `ArgumentError` if `num_rays` is below `min_num_rays(sampling)`.
"""
function check_num_rays(s::AbstractSourceSampling, n::Integer)
    n_min = min_num_rays(s)
    n ≥ n_min || throw(ArgumentError(lazy"$(nameof(typeof(s))) needs at least $n_min rays (passed: $n)"))
    return nothing
end

"""
    source_beams(sampling::AbstractSourceSampling, pos, dir, b1, args..., λ, num_rays, T) -> Vector{<:Beam}

The beams of a source at `pos` [m] along the unit vector `dir`, sampled by `sampling` with
`num_rays` rays of the wavelength `λ` [m] and the number type `T`. `b1` is the unit sampling
reference vector normal to `dir`, see [`BeamletOptics.sampling_basis`](@ref); it sets the
azimuth of the pattern about `dir`. `args` depend on the source: the `diameter` [m] of a
[`CollimatedSource`](@ref), none for a [`PointSource`](@ref). The first beam is the center ray,
which starts at `pos` along `dir`.

Throws an `ArgumentError` if `num_rays` is below [`BeamletOptics.min_num_rays`](@ref) of the
sampling or for [`BeamletOptics.NoSampling`](@ref), and an `ErrorException` if the sampling type
has no method. Used by the constructors of the sources and by [`num_rays!`](@ref).
"""
source_beams(s::AbstractSourceSampling, args...) = error(lazy"source_beams not defined for $(typeof(s))")

source_beams(::NoSampling, args...) =
    throw(ArgumentError("the rays of a source that wraps given beams can not be regenerated"))

#=
CollimatedSource: rings and sunflower on a disc
=#

"""
    source_beams(sampling::Union{DiscRings, DiscSunflower, DiscLine}, pos, dir, b1, diameter, λ, num_rays, T)

Beams of a [`CollimatedSource`](@ref) whose disc with the `diameter` [m] is centered at `pos`
[m] normal to the unit vector `dir`, with the wavelength `λ` [m], sampled by `sampling` with
`num_rays` rays. `b1` is the unit sampling reference vector normal to `dir`, see
`sampling_basis`. For `DiscLine`, the rays lie on the diameter of the disc along `b1`.
"""
function source_beams(s::DiscRings, pos, dir, b1, diameter, λ, num_rays::Int, ::Type{T}) where {T}
    check_num_rays(s, num_rays)
    beams = Vector{Beam{T, Ray{T}}}()
    # Center ray
    push!(beams, Beam(Ray(pos, dir, λ)))
    # Rings
    radii = LinRange(0, diameter / 2, s.num_rings)[2:end]
    counts = rays_per_ring(radii * 2π, num_rays - 1)
    for (r, n) in zip(radii, counts)
        iszero(n) && continue
        R = rotate3d(dir, 2 * T(π) / n)
        offset = b1 * r
        for _ in 1:n
            push!(beams, Beam(pos + offset, dir, λ))
            offset = R * offset
        end
    end
    return beams
end

function source_beams(s::DiscSunflower, pos, dir, b1, diameter, λ, num_rays::Int, ::Type{T}) where {T}
    check_num_rays(s, num_rays)
    R = diameter / 2
    beams = Vector{Beam{T, Ray{T}}}(undef, num_rays)
    b2 = normal3d(dir, b1)
    for k in 0:(num_rays - 1)
        # Equal-area radius
        r = R * √sunflower_fraction(k, num_rays, T)
        φ = k * T(golden_angle)
        x = r * cos(φ) * b1 + r * sin(φ) * b2
        beams[k + 1] = Beam(pos + x, dir, λ)
    end
    return beams
end

function source_beams(s::DiscLine, pos, dir, b1, diameter, λ, num_rays::Int, ::Type{T}) where {T}
    check_num_rays(s, num_rays)
    R = diameter / 2
    return Beam{T, Ray{T}}[Beam(pos + R * c * b1, dir, λ) for c in fan_coordinates(num_rays, T)]
end

#=
PointSource: cones of rings and sunflower on a cap
=#

"""
    source_beams(sampling::Union{ConeRings, ConeSunflower, ConeFan}, pos, dir, b1, λ, num_rays, T)

Beams of a [`PointSource`](@ref) at `pos` [m] around the unit vector `dir` with the wavelength
`λ` [m], sampled by `sampling` (which holds the half spread angle) with `num_rays` rays. `b1` is
the unit sampling reference vector normal to `dir`, see `sampling_basis`. For `ConeFan`, the rays
lie in the plane spanned by `dir` and `b1`. The angles are computed in `T`.
"""
function source_beams(s::ConeRings, pos, dir, b1, λ, num_rays::Int, ::Type{T}) where {T}
    check_num_rays(s, num_rays)
    beams = Vector{Beam{T, Ray{T}}}()
    # Center ray
    push!(beams, Beam(Ray(pos, dir, λ)))
    # Cones
    for d in cone_ring_directions(dir, b1, T(s.θ), s.num_rings, num_rays - 1)
        push!(beams, Beam(pos, d, λ))
    end
    return beams
end

function source_beams(s::ConeSunflower, pos, dir, b1, λ, num_rays::Int, ::Type{T}) where {T}
    check_num_rays(s, num_rays)
    beams = Vector{Beam{T, Ray{T}}}(undef, num_rays)
    b2 = normal3d(dir, b1)
    one_minus_cosθ = 1 - cos(T(s.θ))
    for k in 0:(num_rays - 1)
        # Equal solid angle
        cosϑ = 1 - sunflower_fraction(k, num_rays, T) * one_minus_cosθ
        sinϑ = sqrt(max(0, (1 - cosϑ) * (1 + cosϑ)))
        φ = k * T(golden_angle)
        cdir = cosϑ * dir + sinϑ * (cos(φ) * b1 + sin(φ) * b2)
        beams[k + 1] = Beam(pos, cdir, λ)
    end
    return beams
end

function source_beams(s::ConeFan, pos, dir, b1, λ, num_rays::Int, ::Type{T}) where {T}
    check_num_rays(s, num_rays)
    angles = T(s.θ) * fan_coordinates(num_rays, T)
    return Beam{T, Ray{T}}[Beam(pos, cos(φ) * dir + sin(φ) * b1, λ) for φ in angles]
end

# The fewest rays of a sampling, see `min_num_rays`
min_num_rays(::NoSampling) = nothing
min_num_rays(s::Union{DiscRings, ConeRings}) = 20 * s.num_rings
min_num_rays(::AbstractSourceSampling) = 1

# The numerical aperture of the half spread angle of a PointSource sampling
numerical_aperture(s::AbstractSourceSampling) = error(lazy"numerical_aperture not defined for $(typeof(s))")
numerical_aperture(s::Union{ConeRings, ConeSunflower, ConeFan}) = numerical_aperture(s.θ)
