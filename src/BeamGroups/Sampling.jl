#=
Sampling of the rays of the sources: the reference vector of the azimuthal sampling, the sampling
types stored in a source and the beams they generate, see `source_beams`
=#

"""
    sampling_basis(dir, basis, T)

Returns the unit-length reference vector that seeds the azimuthal sampling of a beam group.
`dir` must already be normalized.

If `basis` is `nothing`, a vector normal to `dir` is picked deterministically via
[`normal3d`](@ref). Otherwise `basis` is projected into the plane normal to `dir`, which
lets a caller rotate the sampling pattern of a source about its own axis.

Throws if `basis` has no significant component in that plane, i.e. if it is zero or
parallel to `dir`. The test is relative to `norm(basis)` because projecting an
unnormalized vector leaves a residual that scales with its length; an absolute
tolerance would let a parallel `basis` through and yield `NaN` sampling.
"""
sampling_basis(dir::AbstractVector, ::Nothing, ::Type{<:Real}) = normal3d(dir)

function sampling_basis(dir::AbstractVector, basis::AbstractVector, ::Type{T}) where {T <: Real}
    b1 = basis - dot(basis, dir) * dir
    if norm(b1) ≤ sqrt(eps(T)) * norm(basis)
        throw(ErrorException("Source `basis` must not be zero or parallel to `dir`"))
    end
    return normalize(b1)
end

"""
    _GOLDEN_ANGLE

The golden angle `π(3 - √5) ≈ 2.39996` rad, i.e. the azimuthal increment of the sunflower
(Fibonacci) sampling shared by [`UniformDiscSource`](@ref) and [`UniformPointSource`](@ref).
"""
const _GOLDEN_ANGLE = π * (3 - √5)

"""
    AbstractSampling

How the rays of a [`CollimatedSource`](@ref) or [`PointSource`](@ref) are sampled. The source
stores its sampling in the field `sampling`, such that [`set_num_rays!`](@ref) can regenerate its
rays with another number of rays. The constructors of the sources generate their beams with the
same methods.

# Implementation reqs.

A subtype `S <: AbstractSampling` implements:

- `source_beams(s::S, pos, dir, b1, args..., λ, num_rays::Int, T) -> Vector{<:Beam}`: the beams
  of the source, see [`BeamletOptics.source_beams`](@ref). `args` are the arguments of the source
  type, i.e. the `diameter` [m] for a `CollimatedSource` and none for a `PointSource`. The
  sampling pattern must depend on the pose only via `pos`, `dir` and `b1`, i.e. it must move
  rigidly with the source. Throws an `ErrorException` if `num_rays` is below
  `min_num_rays(s)`.
- `min_num_rays(s::S) -> Union{Nothing, Int}` (optional): the fewest rays of the sampling,
  default `1`, see [`BeamletOptics.min_num_rays`](@ref).

Parameters of the sampling other than the number of rays (e.g. the number of rings or the
half spread angle) are fields of the subtype.

# Implementations

- [`BeamletOptics.NoSampling`](@ref): given beams, can not be regenerated
- [`BeamletOptics.DiscRings`](@ref), [`BeamletOptics.DiscSunflower`](@ref): `CollimatedSource`
- [`BeamletOptics.ConeRings`](@ref), [`BeamletOptics.ConeSunflower`](@ref): `PointSource`
"""
abstract type AbstractSampling end

"""
    NoSampling()

The beams of the source were given to its constructor, e.g.
`CollimatedSource(beams, diameter, pos, dir)`. They can not be regenerated:
[`set_num_rays!`](@ref) throws an `ArgumentError` and `min_num_rays` returns `nothing`.
"""
struct NoSampling <: AbstractSampling end

"""
    DiscRings(num_rings::Int)

Concentric rings of rays around a center ray on the disc of a [`CollimatedSource`](@ref). The
outermost ring lies on the edge of the disc; the rays are distributed over the rings in
proportion to their circumference. Needs at least `20 num_rings` rays.
"""
struct DiscRings <: AbstractSampling
    num_rings::Int
end

"""
    DiscSunflower()

Sunflower (Fibonacci) sampling of the disc of a [`UniformDiscSource`](@ref): every ray
represents the same area of the disc. Needs at least one ray.
"""
struct DiscSunflower <: AbstractSampling end

"""
    ConeRings(num_rings::Int, θ::Float64)

Concentric cones of rays around a center ray of a [`PointSource`](@ref), the outermost with the
half spread angle `θ` [rad]. The rays are distributed over the cones in proportion to their
circumference. Needs at least `20 num_rings` rays.
"""
struct ConeRings <: AbstractSampling
    num_rings::Int
    θ::Float64
end

"""
    ConeSunflower(θ::Float64)

Sunflower (Fibonacci) sampling of the spherical cap with the half spread angle `θ` [rad] of a
[`UniformPointSource`](@ref): every ray represents the same solid angle. Needs at least one ray.
"""
struct ConeSunflower <: AbstractSampling
    θ::Float64
end

"""
    source_beams(sampling::AbstractSampling, pos, dir, b1, args..., λ, num_rays, T) -> Vector{<:Beam}

The beams of a source at `pos` [m] along the unit vector `dir`, sampled by `sampling` with
`num_rays` rays of the wavelength `λ` [m] and the number type `T`. `b1` is the unit sampling
reference vector normal to `dir`, see [`BeamletOptics.sampling_basis`](@ref); it sets the
azimuth of the pattern about `dir`. `args` depend on the source: the `diameter` [m] of a
[`CollimatedSource`](@ref), none for a [`PointSource`](@ref).

Throws an `ErrorException` if `num_rays` is below [`BeamletOptics.min_num_rays`](@ref) of the
sampling, and an `ArgumentError` for [`BeamletOptics.NoSampling`](@ref). Used by the
constructors of the sources and by [`set_num_rays!`](@ref).
"""
function source_beams end

source_beams(::NoSampling, args...) =
    throw(ArgumentError("the rays of a source that wraps given beams can not be regenerated"))

#=
CollimatedSource: rings and sunflower on a disc
=#

"""
    source_beams(sampling::Union{DiscRings, DiscSunflower}, pos, dir, b1, diameter, λ, num_rays, T)

Beams of a [`CollimatedSource`](@ref) whose disc with the `diameter` [m] is centered at `pos`
[m] normal to the unit vector `dir`, with the wavelength `λ` [m], sampled by `sampling` with
`num_rays` rays. `b1` is the unit sampling reference vector normal to `dir`, see
`sampling_basis`.
"""
function source_beams(s::DiscRings, pos, dir, b1, diameter, λ, num_rays::Int, ::Type{T}) where {T}
    num_rings = s.num_rings
    if num_rays < num_rings * 20
        throw(ErrorException("No. of rays should be atleast 20x no. of rings (passed: $num_rays, req: $(num_rings*20))"))
    end
    # define buffer
    beams = Vector{Beam{T, Ray{T}}}()
    push!(beams, Beam(Ray(pos, dir, λ)))
    num_rays -= 1
    # setup concentric beam ring radii
    r_max = diameter / 2
    radii = LinRange(0, r_max, num_rings)[2:end]
    # calculate total accumulated circumference of all rings
    circm = radii * 2π
    total = sum(circm)
    ds = total / num_rays
    # calculate number of rays per ring
    n_rays = round.(Int, circm / ds)
    # correct n_rays to match num_rays
    n_rays[end] += (num_rays - sum(n_rays))
    # Generate beam rings
    for (i, r) in enumerate(radii)
        numEl = n_rays[i]
        if iszero(numEl)
            continue
        end
        dphi = 2π / numEl
        RotMat = rotate3d(dir, dphi)
        helper = b1 * r
        for _ in 1:numEl
            push!(beams, Beam(pos + helper, dir, λ))
            helper = RotMat * helper
        end
    end
    return beams
end

function source_beams(::DiscSunflower, pos, dir, e1, diameter, λ, num_rays::Int, ::Type{T}) where {T}
    if num_rays < 1
        throw(ErrorException("No. of rays must be at least 1 (passed: $num_rays)"))
    end
    R = diameter / 2
    beams = Vector{Beam{T, Ray{T}}}(undef, num_rays)
    # orthogonal basis in the pupil plane
    e2 = normal3d(dir, e1)
    for k in 0:(num_rays - 1)
        ρ = √((k + 0.5) / num_rays)     # equal-area radius
        φ = k * _GOLDEN_ANGLE
        r = R * ρ
        x = r * cos(φ) * e1 + r * sin(φ) * e2
        beams[k + 1] = Beam(pos + x, dir, λ)
    end
    return beams
end

#=
PointSource: cones of rings and sunflower on a cap
=#

"""
    source_beams(sampling::Union{ConeRings, ConeSunflower}, pos, dir, b1, λ, num_rays, T)

Beams of a [`PointSource`](@ref) at `pos` [m] around the unit vector `dir` with the wavelength
`λ` [m], sampled by `sampling` (which holds the half spread angle) with `num_rays` rays. `b1` is
the unit sampling reference vector normal to `dir`, see `sampling_basis`.
"""
function source_beams(s::ConeRings, pos, dir, b1, λ, num_rays::Int, ::Type{T}) where {T}
    num_rings, θ = s.num_rings, s.θ
    if num_rays < num_rings * 20
        throw(ErrorException("No. of rays should be atleast 20x no. of rings (passed: $num_rays, req: $(num_rings*20))"))
    end
    # define basis vectors
    b2 = normal3d(dir, b1)
    θ_NA = LinRange(0, θ, num_rings)
    # define buffer
    beams = Vector{Beam{T, Ray{T}}}()
    push!(beams, Beam(Ray(pos, dir, λ)))
    num_rays -= 1
    # calculate total accumulated circumference of all rings
    ndirs = [rotate3d(b2, step(θ_NA) * i) * dir for i in eachindex(θ_NA[2:end])]
    circm = norm.(ndirs .- dot.(ndirs, Ref(dir)) .* Ref(dir)) .* Ref(2π)
    total = sum(circm)
    ds = total / num_rays
    # calculate number of rays per ring
    n_rays = round.(Int, circm / ds)
    # correct n_rays to match num_rays
    n_rays[end] += (num_rays - sum(n_rays))
    for (i, ndir) in enumerate(ndirs)
        numEl = n_rays[i]
        if iszero(numEl)
            continue
        end
        dphi = 2π / numEl
        RotMat = rotate3d(dir, dphi)
        cdir = ndir
        for _ in 1:numEl
            push!(beams, Beam(pos, cdir, λ))
            # rotate vector (not-thread safe!)
            cdir = RotMat * cdir
        end
    end
    return beams
end

function source_beams(s::ConeSunflower, pos, dir, e1, λ, num_rays::Int, ::Type{T}) where {T}
    if num_rays < 1
        throw(ErrorException("No. of rays must be at least 1 (passed: $num_rays)"))
    end
    beams = Vector{Beam{T, Ray{T}}}(undef, num_rays)
    # orthogonal basis normal to the cone axis
    e2 = normal3d(dir, e1)
    one_minus_cosθ = 1 - cos(s.θ)
    for k in 0:(num_rays - 1)
        cosϑ = 1 - (k + 0.5) / num_rays * one_minus_cosθ    # equal solid angle
        sinϑ = sqrt(max(0, (1 - cosϑ) * (1 + cosϑ)))
        φ = k * _GOLDEN_ANGLE
        cdir = cosϑ * dir + sinϑ * (cos(φ) * e1 + sin(φ) * e2)
        beams[k + 1] = Beam(pos, cdir, λ)
    end
    return beams
end

# The fewest rays of a sampling, see `min_num_rays`
min_num_rays(::NoSampling) = nothing
min_num_rays(s::Union{DiscRings, ConeRings}) = 20 * s.num_rings
min_num_rays(::AbstractSampling) = 1
