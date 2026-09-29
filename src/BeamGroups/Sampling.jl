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

How a [`CollimatedSource`](@ref) or [`PointSource`](@ref) sampled its rays: stored in the source,
such that [`set_num_rays!`](@ref) can regenerate them with another number of rays. The methods of
`source_beams` generate the beams of each sampling, the constructors use them as well.
"""
abstract type AbstractSampling end

"""The beams were given to the constructor, they can not be regenerated."""
struct NoSampling <: AbstractSampling end

"""Concentric rings of rays around a center ray on a disc, see [`CollimatedSource`](@ref)."""
struct DiscRings <: AbstractSampling
    num_rings::Int
end

"""Sunflower (Fibonacci) sampling of a disc, see [`UniformDiscSource`](@ref)."""
struct DiscSunflower <: AbstractSampling end

"""Concentric cones of rays with the half spread angle `θ`, see [`PointSource`](@ref)."""
struct ConeRings <: AbstractSampling
    num_rings::Int
    θ::Float64
end

"""Sunflower (Fibonacci) sampling of the cap with the half spread angle `θ`, see [`UniformPointSource`](@ref)."""
struct ConeSunflower <: AbstractSampling
    θ::Float64
end

"""
    source_beams(sampling::AbstractSampling, pos, dir, b1, args..., λ, num_rays, T) -> Vector{<:Beam}

The beams of a source at `pos` along the unit vector `dir`, sampled by `sampling` with `num_rays`
rays of the wavelength `λ`. `b1` is the unit sampling reference vector normal to `dir`, see
[`sampling_basis`](@ref). The arguments between `b1` and `λ` depend on the sampling, e.g. the
`diameter` of a disc. Used by the constructors of the sources and by [`set_num_rays!`](@ref).
"""
function source_beams end

source_beams(::NoSampling, args...) =
    throw(ArgumentError("the rays of a source that wraps given beams can not be regenerated"))

#=
CollimatedSource: rings and sunflower on a disc
=#

"""
    source_beams(sampling, pos, dir, b1, diameter, λ, num_rays, T)

Beams of a [`CollimatedSource`](@ref) at `pos` along the unit vector `dir` with the `diameter`
and the wavelength `λ`, sampled by `sampling` with `num_rays` rays. `b1` is the unit sampling
reference vector normal to `dir`, see `sampling_basis`.
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
    source_beams(sampling, pos, dir, b1, λ, num_rays, T)

Beams of a [`PointSource`](@ref) at `pos` around the unit vector `dir` with the wavelength `λ`,
sampled by `sampling` (which holds the half spread angle) with `num_rays` rays. `b1` is the unit
sampling reference vector normal to `dir`, see `sampling_basis`.
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
