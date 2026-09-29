"""
    _sampling_basis(dir, basis, T)

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
_sampling_basis(dir::AbstractVector, ::Nothing, ::Type{<:Real}) = normal3d(dir)

function _sampling_basis(dir::AbstractVector, basis::AbstractVector, ::Type{T}) where {T <: Real}
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
    _group_orientation(dir, e1, T) -> SMatrix{3,3,T,9}

Returns the orientation matrix of a beam group with the columns `(e1, dir, e1 × dir)`, i.e.
local x = sampling reference vector `e1`, local y = optical axis `dir`. Both inputs are
normalized, `e1` must be normal to `dir`.
"""
function _group_orientation(dir::AbstractVector, e1::AbstractVector, ::Type{T}) where {T <: Real}
    d = normalize(SVector{3, T}(dir))
    x = normalize(SVector{3, T}(e1))
    return SMatrix{3, 3, T, 9}(hcat(x, d, cross(x, d)))
end

"""
    _check_orientation(M, T) -> SMatrix{3,3,T,9}

Validates that `M` is a right-handed orthonormal 3x3 matrix, i.e. `MᵀM ≈ I` and `det M ≈ 1`
within `√eps(T)`, and converts it into an `SMatrix{3,3,T,9}`. Throws an `ArgumentError` otherwise.
"""
function _check_orientation(M::AbstractMatrix, ::Type{T}) where {T <: Real}
    if size(M) != (3, 3)
        throw(ArgumentError("Beam group orientation must be a 3x3 matrix"))
    end
    S = SMatrix{3, 3, T, 9}(M)
    tol = sqrt(eps(T))
    # determinant as triple product of the columns
    detS = dot(S[:, 1], cross(S[:, 2], S[:, 3]))
    if !isapprox(transpose(S) * S, SMatrix{3, 3, T, 9}(I); atol = tol) ||
       !isapprox(detS, one(T); atol = tol)
        throw(ArgumentError("Beam group orientation must be a right-handed orthonormal matrix"))
    end
    return S
end

"""
    _AbstractSampling

How a [`CollimatedSource`](@ref) or [`PointSource`](@ref) sampled its rays: stored in the source,
such that [`set_num_rays!`](@ref) can regenerate them with another number of rays. The methods of
`_source_beams` generate the beams of each sampling, the constructors use them as well.
"""
abstract type _AbstractSampling end

"""The beams were given to the constructor, they can not be regenerated."""
struct _NoSampling <: _AbstractSampling end

"""Concentric rings of rays around a center ray on a disc, see [`CollimatedSource`](@ref)."""
struct _DiscRings <: _AbstractSampling
    num_rings::Int
end

"""Sunflower (Fibonacci) sampling of a disc, see [`UniformDiscSource`](@ref)."""
struct _DiscSunflower <: _AbstractSampling end

"""Concentric cones of rays with the half spread angle `θ`, see [`PointSource`](@ref)."""
struct _ConeRings <: _AbstractSampling
    num_rings::Int
    θ::Float64
end

"""Sunflower (Fibonacci) sampling of the cap with the half spread angle `θ`, see [`UniformPointSource`](@ref)."""
struct _ConeSunflower <: _AbstractSampling
    θ::Float64
end

_source_beams(::_NoSampling, args...) =
    throw(ArgumentError("the rays of a source that wraps given beams can not be regenerated"))

"""
    set_num_rays!(group, n)

Regenerates the rays of the source `group` with `n` rays: a [`CollimatedSource`](@ref) or
[`PointSource`](@ref), including [`UniformDiscSource`](@ref) and [`UniformPointSource`](@ref). The
rays are sampled as by the constructor of the source (concentric rings with the same `num_rings`,
or the sunflower pattern), in the current pose of the source and at its wavelength, i.e. the
source equals a new one with `num_rays = n` at its position and orientation. The previous beams,
including their traced rays, are replaced: solve the system again afterwards.

Ring sources need `n ≥ 20 num_rings` and throw an `ErrorException` otherwise, like their
constructors. A source built from given beams, e.g. `CollimatedSource(beams, diameter, pos, dir)`,
can not be regenerated and throws an `ArgumentError`.
"""
function set_num_rays! end

"""
    source_wavelength(bg::AbstractBeamGroup)

The wavelength of the source `bg` [m], i.e. of its first ray.
"""
source_wavelength(bg::AbstractBeamGroup) = wavelength(first(rays(first(beams(bg)))))

include("PointSource.jl")
include("CollimatedSource.jl")

"""
    min_num_rays(src) -> Union{Nothing, Int}

The fewest rays with which [`set_num_rays!`](@ref) regenerates the source `src`: `20 num_rings`
for ring sources, `1` for sunflower sources, and `nothing` if the rays of `src` can not be
regenerated, e.g. a source built from given beams or any other beam or beam group.
"""
min_num_rays(src) = nothing
min_num_rays(src::Union{CollimatedSource, PointSource}) = min_num_rays(src.sampling)
min_num_rays(::_NoSampling) = nothing
min_num_rays(s::Union{_DiscRings, _ConeRings}) = 20 * s.num_rings
min_num_rays(::_AbstractSampling) = 1
include("AstigmaticBeamGroup.jl")
include("AstigmaticBeamletSources.jl")
include("BeamletDecompositions.jl")
