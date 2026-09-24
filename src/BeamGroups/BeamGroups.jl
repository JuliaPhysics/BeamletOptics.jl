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

include("PointSource.jl")
include("CollimatedSource.jl")
include("AstigmaticBeamGroup.jl")
include("AstigmaticBeamletSources.jl")
include("BeamletDecompositions.jl")
