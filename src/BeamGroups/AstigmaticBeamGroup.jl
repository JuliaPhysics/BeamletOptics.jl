"""
    AstigmaticBeamGroup{T, R} <: AbstractBeamGroup{T, R}

A generic container for groups of [`AstigmaticGaussianBeamlet`](@ref)s.

# Fields

- `beams`: vector of all beamlets
- `center`: source position, pivot for rotations
- `orientation`: right-handed orthonormal matrix, columns are the sampling reference vector, the central source direction and their cross product, see [`AbstractBeamGroup`](@ref)
"""
mutable struct AstigmaticBeamGroup{T, R <: AbstractRay{T}} <: AbstractBeamGroup{T, R}
    beams::Vector{AstigmaticGaussianBeamlet{T}}
    center::Point3{T}
    orientation::SMatrix{3, 3, T, 9}
end

"""
    AstigmaticBeamGroup(beams, pos, dir::AbstractVector)
    AstigmaticBeamGroup(beams, pos, orientation::AbstractMatrix)

Wraps existing `beams` into an [`AstigmaticBeamGroup`](@ref) with the source position `pos`.
The group orientation is either derived from the central direction `dir` (the sampling
reference vector is then picked deterministically via [`normal3d`](@ref)), or passed
explicitly as a right-handed orthonormal 3x3 `orientation` matrix whose second column is the
central direction. An invalid `orientation` throws an `ArgumentError`.
"""
function AstigmaticBeamGroup(beams::Vector{AstigmaticGaussianBeamlet{T}}, pos, dir::AbstractVector) where {T}
    d = normalize(dir)
    M = _group_orientation(d, _sampling_basis(d, nothing, T), T)
    _check_kinematic_members(beams)
    return AstigmaticBeamGroup{T, PolarizedRay{T}}(beams, Point3{T}(pos), M)
end

function AstigmaticBeamGroup(beams::Vector{AstigmaticGaussianBeamlet{T}}, pos, orientation::AbstractMatrix) where {T}
    _check_kinematic_members(beams)
    return AstigmaticBeamGroup{T, PolarizedRay{T}}(beams, Point3{T}(pos), _check_orientation(orientation, T))
end
