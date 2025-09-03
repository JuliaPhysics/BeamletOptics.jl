"""
    AbstractJonesPolarizer <: AbstractObject

Represents infinitesimally thin components that change the polarization state of incoming [`PolarizedRay`](@ref)s via global Jones matrix calculus.
Rather than using the generic Yun ray tracing scheme as referred to in the `PolarizedRay` docs, this element interacts with
the global E-field vector `E0` by using a [`GlobalJonesBasis`](@ref) and projecting the entries into the transverse plane defined
by the incoming ray direction and orthogonal E-field vector. This approach is partially inspired by the publication:

**Jan Korger et al., "The polarization properties of a tilted polarizer," Opt. Express 21, 27032-27042 (2013)**

# Implementation reqs.

Subtypes of `AbstractJonesPolarizer` should implement all supertype requirements.

## Interaction logic

The [`GlobalJonesBasis`](@ref) tracks the rotation in 3D-space via the [`orientation`](@ref) of the attached [`AbstractShape`](@ref).
The polarization matrix `P` is calculated by projecting the previous matrix into the incoming orthogonal plane of polarization.
Refer to the [`_calculate_global_E0`](@ref) implementation for more information.

!!! info
    The validity of this approach is still under consideration for non-normal incidence.
"""
abstract type AbstractJonesPolarizer{T, S} <: AbstractObject{T, S} end

"""
    interact3d(AbstractSystem, AbstractJonesPolarizer, Beam, Ray)

Non‑polarized [`Ray`](@ref)s pass through an [`AbstractJonesPolarizer`](@ref). without modification.
"""
function interact3d(::AbstractSystem, ::AbstractJonesPolarizer, ::Beam{T,R},
        ray::R) where {T<:Real, R<:Ray{T}}
    pos = position(ray) + length(ray) * direction(ray)
    return BeamInteraction{T,R}(nothing,
        Ray{T}(pos, direction(ray), nothing, wavelength(ray), refractive_index(ray)))
end

abstract type AbstractJonesMatrix{T} <: AbstractMatrix{T} end

# Required methods for AbstractArray
Base.size(A::AbstractJonesMatrix) = size(A.data)
Base.getindex(A::AbstractJonesMatrix, i::Int, j::Int) = getindex(A.data, i, j)
Base.setindex!(A::AbstractJonesMatrix, v, i::Int, j::Int) = setindex!(A.data, v, i, j)

"Getter fct. for the static array in the `AbstractJonesMatrix`"
static_data(A::AbstractJonesMatrix) = A.data

"""
    LocalJonesBasis

Stores the s-p-basis Jones matrix coefficients. Must be defined for x-y-aligned elements
where z is the optical axis.
"""
struct LocalJonesBasis{T} <: AbstractJonesMatrix{T}
    data::SMatrix{3, 3, T, 9}
end

function SPBasis(j11::Number, j12::Number, j21::Number, j22::Number)
    T = promote_type(
        typeof(j11),
        typeof(j12),
        typeof(j21),
        typeof(j22),
    )
    return LocalJonesBasis{T}(SMatrix{3, 3, T, 9}(j11, j12, 0, j21, j22, 0, 0, 0, 1))
end

"""
    GlobalJonesBasis <: AbstractJonesMatrix

Stores the Jones matrix entries for a polarizing optical element that is aligned
with the global y-axis as the optical axis.
"""
struct GlobalJonesBasis{T} <: AbstractJonesMatrix{T}
    data::SMatrix{3, 3, T, 9}
end

# Generic path for any AbstractMatrix
function GlobalJonesBasis(J::AbstractMatrix)
    size(J) != (3, 3) && throw(ArgumentError("GlobalJonesBasis expects a 3×3 matrix"))
    T = eltype(J)
    return GlobalJonesBasis{T}(SMatrix{3,3,T,9}(J))
end

# Data type conversion constructor
GlobalJonesBasis{T}(J::GlobalJonesBasis) where {T} =
    GlobalJonesBasis{T}(SMatrix{3,3,T,9}(static_data(J)))

XYBasis(j11::Number, j12::Number, j21::Number, j22::Number) = GlobalJonesBasis(@SArray([j11 j12 0; j21 j22 0; 0 0 1]))
XZBasis(j11::Number, j12::Number, j21::Number, j22::Number) = GlobalJonesBasis(@SArray([j11 0 j12; 0 1 0; j21 0 j22]))
YZBasis(j11::Number, j12::Number, j21::Number, j22::Number) = GlobalJonesBasis(@SArray([1 0 0; 0 j22 j21; 0 j12 j11]))

"""
    _calculate_global_E0(in_dir, out_dir, normal, J)

Calculates the resulting polarization matrix as per the publication by Yun et al. for each surface interaction.
If the `normal`- vector and the `in`- and `out`-directions of propagation are in parallel, an arbitrary basis
is chosen for the s- and p-components.

# Arguments
- `in_dir`: propagation direction before surface interaction
- `out_dir`: propagation direction after surface interaction
- `normal`: surface normal at the point of intersection
- `J`: Jones matrix extended to 3x3, e.g. [-rₛ 0 0; 0 rₚ 0; 0 0 1] for reflection
"""
function _calculate_global_E0(in_dir::AbstractArray, out_dir::AbstractArray, normal::AbstractArray, J::LocalJonesBasis)
    # Choose basis vectors
    if !isparallel3d(in_dir, out_dir)
        v = out_dir
    else
        v = normal
    end
    # test if in-dir and normal are parallel
    if isparallel3d(in_dir, normal)
        # Does this really always work for normal s-p-incidence?
        v = normal3d(in_dir)
    end
    # Calculate support vector
    s = cross(in_dir, v)
    s = normalize(s)
    # Calculate transforms
    p1 = cross(in_dir, s)
    O_in = vcat(s', p1', in_dir')
    # Fallback method as per eq. 17
    if isparallel3d(in_dir, out_dir) && !(in_dir ≈ -out_dir)
        O_out = hcat(s, p1, in_dir)
    else
        p2 = cross(out_dir, s)
        O_out = hcat(s, p2, out_dir)
    end
    # Calculate new E0
    P = O_out * J * O_in
    return P
end

# unwrap the Local/GlobalJonesBasis types to allow static array optimization to happen
function _calculate_global_E0(in_dir::AbstractArray, out_dir::AbstractArray, normal::AbstractArray, J::AbstractJonesMatrix)
    _calculate_global_E0(in_dir, out_dir, normal, static_data(J))
end

function _calculate_global_E0(::AbstractObject, ray::PolarizedRay, out_dir::AbstractArray, J::LocalJonesBasis)
    # Update Jones matrix according to global object orientation
    in_dir = direction(ray)
    normal = normal3d(intersection(ray))
    E0 = polarization(ray)
    P = _calculate_global_E0(in_dir, out_dir, normal, J)
    return P*E0
end

function _calculate_global_E0(object::AbstractObject, ray::PolarizedRay, out_dir::AbstractArray, J::GlobalJonesBasis)
    in_dir = direction(ray)
    E0 = polarization(ray)
    # Transform Jones matrix according to global object orientation
    R = orientation(object)
    P = R * J * transpose(R)

    Q_in = I - in_dir * transpose(in_dir)
    Q_out = I - out_dir * transpose(out_dir)
    P = Q_out * P * Q_in
    return P * E0
end

include("PolarizationFilter.jl")
include("Waveplates.jl")
