"""
    PolarizedRay{T} <: AbstractRay{T}

A ray type to model the propagation of an electric field vector based on the publication:

**Yun, Garam, Karlton Crabtree, and Russell A. Chipman. "Three-dimensional polarization ray-tracing calculus I: definition and diattenuation." Applied optics 50.18 (2011): 2855-2865.**

The geometrical ray description is identical to the standard [`Ray`](@ref). The polarization interaction can be described in local s-p-coordinates
but must be transformed into global coordinates using the method described in the publication above, see also [`_calculate_global_E0`](@ref).

# Fields

- `pos`: a point in R³ that describes the `Ray` origin
- `dir`: a normalized vector in R³ that describes the `Ray` direction
- `intersection`: refer to [`Intersection`](@ref)
- `λ`: wavelength in [m]
- `n`: refractive index along the beam path
- `E0`: complex-valued 3-tuple to represent the electric field in global coordinates

# Jones matrices

In local coordinates the Jones matrices in the case of reflection/refraction are defined as

- reflection: [-rₛ 0; 0 rₚ]
- transmission: [tₛ 0; 0 tₚ]

where r and t are the complex-valued Fresnel coefficients (see also [`fresnel_coefficients`](@ref)).

# Additional information

!!! warning "Field vector"
    It is assumed that the electric field vector ``E_0`` stays orthogonal to the direction of propagation throughout the optical system.

!!! warning "Intensity"
    E0 can not be converted into an [`intensity`](@ref) value, since a single `PolarizedRay` can not directly model the change in intensity during imaging by an optical system.
"""
mutable struct PolarizedRay{T} <: AbstractRay{T}
    pos::Point3{T}
    dir::Point3{T}
    intersection::Nullable{Intersection}
    λ::T
    n::T
    E0::Point3{Complex{T}}
end

polarization(ray::PolarizedRay) = ray.E0
polarization!(ray::PolarizedRay, new) = (ray.E0 = new)

"""
    PolarizedRay(pos, dir, λ = 1000e-9, E0 = [1, 0, 0])

1 V/m in x-dir.
"""
function PolarizedRay(pos::AbstractArray{P},
        dir::AbstractArray{D},
        λ = 1000e-9,
        E0 = [electric_field(1), 0, 0]) where {P <: Real, D <: Real}
    F = promote_type(P, D)
    # Test if E0 is orthogonal to dir. of propagation
    if !isapprox(dot(dir, E0), 0, atol=1e-14)
        error(lazy"Ray dir. and E0 must be orthogonal (n=$(dot(dir, E0)))")
    end
    return PolarizedRay{F}(
        Point3{F}(pos),
        normalize(Point3{F}(dir)),
        nothing,
        λ,
        F(1),
        E0)
end

#FIXME: revert to static array for speed
abstract type AbstractJonesMatrix{T} <: AbstractMatrix{T} end

# Required methods for AbstractArray
Base.size(A::AbstractJonesMatrix) = size(A.data)
Base.getindex(A::AbstractJonesMatrix, i::Int, j::Int) = A.data[i, j]
Base.setindex!(A::AbstractJonesMatrix, v, i::Int, j::Int) = (A.data[i, j] = v)

struct SPBasis{T} <: AbstractJonesMatrix{T}
    data::Matrix{T}
end

function SPBasis(j11::Number, j12::Number, j21::Number, j22::Number)
    T = promote_type(
        typeof(j11),
        typeof(j12),
        typeof(j21),
        typeof(j22),
    )
    return SPBasis{T}(T.([j11 j12 0; j21 j22 0; 0 0 1]))
end

struct XYBasis{T} <: AbstractJonesMatrix{T}
    data::Matrix{T}
end

function XYBasis(j11::Number, j12::Number, j21::Number, j22::Number)
    T = promote_type(
        typeof(j11),
        typeof(j12),
        typeof(j21),
        typeof(j22),
    )
    return XYBasis{T}(T.([j11 j12 0; j21 j22 0; 0 0 1]))
end

"""
    _calculate_global_E0(in_dir, out_dir, J, E0)

Calculates the resulting polarization vector as per the publication by Yun et al. for each surface interaction.
If the `in`- and `out`-directions of propagation are parallel, an arbitrary basis is chosen for the s- and p-components.

# Arguments
- `in_dir`: propagation direction before surface interaction
- `out_dir`: propagation direction after surface interaction
- `J`: Jones matrix extended to 3x3, e.g. [-rₛ 0 0; 0 rₚ 0; 0 0 1] for reflection
- `E0`: Polarization vector before surface interaction
"""
function _calculate_global_E0(in_dir::AbstractArray, out_dir::AbstractArray, normal::AbstractArray, J::AbstractMatrix, E0::AbstractArray)
    # Choose basis vectors
    if !isparallel3d(in_dir, out_dir)
        v = out_dir
    else
        v = normal
    end
    # test if in-dir and normal are parallel
    if isparallel3d(in_dir, normal)
        # v = normal3d(in_dir)
        @info round.(J, digits=3)
        return J * E0
    end
    # Calculate support vector
    s = cross(in_dir, v)
    s = normalize(s)
    # Calculate transforms
    p1 = cross(in_dir, s)
    p2 = cross(out_dir, s)
    O_in = [s'; p1'; in_dir']
    O_out = [s p2 out_dir]
    # Calculate new E0
    return O_out * J * O_in * E0
end

function _calculate_global_E0(object::AbstractObject, ray::PolarizedRay, out_dir::AbstractArray, J::XYBasis)
    in_dir = direction(ray)
    normal = normal3d(intersection(ray))
    E0 = polarization(ray)
    # Update Jones matrix according to global object orientation
    R = orientation(object)
    J_glbl = R * J * transpose(R)
    E0_local = _calculate_global_E0(in_dir, out_dir, normal, J_glbl, E0)
    # Retransform
    E0_global = transpose(R) * E0_local
    return E0_global
end

function _calculate_global_E0(::AbstractObject, ray::PolarizedRay, out_dir::AbstractArray, J::SPBasis)
    # Update Jones matrix according to global object orientation
    in_dir = direction(ray)
    normal = normal3d(intersection(ray))
    E0 = polarization(ray)
    return _calculate_global_E0(in_dir, out_dir, normal, J, E0)
end