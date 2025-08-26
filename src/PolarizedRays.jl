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
    intersection::Nullable{Intersection{T}}
    λ::T
    n::T
    E0::Point3{Complex{T}}
    function PolarizedRay{T}(
            pos::AbstractArray{P},
            dir::AbstractArray{D},
            int::Nullable{Intersection},
            λ::L,
            n::N,
            E0::AbstractArray{<:Union{E, Complex{E}}}
        ) where {T, P, D, L, N, E}
        M = promote_type(T, P, D, L, N, E)
        # This test is very important and must be performed for each pol. ray
        if !isorthogonal3d(dir, E0; atol=1e-14)
            throw(ErrorException("Ray dir. and E0 must be orthogonal."))
        end
        return new{M}(
            Point3{M}(pos),
            Point3{M}(dir),
            int,
            M(λ),
            M(n),
            Point3{Complex{M}}(E0)
        )
    end
end

polarization(ray::PolarizedRay) = ray.E0
polarization!(ray::PolarizedRay, new) = (ray.E0 = new)

islinear(ray::PolarizedRay) = islinear(polarization(ray))
iscircular(ray::PolarizedRay) = iscircular(polarization(ray))
iselliptical(ray::PolarizedRay) = iselliptical(polarization(ray))

"""
    PolarizedRay(pos, dir, λ = 1000e-9, E0 = [1, 0, 0])

1 V/m in x-dir.
"""
function PolarizedRay(
        pos::AbstractArray{P},
        dir::AbstractArray{D},
        λ::L = 1000e-9,
        E0::AbstractArray{<:Union{E, Complex{E}}} = [electric_field(1), 0, 0]
    ) where {P <: Real, D <: Real, L <: Real, E<:Real}
    F = promote_type(P, D, L, E)
    return PolarizedRay{F}(
        Point3{F}(pos),
        normalize(Point3{F}(dir)),
        nothing,
        F(λ),
        F(1),
        Point3{Complex{F}}(E0)
    )
end

abstract type AbstractJonesMatrix{T} <: AbstractMatrix{T} end

# Required methods for AbstractArray
Base.size(A::AbstractJonesMatrix) = size(A.data)
Base.getindex(A::AbstractJonesMatrix, i::Int, j::Int) = A.data[i, j]
Base.setindex!(A::AbstractJonesMatrix, v, i::Int, j::Int) = (A.data[i, j] = v)

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
