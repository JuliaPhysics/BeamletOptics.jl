"""
    PolarizedRay{T, S} <: AbstractRay{T, S}

A ray type to model the propagation of an electric field vector based on the publication:

**Yun, Garam, Karlton Crabtree, and Russell A. Chipman. "Three-dimensional polarization ray-tracing calculus I: definition and diattenuation." Applied Optics 50.18 (2011): 2855-2865.**

The geometrical ray trajectory is stored in `segment::S` (where `S <: AbstractSegment{T}`). The polarization interaction can be described in local s-p-coordinates
but must be transformed into global coordinates using the method described in the publication above, see also [`_calculate_global_E0`](@ref).

# Fields
- `segment::S`: the trajectory segment of the ray
- `λ::T`: wavelength in [m]
- `n::T`: refractive index along the beam path
- `E0::Point3{Complex{T}}`: complex-valued 3-tuple to represent the electric field in global coordinates

# Jones matrices
In local coordinates the Jones matrices in the case of reflection/refraction are defined as
- reflection: [-rₛ 0; 0 rₚ]
- transmission: [tₛ 0; 0 tₚ]
where r and t are the complex-valued Fresnel coefficients (see also [`fresnel_coefficients`](@ref)).
"""
struct PolarizedRay{T <: Real, S <: AbstractSegment{T}} <: AbstractRay{T, S}
    segment::S
    λ::T
    n::T
    E0::Point3{Complex{T}}

    function PolarizedRay{T, S}(
            seg::S,
            λ::Real,
            n::Real,
            E0::AbstractArray{<:Union{E, Complex{E}}}
        ) where {T <: Real, S <: AbstractSegment{T}, E <: Real}
        if !isorthogonal3d(direction(seg), E0; atol=1e-10)
            error("Ray dir. and E0 must be orthogonal (dot product: $(dot(direction(seg), E0)))")
        end
        return new{T, S}(
            seg,
            T(λ),
            T(n),
            Point3{Complex{T}}(E0)
        )
    end
end

polarization(ray::PolarizedRay) = ray.E0

islinear(ray::PolarizedRay) = islinear(polarization(ray))
iscircular(ray::PolarizedRay) = iscircular(polarization(ray))
iselliptical(ray::PolarizedRay) = iselliptical(polarization(ray))

"""
    PolarizedRay(segment, λ, n, E0)

Constructs a `PolarizedRay` from an existing `segment`.
"""
function PolarizedRay(
        seg::S,
        λ::L,
        n::N,
        E0::AbstractArray{<:Union{E, Complex{E}}}
    ) where {T <: Real, S <: AbstractSegment{T}, L <: Real, N <: Real, E <: Real}
    F = promote_type(T, L, N, E)
    if F != T
        # Re-promote segment if needed
        seg_promoted = S <: OpenSegment ? OpenSegment(position(seg), direction(seg), F(accumulated_opl(seg))) :
                       LineSegment(position(seg), direction(seg), Intersection(F(length(seg.intersection)), Point3{F}(position(seg.intersection)), Point3{F}(normal3d(seg.intersection))), F(accumulated_opl(seg)))
        return PolarizedRay{F, typeof(seg_promoted)}(seg_promoted, F(λ), F(n), Point3{Complex{F}}(E0))
    else
        return PolarizedRay{T, S}(seg, T(λ), T(n), Point3{Complex{T}}(E0))
    end
end

"""
    PolarizedRay(ray::PolarizedRay, segment)

Resolves or replaces the trajectory segment of an existing `ray` with a new `segment`.
"""
PolarizedRay{T, S}(ray::PolarizedRay{T}, seg::S) where {T <: Real, S <: AbstractSegment{T}} = PolarizedRay{T, S}(seg, ray.λ, ray.n, ray.E0)
function PolarizedRay(ray::PolarizedRay{T}, seg::S) where {T <: Real, S <: AbstractSegment{T}}
    return PolarizedRay(seg, ray.λ, ray.n, ray.E0)
end

function (::Type{PolarizedRay{T}})(
        pos::AbstractArray{P},
        dir::AbstractArray{D},
        λ::L = 1000e-9,
        n::N = 1.0,
        E0::AbstractArray = [electric_field(1), 0, 0]
    ) where {T <: Real, P <: Real, D <: Real, L <: Real, N <: Real}
    seg = OpenSegment(pos, dir, zero(T))
    return PolarizedRay{T, typeof(seg)}(
        seg,
        T(λ),
        T(n),
        Point3{Complex{T}}(E0)
    )
end

function PolarizedRay(
        pos::AbstractArray{P},
        dir::AbstractArray{D},
        λ::L = 1000e-9,
        E0::AbstractArray{<:Union{E, Complex{E}}} = [electric_field(1), 0, 0]
    ) where {P <: Real, D <: Real, L <: Real, E <: Real}
    F = promote_type(P, D, L, E)
    return PolarizedRay{F}(pos, dir, F(λ), F(1), E0)
end

function PolarizedRay(
        pos::AbstractArray{P},
        dir::AbstractArray{D},
        λ::L,
        n::N,
        E0::AbstractArray{<:Union{E, Complex{E}}}
    ) where {P <: Real, D <: Real, L <: Real, N <: Real, E <: Real}
    F = promote_type(P, D, L, N, E)
    return PolarizedRay{F}(pos, dir, F(λ), F(n), E0)
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

function _calculate_global_E0(::AbstractObject, ray::PolarizedRay, out_dir::AbstractArray, J::LocalJonesBasis, int::Nullable{<:AbstractIntersection} = nothing)
    # Update Jones matrix according to global object orientation
    in_dir = direction(ray)
    normal = isnothing(int) ? (isnothing(intersection(ray)) ? Point3(0.0, 0.0, 1.0) : normal3d(intersection(ray))) : normal3d(int)
    E0 = polarization(ray)
    P = _calculate_global_E0(in_dir, out_dir, normal, J)
    return P*E0
end
