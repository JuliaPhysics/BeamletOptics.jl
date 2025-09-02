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