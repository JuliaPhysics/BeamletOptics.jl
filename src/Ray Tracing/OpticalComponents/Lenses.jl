"""
    AbstractRefractiveOptic <: AbstractObject

A generic type to represent `AbstractObject`s which refract incoming rays. The main function of `interact3d` should be akin to [`refraction3d`](@ref).

# Implementation reqs.

Subtypes of `AbstractRefractiveOptic` should implement all supertype reqs. as well as:

# Fields

- `n::Function`: a function which returns the [`RefractiveIndex`](@ref) for a wavelength λ
"""
abstract type AbstractRefractiveOptic{T, S <: AbstractShape{T}, F} <: AbstractObject{T, S} end

refractive_index(object::AbstractRefractiveOptic) = object.n
refractive_index(object::AbstractRefractiveOptic{<:Any, <:Any, <:RefractiveIndex}, λ::Real)::Float64 = object.n(λ)

"""
    interact3d(AbstractSystem, AbstractRefractiveOptic, Beam, Ray)

Implements the refraction of a [`Ray`](@ref) at an optical surface. The "outside" ref. index is obtained from the `system` unless specified otherwise.
At the critical angle, total internal reflection occurs (see [`refraction3d`](@ref)).
"""
function interact3d(system::AbstractSystem,
        optic::AbstractRefractiveOptic,
        ::Beam{T, R},
        ray::R) where {T <: Real, R <: Ray{T}}
    # Check dir. of ray and surface normal
    normal = normal3d(intersection(ray))
    lambda = wavelength(ray)
    if dot(direction(ray), normal) < 0
        # Entering optic
        n1 = refractive_index(ray)
        n2 = refractive_index(optic, lambda)
        # Hint to test optic again
        hint = Hint(optic)
    else
        # Exiting optic
        n1 = refractive_index(optic, lambda)
        n2 = refractive_index(system, lambda)
        hint = nothing
        # Flip normal for refraction3d
        normal = -normal
    end
    # Calculate new dir. and pos.
    ndir, TIR = refraction3d(direction(ray), normal, n1, n2)
    npos = position(ray) + length(ray) * direction(ray)
    # In case of TIR, update hint and n2
    if TIR
        hint = Hint(optic)
        n2 = refractive_index(optic, lambda)
    end
    return BeamInteraction{T, R}(hint,
        Ray{T}(npos, ndir, nothing, wavelength(ray), n2))
end

"""
    interact3d(AbstractSystem, AbstractRefractiveOptic, Beam, PolarizedRay)

Implements the refraction of a [`PolarizedRay`](@ref) at an uncoated optical surface. The "outside" ref. index is obtained from the `system` unless specified otherwise.
Reflection and transmission values are calculated via the [`fresnel_coefficients`](@ref). Stray light is not tracked.
In the case of total internal reflection, only the reflected light is traced.
"""
function interact3d(system::AbstractSystem, optic::AbstractRefractiveOptic,
        ::Beam{T, R}, ray::R) where {T <: Real, R <: PolarizedRay{T}}
    lambda = wavelength(ray)
    normal = normal3d(intersection(ray))
    raypos = position(ray) + length(ray) * direction(ray)
    if dot(direction(ray), normal) < 0
        # Entering optic
        n1 = refractive_index(ray)
        n2 = refractive_index(optic, lambda)
        # Hint to test optic again
        hint = Hint(optic)
    else
        # Exiting optic
        n1 = refractive_index(optic, lambda)
        n2 = refractive_index(system, lambda)
        hint = nothing
        # Flip normal for refraction3d
        normal = -normal
    end
    # Calculate (and correct into 1. quadrant) the angle of incidence
    θi = angle3d(direction(ray), -normal)
    # Get Fresnel coefficients
    rs, rp, ts, tp = fresnel_coefficients(θi, n2 / n1)
    # Optical interaction
    if is_internally_reflected(rp, rs)
        # Update hint and outgoing ref. index
        hint = Hint(optic)
        n2 = refractive_index(optic, lambda)
        # Calculate reflection
        new_dir = reflection3d(direction(ray), normal)
        J = [-rs 0 0; 0 rp 0; 0 0 1]
    else
        # Calculate refraction
        new_dir, ~ = refraction3d(direction(ray), normal, n1, n2)
        J = [ts 0 0; 0 tp 0; 0 0 1]
    end
    # Calculate new polarization
    E0 = _calculate_global_E0(direction(ray), new_dir, J, polarization(ray))
    return BeamInteraction{T, R}(
        hint, PolarizedRay{T}(raypos, new_dir, nothing, wavelength(ray), n2, E0))
end

"""
    Lens{T, S <: AbstractShape{T}, N <: RefractiveIndex} <: AbstractRefractiveOptic{T, S, N}

Represents an uncoated `Lens` with a homogeneous [`RefractiveIndex`](@ref) `n = n(λ)`.
Refer to the [`SphericalLens`](@ref) constructor for more information on how to generate lenses.

# Fields

- `shape`: geometry of the lens, refer to [`AbstractShape`](@ref) for more information
- `n`: [`RefractiveIndex`](@ref) function that returns n(λ)

# Additional information

!!! info "Refractive index"
    The chromatic dispersion of the lens is represented by a λ-dependent function for `n`
    and must be provided by the user. For testing purposes, an anonymous function, e.g. λ -> 1.5
    can be passed such that the lens has the same refractive index for all wavelengths.
"""
struct Lens{T, S <: AbstractShape{T}, N <: RefractiveIndex} <: AbstractRefractiveOptic{T, S, N}
    shape::S
    n::N
    function Lens(shape::S, n::N) where {T<:Real, S<:AbstractShape{T}, N<:RefractiveIndex}
        test_refractive_index_function(n)
        return new{T, S, N}(shape, n)
    end
end

thickness(l::Lens) = thickness(shape(l))