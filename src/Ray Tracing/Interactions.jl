
"""
    AbstractReflectiveOptic <: AbstractObject

A generic type to represent `AbstractObject`s which reflect incoming rays. The main function of `interact3d` should be akin to [`reflection3d`](@ref).
"""
abstract type AbstractReflectiveOptic <: AbstractObject end

"""
    interact3d(AbstractReflectiveOptic, Ray)

Implements the reflection of a [`Ray`](@ref) via the normal at the intersection point on an optical surface.
"""
function interact3d(::AbstractSystem,
        ::AbstractReflectiveOptic,
        ::Beam{T, R},
        ray::R) where {T <: Real, R <: Ray{T}}
    normal = intersection(ray).n
    npos = position(ray) + length(ray) * direction(ray)
    ndir = reflection3d(direction(ray), normal)
    return BeamInteraction{T, R}(nothing, Ray{T}(uuid4(), npos, ndir, nothing, wavelength(ray), refractive_index(ray)))
end

"""
    interact3d(AbstractReflectiveOptic, PolarizedRay)

Implements the ideal reflection of a [`PolarizedRay`](@ref) via the normal at the intersection point on an optical surface.
A Jones matrix of [-1 0 0; 0 1 0] is assumed as per Peatross (2015, 2023 Ed. p. 154) and Yun et al. (see [`PolarizedRay`](@ref) for more information).
"""
function interact3d(::AbstractSystem,
    ::AbstractReflectiveOptic,
    ::Beam{T, R},
    ray::R) where {T <: Real, R <: PolarizedRay{T}}
    normal = normal3d(intersection(ray))
    npos = position(ray) + length(ray) * direction(ray)
    ndir = reflection3d(direction(ray), normal)
    # Jones reflection matrix
    J = @SArray [-1 0 0; 0 1 0; 0 0 1]
    E0 = _calculate_global_E0(direction(ray), ndir, J, electric_field(ray))
    return BeamInteraction{T, R}(nothing, PolarizedRay{T}(uuid4(), npos, ndir, nothing, wavelength(ray), refractive_index(ray), E0))
end

"""
    Mirror{S <: AbstractShape} <: AbstractReflectiveOptic

Concrete implementation of a perfect mirror with arbitrary shape.
"""
struct Mirror{S <: AbstractShape} <: AbstractReflectiveOptic
    id::UUID
    shape::S
end

"""
    AbstractRefractiveOptic <: AbstractObject

A generic type to represent `AbstractObject`s which refract incoming rays. The main function of `interact3d` should be akin to [`refraction3d`](@ref).

# Implementation reqs.

Subtypes of `AbstractRefractiveOptic` should implement all supertype reqs. as well as:

# Fields

- `n::Function`: a function which returns the refractive index for a wavelength λ
"""
abstract type AbstractRefractiveOptic{T} <: AbstractObject end

refractive_index(object::AbstractRefractiveOptic) = object.n
refractive_index(object::AbstractRefractiveOptic{<: Function}, λ::Real)::Float64 = object.n(λ)

"""
    interact3d(AbstractRefractiveOptic, Ray)

Implements the refraction of a [`Ray`](@ref) at an optical surface. The "outside" ref. index is always assumed to be 1.
At the critical angle, total internal reflection occurs (see [`refraction3d`](@ref)).
"""
function interact3d(::AbstractSystem,
        object::AbstractRefractiveOptic,
        ::Beam{T, R},
        ray::R) where {T <: Real, R <: Ray{T}}
    # Check dir. of ray and surface normal
    normal = normal3d(intersection(ray))
    λ = wavelength(ray)
    if dot(direction(ray), normal) < 0
        # "Outside prism"
        n1 = 1.0
        n2 = refractive_index(object, λ)
    else
        # "Inside prism"
        n1 = refractive_index(object, λ)
        n2 = 1.0
        normal *= -1
    end
    # Calculate new dir. and pos.
    ndir = refraction3d(direction(ray), normal, n1, n2) # FIXME: in case of TIR, n2 is not correctly set to object ref. index
    npos = position(ray) + length(ray) * direction(ray)
    # Hint is the current object ID
    return BeamInteraction{T, R}(id(object),
        Ray{T}(uuid4(), npos, ndir, nothing, wavelength(ray), n2))
end

"""
    interact3d(AbstractRefractiveOptic, PolarizedRay)

Implements the refraction of a [`PolarizedRay`](@ref) at an uncoated optical surface. The "outside" ref. index is always assumed to be 1.
Reflection and transmission values are calculated via the [`fresnel_coefficients`](@ref). Stray light is not tracked.
In the case of total internal reflection, only the reflected light is traced.
"""
function interact3d(system::AbstractSystem, optic::AbstractRefractiveOptic, ::Beam{T,R}, ray::R) where {T<:Real,R<:PolarizedRay{T}}
    lambda = wavelength(ray)
    normal = normal3d(intersection(ray))
    raypos = position(ray) + length(ray) * direction(ray)
    if dot(direction(ray), normal) < 0
        # Entering optic
        n1 = refractive_index(ray)
        n2 = refractive_index(optic, lambda)
        # Hint to test optic again
        hint = id(optic)
    else
        # Exiting optic
        n1 = refractive_index(optic, lambda)
        n2 = refractive_index(system)
        hint = nothing
        # Flip normal for refraction3d
        normal = -normal
    end
    # Calculate (and correct into 1. quadrant) the angle of incidence
    θi = angle3d(direction(ray), -normal)
    # Get Fresnel coefficients
    rs, rp, ts, tp = fresnel_coefficients(θi, n2/n1)
    # Optical interaction
    if is_internally_reflected(rp, rs)
        # Update hint and outgoing ref. index
        hint = id(optic)
        n2 = refractive_index(optic, lambda)
        # Calculate reflection
        new_dir = reflection3d(direction(ray), normal)
        J = [-rs 0 0; 0 rp 0; 0 0 1]
    else
        # Calculate refraction
        new_dir = refraction3d(direction(ray), normal, n1, n2)
        J = [ts 0 0; 0 tp 0; 0 0 1]
    end
    # Calculate new polarization
    E0 = _calculate_global_E0(direction(ray), new_dir, J, electric_field(ray))
    return BeamInteraction{T, R}(hint, PolarizedRay{T}(uuid4(), raypos, new_dir, nothing, wavelength(ray), n2, E0))
end

"""
    Lens{S <: AbstractShape, T <: Function} <: AbstractRefractiveOptic{T}

Represents an uncoated `Lens` with a homogeneous refractive index `n = n(λ)`.
Refer to the [`SphericalLens`](@ref) constructor for more information on how to generate lenses.

# Fields

- `id`: lens ID (uuid4)
- `shape`: geometry of the lens, refer to [`AbstractShape`](@ref) for more information
- `n`: **single-argument** function that returns n(λ)

# Additional information

!!! info "Refractive index"
    The chromatic dispersion of the lens is represented by a λ-dependent function for `n`
    and must be provided by the user. For testing purposes, an anonymous function, e.g. λ -> 1.5
    can be passed such that the lens has the same refractive index for all wavelengths.
"""
struct Lens{S <: AbstractShape, T <: Function} <: AbstractRefractiveOptic{T}
    id::UUID
    shape::S
    n::T # FIXME: constructor check that n=n(λ), # args
end

SphericalLens(r1, r2, l, d, n) = SphericalLens(r1, r2, l, d, λ -> n)

"""
    SphericalLens(r1, r2, l, d=1inch, n=λ->1.5)

Creates a spherical [`Lens`](@ref) based on:

- `r1`: front radius
- `r2`: back radius
- `l`: lens thickness
- `d`: lens diameter, default is one inch
- `n`: refractive index as a function of λ

# Surface types

- `r1/2 > 0`: convex
- `r1/2 < 0`: concave
- `r2 = Inf`: planar

# Hint: thin lenses

If `l` is set to zero, an ideal thin [`SphericalLens`](@ref) will be created. However, note that the actual lens thickness will be different from zero.
"""
function SphericalLens(r1::Real, r2::Real, l::Real, d::Real = 1inch, n::Function = λ -> 1.5)
    # # Test for thin lens
    if iszero(l)
        shape = ThinLensSDF(r1, r2, d)
    elseif isinf(r1) && isinf(r2)# Test for cylinder lens. FIXME: This is incomplete!
        shape = CylinderSDF(d / 2, l / 2)
    elseif isinf(r2) # Test for plano lens
        if r1 > 0
            shape = PlanoConvexLensSDF(r1, l, d)
        elseif r1 < 0
            shape = PlanoConcaveLensSDF(abs(r1), l, d)
        end
    elseif r1 > 0 && r2 > 0 # Test for bi-convex/concave or meniscus
        shape = BiConvexLensSDF(r1, r2, l, d)
    elseif r1 < 0 && r2 < 0
        shape = BiConcaveLensSDF(abs(r1), abs(r2), l, d)
    elseif r1 > 0 && r2 < 0
        shape = ConvexConcaveLensSDF(r1, abs(r2), l, d)
    else
        throw(DomainError("Could not find suitable lens SDF for the given parameters"))
    end
    # Create lens
    return Lens(uuid4(), shape, n)
end

"""
    ThinLens(R1::Real, R2::Real, d::Real, n::Function)

Directly creates an ideal spherical thin [`Lens`](@ref) with radii of curvature `R1` and `R2` and diameter `d`
and refractive index `n`.
"""
function ThinLens(R1::Real, R2::Real, d::Real, n::Function)
    shape = ThinLensSDF(R1, R2, d)
    return Lens(uuid4(), shape, n)
end
ThinLens(R1::Real, R2::Real, d::Real, n::Real) = ThinLens(R1, R2, d, x -> n)

"""
    Prism{S <: AbstractShape, T <: Function} <: AbstractRefractiveOptic{T}

Essentially represents the same functionality as [`Lens`](@ref).
Refer to its documentation.
"""
struct Prism{S <: AbstractShape, T <: Function} <: AbstractRefractiveOptic{T}
    id::UUID
    shape::S
    n::T
end

#=
Asphere related constructors
=#

"""
    PlanoConvexAsphericalLens(radius, conic_constant, even_coefficients, d, t, n, md = d)

Creates a plano-aspherical [`Lens`](@ref) with convex shape based on:

- `radius`: front radius
- `conic_constant`: Conic constant of the aspherical surface
- `even_coefficients`: A vector of the even aspheric coefficients
- `d`: lens diameter
- `t`: lens thickness
- `n`: refractive index as a function of λ

!!! note
    Aspheric lenses are somewhat experimental at the moment. Use this feature with some caution.
    Future versions of this package will offer a convenience constructor for abitrary mixed
    lenses, e.g. bi-aspheres, aspheric-spheric lenses, etc. as well as odd aspheres and extended
    aspheres.
"""
function PlanoConvexAsphericalLens(radius::Real, conic_constant::Real, even_coefficients::Vector{<:Real}, d::Real, t::Real, n::Real)
    return PlanoConvexAsphericalLens(radius, conic_constant, even_coefficients, d, t, x -> n)
end

function PlanoConvexAsphericalLens(radius::Real, conic_constant::Real, even_coefficients::Vector{<:Real}, d::Real, t::Real, n::Function)
    shape = PlanoConvexAsphericalLensSDF(radius, t, d, conic_constant, even_coefficients)

    return Lens(uuid4(), shape, n)
end

"""
    PlanoConcaveAsphericalLens(radius, conic_constant, even_coefficients, d, t, n, md = d)

Creates a plano-aspherical [`Lens`](@ref) with convex shape based on:

- `radius`: front radius
- `conic_constant`: Conic constant of the aspherical surface
- `even_coefficients`: A vector of the even aspheric coefficients
- `d`: lens diameter
- `t`: lens thickness
- `n`: refractive index as a function of λ
- `md`: mechanical diameter of the lens (defaults to `d`). If this is set to a value `md` > `d`
        an outer flat section will be added to the lens. This can be used to model more realistic
        lenses where this flat section is present for mounting purposes but is also an active
        optical region for extreme lenses (HUDs, AR-wearables, etc.)

!!! note
    Aspheric lenses are somewhat experimental at the moment. Use this feature with some caution.
    Future versions of this package will offer a convenience constructor for abitrary mixed
    lenses, e.g. bi-aspheres, aspheric-spheric lenses, etc. as well as odd aspheres and extended
    aspheres.
"""
function PlanoConcaveAsphericalLens(radius::Real, conic_constant::Real, even_coefficients::Vector{<:Real}, d::Real, t::Real, n::Real, md::Real = d)
    return PlanoConcaveAsphericalLens(radius, conic_constant, even_coefficients, d, t, x -> n, md)
end

function PlanoConcaveAsphericalLens(radius::Real, conic_constant::Real, even_coefficients::Vector{<:Real}, d::Real, t::Real, n::Function, md::Real = d)
    shape = PlanoConcaveAsphericalLensSDF(radius, t, d, conic_constant, even_coefficients, md)

    return Lens(uuid4(), shape, n)
end


#=
Implements photodetector, efield calculation during solve_system!
=#
abstract type AbstractDetector <: AbstractObject end

"""
    Photodetector{S <: AbstractShape, T} <: AbstractDetector

Represents a **flat** rectangular or quadratic surface in R³ that is the active surface of a photodetector.
The active surface is discretized in the local R² x-y-coordinate system.
Field contributions Eᵢ are added by the corresponding [`interact3d`](@ref) method.

# Fields

- `id`: detector ID (uuid4)
- `shape`: geometry of the active surface, must represent 2D-`field` in `x` any `y` dimensions
- `x`: linear range of local x-coordinates
- `y`: linear range of local y-coordinates
- `field`: `size(x)` by `size(y)` matrix of complex values to store superposition E₀

# Additional information

!!! warning "Reset behavior"
    The `Photodetector` must be reset between each call of [`solve_system!`](@ref) in order to
    overwrite previous results using the [`reset_photodetector!`](@ref) function.
    Otherwise, the current result will be added onto the previous result.

!!! info "Supported beams"
    Currently, only the [`GaussianBeamlet`] is supported.
"""
mutable struct Photodetector{S <: AbstractShape, T} <: AbstractDetector
    const id::UUID
    const shape::S
    x::LinRange{T, Int64}
    y::LinRange{T, Int64}
    field::Matrix{Complex{T}}
end

function Photodetector(scale::T, n::Int) where {T}
    sz = 0.5
    vertices = [sz 0 sz
        sz 0 -sz
        -sz 0 -sz
        -sz 0 sz]
    faces = [1 2 4
        2 3 4]
    x = y = LinRange(-sz, sz, n) * scale
    field = zeros(Complex{T}, n, n)
    shape = Mesh{T}(uuid4(),
        vertices .* scale,
        faces,
        Matrix{T}(I, 3, 3),
        T.([0, 0, 0]),
        scale)
    return Photodetector{typeof(shape), T}(uuid4(), shape, x, y, field)
end

function interact3d(::AbstractSystem, ::Photodetector, ::B, ::Ray) where {B <: AbstractBeam}
    @warn "Photodetection for $B not implemented"
    return nothing
end

function interact3d(::AbstractSystem,
        pd::Photodetector,
        gauss::GaussianBeamlet,
        ray_id::Int)
    # Add efield contribution to pd.field
    ray = gauss.chief.rays[ray_id]
    l = length(gauss, opl=true)
    T = transpose(orientation(shape(pd)))
    P = ray.pos + ray.dir * length(ray.intersection)
    shape_pos = position(shape(pd))
    Threads.@threads for j in eachindex(pd.y)     # row column major order?
        y = pd.y[j]
        @inbounds for i in eachindex(pd.x)
            x = pd.x[i]
            # Transform point p on PD into world coords
            p = Point3(T[1, 1] * x + T[1, 3] * y + shape_pos[1],
                T[2, 1] * x + T[2, 3] * y + shape_pos[2],
                T[3, 1] * x + T[3, 3] * y + shape_pos[3])
            r = line_point_distance3d(ray, p)
            c = sqrt(sum(x -> (x[1] - x[2])^2, zip(P, p)))
            z = sqrt(abs(c^2 - r^2)) # abs to protect against small neg. values
            # Correct sign of z
            if isinfrontof(p, P, ray.dir)
                z = -z
            end
            pd.field[i, j] += electric_field(gauss, r, l + z)
        end
    end
    return nothing
end

intensity(pd::Photodetector) = intensity.(pd.field)

"""
    optical_power(pd::Photodetector)

Calculates the total optical power on `pd` in [W] by integration over the local intensity.
"""
optical_power(pd::Photodetector) = trapz((pd.x, pd.y), intensity(pd))

"""Resets the values currently stored in `pd.field` to zero"""
reset_photodetector!(pd::Photodetector{S, T}) where {S, T} = (pd.field .= zero(Complex{T}))

"""
    photodetector_resolution!(pd::Photodetector, n::Int)

Sets the resolution of `pd` to `n` × `n`. Note that this resets the current `pd.field`.
"""
function photodetector_resolution!(pd::Photodetector{S, T}, n::Int) where {S, T}
    pd.x = LinRange(pd.x.start, pd.x.stop, n)
    pd.y = LinRange(pd.y.start, pd.y.stop, n)
    pd.field = zeros(Complex{T}, n, n)
    return nothing
end

#=
Implements thin beam splitter, beam spawning
=#
abstract type AbstractBeamSplitter <: AbstractObject end

"""Models a generic beam splitter"""
struct BeamSplitter{S <: AbstractShape, T <: Real} <: AbstractBeamSplitter
    id::UUID
    shape::S
    reflectance::T
    transmittance::T
end

reflectance(bs::BeamSplitter) = bs.reflectance
transmittance(bs::BeamSplitter) = bs.transmittance

"""
    ThinBeamSplitter(scale::T, reflectance::Real=0.5) where {T}

Creates a zero-thickness, lossless, non-polarizing quadratic rectangle beam splitter where

- `scale`: is the edge length
- `reflectance`: determines how much light is **reflected**, i.e. 0.7 for a 70:30 splitter

## Reflectance

The input value for the `reflectance` R is normed such that R² + T² = 1, where T is the `transmittance`.
The transmittance is calculated via T = √(1 - R²).

## Phase shift

Note that the reflection phase shift θᵣ ∈ [0, π] is not modeled here for simplicity, since in practice it will have no effect on the interference at the detector.
"""
function ThinBeamSplitter(scale::T, reflectance::Real = 0.5) where {T}
    if reflectance ≥ 1 || isapprox(reflectance, 0)
        error("Splitting ratio ∈ (0, 1)!")
    end
    sz = 0.5
    vertices = [sz 0 sz
        sz 0 -sz
        -sz 0 -sz
        -sz 0 sz]
    faces = [1 2 4
        2 3 4]
    shape = Mesh{T}(uuid4(),
        vertices .* scale,
        faces,
        Matrix{T}(I, 3, 3),
        T.([0, 0, 0]),
        scale)
    Reflected = sqrt(reflectance)
    Transmitted = sqrt(1 - Reflected^2)
    return BeamSplitter(uuid4(), shape, Reflected, Transmitted)
end

Base.isvalid(bs::BeamSplitter) = reflectance(bs)^2 + transmittance(bs)^2 ≈ 1

@inline function _beamsplitter_transmitted_beam(::BeamSplitter, ray::Ray)
    pos = position(ray) + length(ray) * direction(ray)
    dir = direction(ray)
    return Beam(Ray(pos, dir, wavelength(ray)))
end

@inline function _beamsplitter_reflected_beam(::BeamSplitter, ray::Ray)
    normal = normal3d(intersection(ray))
    pos = position(ray) + length(ray) * direction(ray)
    dir = reflection3d(direction(ray), normal)
    return Beam(Ray(pos, dir, wavelength(ray)))
end

function _beamsplitter_transmitted_beam(bs::BeamSplitter, ray::PolarizedRay)
    J = @SArray [transmittance(bs) 0 0; 0 transmittance(bs) 0; 0 0 1]
    pos = position(ray) + length(ray) * direction(ray)
    dir = direction(ray)
    E0 =  _calculate_global_E0(dir, dir, J, electric_field(ray))
    return Beam(PolarizedRay(pos, dir, wavelength(ray), E0))
end

function _beamsplitter_reflected_beam(bs::BeamSplitter, ray::PolarizedRay)
    J = @SArray [-reflectance(bs) 0 0; 0 reflectance(bs) 0; 0 0 1]
    normal = normal3d(intersection(ray))
    pos = position(ray) + length(ray) * direction(ray)
    in_dir = direction(ray)
    out_dir = reflection3d(in_dir, normal)
    E0 =  _calculate_global_E0(in_dir, out_dir, J, electric_field(ray))
    return Beam(PolarizedRay(pos, out_dir, wavelength(ray), E0))
end

function interact3d(::AbstractSystem, bs::BeamSplitter, beam::Beam{T, R}, ray::R) where {T<:Real, R<:AbstractRay{T}}
    # Push transmitted and reflected beams to system
    children!(beam,
        [_beamsplitter_transmitted_beam(bs, ray), _beamsplitter_reflected_beam(bs, ray)])
    # Stop for beam spawning
    return nothing
end

"""
    interact3d(::AbstractSystem, bs::BeamSplitter, gauss::GaussianBeamlet, ray_id::Int)

Models the interaction between a [`BeamSplitter`](@ref) and a [`GaussianBeamlet`](@ref).
For more information refer to [`ThinBeamSplitter`](@ref).
"""
function interact3d(::AbstractSystem, bs::BeamSplitter, gauss::GaussianBeamlet, ray_id::Int)
    # Transmitted gauss
    chief = _beamsplitter_transmitted_beam(bs, rays(gauss.chief)[ray_id])
    waist = _beamsplitter_transmitted_beam(bs, rays(gauss.waist)[ray_id])
    divergence = _beamsplitter_transmitted_beam(bs, rays(gauss.divergence)[ray_id])
    λ = wavelength(gauss)
    w0 = gauss_parameters(gauss, length(gauss))[4]
    E0 = transmittance(bs) * electric_field(gauss) * (beam_waist(gauss) / w0)
    t = GaussianBeamlet(chief, waist, divergence, λ, w0, E0)
    # Reflected gauss
    chief = _beamsplitter_reflected_beam(bs, rays(gauss.chief)[ray_id])
    waist = _beamsplitter_reflected_beam(bs, rays(gauss.waist)[ray_id])
    divergence = _beamsplitter_reflected_beam(bs, rays(gauss.divergence)[ray_id])
    λ = wavelength(gauss)
    w0 = gauss_parameters(gauss, length(gauss))[4]
    E0 = reflectance(bs) * electric_field(gauss) * (beam_waist(gauss) / w0)
    r = GaussianBeamlet(chief, waist, divergence, λ, w0, E0)
    children!(gauss, [t, r])
    return nothing
end
