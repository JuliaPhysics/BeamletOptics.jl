"""
    ThinBeamsplitter <: AbstractBeamsplitter

Represents a 2D beam-splitting device. 

# Fields

- `shape`: 2D [`AbstractShape`](@ref) at which the splitting process occurs (e.g. a 2D-[`Mesh`](@ref))
- `reflectance`: scalar reflection factor
- `transmittance`: scalar transmission factor

!!! warning
    Note that the `transmittance` should be calculated from an input `reflectance` in order
    to ensure that R² + T² = 1.
"""
struct ThinBeamsplitter{T <: Real, S <: AbstractShape{T}} <: AbstractBeamsplitter{T, S}
    shape::S
    reflectance::T
    transmittance::T
end

reflectance(bs::ThinBeamsplitter) = bs.reflectance
transmittance(bs::ThinBeamsplitter) = bs.transmittance

"""
    ThinBeamsplitter(width, height; reflectance=0.5)

Creates a zero-thickness, lossless, non-polarizing 2D rectangular [`ThinBeamsplitter`](@ref) where

- `width`: is the x-dir. edge length in [m]
- `height`: is the z-dir. edge length in [m]
- `reflectance`: kw-arg that determines how much light is **reflected**, i.e. 0.7 for a 70:30 splitter

# Additional information

!!! info "Reflectance"
    The input value for the `reflectance` R is normed such that R² + T² = 1, where T is the `transmittance`.
    The transmittance is calculated via T = √(1 - R²).

!!! warning "Reflection phase jump"
    Note that the reflection phase jump θᵣ is implemented by the individual [`interact3d`](@ref)-methods. Refer to them for more information.
"""
function ThinBeamsplitter(width::Real, height::Real; reflectance::Real = 0.5)
    if reflectance ≥ 1 || reflectance ≈ 0
        error("Splitting ratio ∈ (0, 1)!")
    end
    shape = RectangularFlatMesh(width, height)
    Reflected = sqrt(reflectance)
    Transmitted = sqrt(1 - Reflected^2)
    return ThinBeamsplitter(shape, Reflected, Transmitted)
end

"""
    RoundThinBeamsplitter(diameter; reflectance=0.5)

Creates a zero-thickness, 2D round [`ThinBeamsplitter`](@ref) with the specified `diameter` in [m].
For more information, refer to the [`ThinBeamsplitter`](@ref) constructor.
"""
function RoundThinBeamsplitter(diameter::Real; reflectance::Real = 0.5)
    if reflectance ≥ 1 || reflectance ≈ 0
        error("Splitting ratio ∈ (0, 1)!")
    end
    shape = CircularFlatMesh(diameter/2)
    R = sqrt(reflectance)
    T = sqrt(1 - R^2)
    return ThinBeamsplitter(shape, R, T)
end

ThinBeamsplitter(width::Real; reflectance::Real=0.5) = ThinBeamsplitter(width, width; reflectance)

Base.isvalid(bs::AbstractBeamsplitter) = reflectance(bs)^2 + transmittance(bs)^2 ≈ 1

@inline function _beamsplitter_transmitted_beam(
        ::AbstractBeamsplitter, ::Beam{T, R}, ray::R) where {T <: Real, R <: Ray{T}}
    pos = position(ray) + length(ray) * direction(ray)
    dir = direction(ray)
    return Beam(Ray(pos, dir, wavelength(ray)))
end

@inline function _beamsplitter_reflected_beam(
        ::AbstractBeamsplitter, ::Beam{T, R}, ray::R) where {T <: Real, R <: Ray{T}}
    normal = normal3d(intersection(ray))
    pos = position(ray) + length(ray) * direction(ray)
    dir = reflection3d(direction(ray), normal)
    return Beam(Ray(pos, dir, wavelength(ray)))
end

@inline function _beamsplitter_transmitted_beam(bs::AbstractBeamsplitter, ::Beam{T, R},
        ray::R) where {T <: Real, R <: PolarizedRay{T}}
    J = SPBasis([transmittance(bs) 0 0; 0 transmittance(bs) 0; 0 0 1])
    pos = position(ray) + length(ray) * direction(ray)
    dir = direction(ray)
    E0 = _calculate_global_E0(bs, ray, dir, J)
    return Beam(PolarizedRay(pos, dir, wavelength(ray), E0))
end

@inline function _beamsplitter_reflected_beam(bs::AbstractBeamsplitter, ::Beam{T, R},
    ray::R) where {T <: Real, R <: PolarizedRay{T}}
    J = SPBasis([-reflectance(bs) 0 0; 0 reflectance(bs) 0; 0 0 1])
    normal = normal3d(intersection(ray))
    pos = position(ray) + length(ray) * direction(ray)
    in_dir = direction(ray)
    out_dir = reflection3d(in_dir, normal)
    E0 = _calculate_global_E0(bs, ray, out_dir, J)
    return Beam(PolarizedRay(pos, out_dir, wavelength(ray), E0))
end

function interact3d(::AbstractSystem, bs::ThinBeamsplitter, beam::Beam{T, R},
        ray::R) where {T <: Real, R <: AbstractRay{T}}
    # Push transmitted and reflected beams to system
    children!(beam,
        [_beamsplitter_transmitted_beam(bs, beam, ray), _beamsplitter_reflected_beam(bs, beam, ray)])
    # Stop for beam spawning
    return nothing
end

@inline function _beamsplitter_transmitted_beam(bs::AbstractBeamsplitter, gauss::GaussianBeamlet, ray_id::Int)
    # Transmitted Gaussian (no phase flip)
    chief = _beamsplitter_transmitted_beam(bs, gauss.chief, rays(gauss.chief)[ray_id])
    waist = _beamsplitter_transmitted_beam(bs, gauss.waist, rays(gauss.waist)[ray_id])
    divergence = _beamsplitter_transmitted_beam(bs, gauss.divergence, rays(gauss.divergence)[ray_id])
    λ = wavelength(gauss)
    w0 = gauss_parameters(gauss, length(gauss))[4]
    E0 = transmittance(bs) * electric_field(gauss) * (beam_waist(gauss) / w0)
    return GaussianBeamlet(chief, waist, divergence, λ, w0, E0)
end

@inline function _beamsplitter_reflected_beam(bs::AbstractBeamsplitter, gauss::GaussianBeamlet, ray_id::Int)
    # Reflected Gaussian (no phase flip)
    chief = _beamsplitter_reflected_beam(bs, gauss.chief, rays(gauss.chief)[ray_id])
    waist = _beamsplitter_reflected_beam(bs, gauss.waist, rays(gauss.waist)[ray_id])
    divergence = _beamsplitter_reflected_beam(bs, gauss.divergence, rays(gauss.divergence)[ray_id])
    λ = wavelength(gauss)
    w0 = gauss_parameters(gauss, length(gauss))[4]
    E0 = reflectance(bs) * electric_field(gauss) * (beam_waist(gauss) / w0)
    return GaussianBeamlet(chief, waist, divergence, λ, w0, E0)
end

"""
    interact3d(::AbstractSystem, bs::ThinBeamsplitter, gauss::GaussianBeamlet, ray_id::Int)

Models the interaction between a [`ThinBeamsplitter`](@ref) and a [`GaussianBeamlet`](@ref).

# Reflection phase jump

The reflection phase jump is modeled here as θᵣ = π for simplicity. This is since in practice it will have only a relative effect on the signal at the detector for interferometric setups.
The phase jump is applied to the reflected portion of any incoming beam that faces the [`ThinBeamsplitter`](@ref) normal vector, which assumes that the splitter has an unambigous normal, i.e. a 2D mesh.
This is intended to model the effect of the Fresnel equations without full polarization calculus.
"""
function interact3d(::AbstractSystem, bs::ThinBeamsplitter, gauss::GaussianBeamlet, ray_id::Int)
    # Phase flip
    ray = gauss.chief.rays[ray_id]
    df = dot(direction(ray), normal3d(intersection(ray)))
    if df < 0
        ϕ = π
    else
        ϕ = 0
    end

    t = _beamsplitter_transmitted_beam(bs, gauss, ray_id)
    r = _beamsplitter_reflected_beam(bs, gauss, ray_id)

    # Add conditional phase flip to reflected beam
    r.E0 *= exp(im*ϕ)

    children!(gauss, [t, r])
    return nothing
end
