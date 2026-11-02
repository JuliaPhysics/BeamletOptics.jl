"""
    BeamsplitterCoating{T} <: AbstractSurfaceModel

A lossless, non-polarizing coating that splits a ray into a transmitted and a reflected beam with the
field factors `reflectance` and `transmittance` (R² + T² = 1). It acts between the media of the
boundary, see [`Transition`](@ref): the transmitted beam is refracted into the medium behind it, e.g.
the substrate that it is coated on, the reflected beam stays in the medium before it. Beyond the
critical angle of these media, the ray is totally reflected and not split.

A [`ThinBeamsplitter`](@ref) carries this coating, see [`surface_model`](@ref).
"""
struct BeamsplitterCoating{T <: Real} <: AbstractSurfaceModel
    reflectance::T
    transmittance::T
end

reflectance(coating::BeamsplitterCoating) = coating.reflectance
transmittance(coating::BeamsplitterCoating) = coating.transmittance

"""
    ThinBeamsplitter <: AbstractBeamsplitter

Represents a 2D beam-splitting device: a thin layer whose [`BeamsplitterCoating`](@ref) acts at the
boundary on which it lies, see [`surface_model`](@ref).

# Fields

- `shape`: 2D [`AbstractShape`](@ref) at which the splitting process occurs (e.g. a 2D-[`Mesh`](@ref))
- `coating`: the [`BeamsplitterCoating`](@ref) with the field factors `reflectance` and `transmittance`

!!! warning
    Note that the `transmittance` should be calculated from an input `reflectance` in order
    to ensure that R² + T² = 1.
"""
struct ThinBeamsplitter{T <: Real, S <: AbstractShape{T}} <: AbstractBeamsplitter{T}
    shape::S
    coating::BeamsplitterCoating{T}
end

# The factors in the number type of the `shape`, e.g. of a shape built from dual numbers for
# automatic differentiation w.r.t. its pose, while the reflectance is a `Float64`
function ThinBeamsplitter(shape::AbstractShape{T}, reflectance::Real, transmittance::Real) where {T}
    return ThinBeamsplitter{T, typeof(shape)}(shape, BeamsplitterCoating(T(reflectance), T(transmittance)))
end

surface_model(bs::ThinBeamsplitter) = bs.coating
reflectance(bs::ThinBeamsplitter) = reflectance(surface_model(bs))
transmittance(bs::ThinBeamsplitter) = transmittance(surface_model(bs))

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
    shape = CircularFlatMesh(diameter / 2)
    R = sqrt(reflectance)
    T = sqrt(1 - R^2)
    return ThinBeamsplitter(shape, R, T)
end

function ThinBeamsplitter(width::Real; reflectance::Real = 0.5)
    ThinBeamsplitter(width, width; reflectance)
end

Base.isvalid(bs::AbstractBeamsplitter) = reflectance(bs)^2 + transmittance(bs)^2 ≈ 1


@inline function _beamsplitter_transmitted_beam(
        ::BeamsplitterCoating, ::Beam{T, R}, ray::R) where {T <: Real, R <: Ray{T}}
    pos = hit_point(ray)
    dir = direction(ray)
    return Beam(Ray(pos, dir, wavelength(ray)))
end

@inline function _beamsplitter_reflected_beam(
        ::BeamsplitterCoating, ::Beam{T, R}, ray::R) where {T <: Real, R <: Ray{T}}
    normal = normal3d(intersection(ray))
    pos = hit_point(ray)
    dir = reflection3d(direction(ray), normal)
    return Beam(Ray(pos, dir, wavelength(ray)))
end

@inline function _beamsplitter_transmitted_beam(coating::BeamsplitterCoating, ::Beam{T, R},
        ray::R) where {T <: Real, R <: PolarizedRay{T}}
    J = SPBasis(transmittance(coating), 0, 0, transmittance(coating))
    pos = hit_point(ray)
    dir = direction(ray)
    E0 = _calculate_global_E0(coating, ray, dir, J)
    return Beam(PolarizedRay(pos, dir, wavelength(ray), E0))
end

@inline function _beamsplitter_reflected_beam(coating::BeamsplitterCoating, ::Beam{T, R},
        ray::R) where {T <: Real, R <: PolarizedRay{T}}
    J = SPBasis(-reflectance(coating), 0, 0, reflectance(coating))
    normal = normal3d(intersection(ray))
    pos = hit_point(ray)
    in_dir = direction(ray)
    out_dir = reflection3d(in_dir, normal)
    E0 = _calculate_global_E0(coating, ray, out_dir, J)
    return Beam(PolarizedRay(pos, out_dir, wavelength(ray), E0))
end

"""
    transmit!(child::AbstractRay, ray::AbstractRay, t::Transition)

Turns the first ray `child` of the transmitted beam of a coating, which continues along the incident
`ray`, into the ray behind the boundary of the transition `t`: refracted into the medium behind it. For a
[`PolarizedRay`](@ref), the field is transformed into the new direction as well, and scaled for the
medium behind the boundary: the coating transmits a fraction of the power, and the power of a beam is
proportional to `n cos(θ) |E|²`.
"""
function transmit!(child::AbstractRay, ray::AbstractRay, t::Transition)
    dir, _ = refraction3d(direction(ray), normal3d(t), refractive_index_in(t), refractive_index_out(t))
    transmit_field!(child, ray, dir, t)
    refractive_index!(child, refractive_index_out(t))
    direction!(child, dir)
    return nothing
end

"""
    transmit_field!(child, ray, dir, t::Transition)

Transforms the field of the transmitted `child` into the refracted direction `dir` and scales it for
the medium behind the boundary of `t`; a [`Ray`](@ref) carries no field.
"""
transmit_field!(::AbstractRay, ::AbstractRay, dir, ::Transition) = nothing

function transmit_field!(child::PolarizedRay, ray::PolarizedRay, dir, t::Transition)
    normal = normal3d(t)
    P = _calculate_global_E0(direction(ray), dir, normal3d(intersection(ray)), SPBasis(1, 0, 0, 1))
    cos_in = abs(cos_angle3d(direction(ray), normal))
    cos_out = abs(cos_angle3d(dir, normal))
    scale = sqrt(refractive_index_in(t) * cos_in / (refractive_index_out(t) * cos_out))
    polarization!(child, scale * _transverse(P * polarization(child), dir))
    return nothing
end

"""
    reflect!(child::AbstractRay, t::Transition)

Sets the medium of the first ray `child` of the reflected beam of a coating: the medium before the
boundary of the transition `t`.
"""
reflect!(child::AbstractRay, t::Transition) = refractive_index!(child, refractive_index_in(t))

"""
    interact3d(system, coating::BeamsplitterCoating, t::Transition, beam::Beam, ray)

Splits the `ray` into a transmitted and a reflected child beam, see [`BeamsplitterCoating`](@ref).
"""
function interact3d(::AbstractSystem, coating::BeamsplitterCoating, t::Transition, beam::Beam{T, R},
        ray::R) where {T <: Real, R <: AbstractRay{T}}
    is_totally_reflected(ray, t) && return cross_boundary(t, beam, ray)
    transmitted = _beamsplitter_transmitted_beam(coating, beam, ray)
    reflected = _beamsplitter_reflected_beam(coating, beam, ray)
    transmit!(first(rays(transmitted)), ray, t)
    reflect!(first(rays(reflected)), t)
    # Push transmitted and reflected beams to system
    children!(beam, [transmitted, reflected])
    # Stop for beam spawning
    return nothing
end

@inline function _beamsplitter_transmitted_beam(
        coating::BeamsplitterCoating, gauss::GaussianBeamlet, ray_id::Int)
    # Transmitted Gaussian (no phase flip)
    chief = _beamsplitter_transmitted_beam(coating, gauss.chief, rays(gauss.chief)[ray_id])
    waist = _beamsplitter_transmitted_beam(coating, gauss.waist, rays(gauss.waist)[ray_id])
    divergence = _beamsplitter_transmitted_beam(
        coating, gauss.divergence, rays(gauss.divergence)[ray_id])
    λ = wavelength(gauss)
    w0 = gauss_parameters(gauss, length(gauss))[4]
    E0 = transmittance(coating) * electric_field(gauss) * (beam_waist(gauss) / w0)
    return GaussianBeamlet(chief, waist, divergence, λ, w0, E0)
end

@inline function _beamsplitter_reflected_beam(
        coating::BeamsplitterCoating, gauss::GaussianBeamlet, ray_id::Int)
    # Reflected Gaussian (no phase flip)
    chief = _beamsplitter_reflected_beam(coating, gauss.chief, rays(gauss.chief)[ray_id])
    waist = _beamsplitter_reflected_beam(coating, gauss.waist, rays(gauss.waist)[ray_id])
    divergence = _beamsplitter_reflected_beam(
        coating, gauss.divergence, rays(gauss.divergence)[ray_id])
    λ = wavelength(gauss)
    w0 = gauss_parameters(gauss, length(gauss))[4]
    E0 = reflectance(coating) * electric_field(gauss) * (beam_waist(gauss) / w0)
    return GaussianBeamlet(chief, waist, divergence, λ, w0, E0)
end

@inline function _beamsplitter_transmitted_beam(
        coating::BeamsplitterCoating, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    return AstigmaticGaussianBeamlet(map(
        b -> _beamsplitter_transmitted_beam(coating, b, rays(b)[ray_id]), _component_beams(agb))...)
end

@inline function _beamsplitter_reflected_beam(
        coating::BeamsplitterCoating, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    return AstigmaticGaussianBeamlet(map(
        b -> _beamsplitter_reflected_beam(coating, b, rays(b)[ray_id]), _component_beams(agb))...)
end

"""
    interact3d(system, coating::BeamsplitterCoating, t::Transition, beamlet, ray_id::Int)

Splits a [`GaussianBeamlet`](@ref) or an [`AstigmaticGaussianBeamlet`](@ref) into a transmitted and a
reflected child. Each component ray is refracted into the medium behind the coating or stays in the
medium before it, at its own hit. Beyond the critical angle of the transition `t` of the chief ray, the
`beamlet` is totally reflected and not split.

# Reflection phase jump

The reflection phase jump is modeled here as θᵣ = π for simplicity. This is since in practice it will have only a relative effect on the signal at the detector for interferometric setups.
The phase jump is applied to the reflected portion of any incoming beam that faces the normal vector of the coating, which assumes that the splitter has an unambigous normal, i.e. a 2D mesh.
This is intended to model the effect of the Fresnel equations without full polarization calculus.
"""
function interact3d(system::AbstractSystem, coating::BeamsplitterCoating, t::Transition,
        beamlet::GaussianBeamlets, ray_id::Int)
    ray = rays(beamlet)[ray_id]
    is_totally_reflected(ray, t) && return interact3d(system, Uncoated(), t, beamlet, ray_id)
    tb = _beamsplitter_transmitted_beam(coating, beamlet, ray_id)
    rb = _beamsplitter_reflected_beam(coating, beamlet, ray_id)
    # Each component ray at its own hit, map unrolls the heterogeneous tuples
    map(_component_beams(beamlet), _component_beams(tb), _component_beams(rb)) do b, tc, rc
        ray_b = rays(b)[ray_id]
        t_b = transition(system, ray_b)
        transmit!(first(rays(tc)), ray_b, t_b)
        reflect!(first(rays(rc)), t_b)
    end
    # Conditional phase flip of the reflected beam
    isentering(ray) && phase_shift!(rb, π)
    children!(beamlet, [tb, rb])
    return nothing
end

"""
    phase_shift!(beamlet, ϕ)

Shifts the phase of the `beamlet` by `ϕ` in [rad]: the field of a [`GaussianBeamlet`](@ref), the
polarization of the chief ray of an [`AstigmaticGaussianBeamlet`](@ref).
"""
phase_shift!(g::GaussianBeamlet, ϕ) = electric_field!(g, electric_field(g) * exp(im * ϕ))

function phase_shift!(agb::AstigmaticGaussianBeamlet, ϕ)
    chief_ray = first(rays(agb.c))
    polarization!(chief_ray, polarization(chief_ray) * exp(im * ϕ))
    return nothing
end
