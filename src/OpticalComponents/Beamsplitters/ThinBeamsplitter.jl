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
struct ThinBeamsplitter{T <: Real, S <: AbstractShape{T}} <: AbstractBeamsplitter{T}
    shape::S
    reflectance::T
    transmittance::T
end

# The factors in the number type of the `shape`, e.g. of a shape built from dual numbers for
# automatic differentiation w.r.t. its pose, while the reflectance is a `Float64`
function ThinBeamsplitter(shape::AbstractShape{T}, reflectance::Real, transmittance::Real) where {T}
    return ThinBeamsplitter{T, typeof(shape)}(shape, T(reflectance), T(transmittance))
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
        ::AbstractBeamsplitter, ::Beam{T, R}, ray::R) where {T <: Real, R <: Ray{T}}
    pos = hit_point(ray)
    dir = direction(ray)
    return Beam(Ray(pos, dir, wavelength(ray)))
end

@inline function _beamsplitter_reflected_beam(
        ::AbstractBeamsplitter, ::Beam{T, R}, ray::R) where {T <: Real, R <: Ray{T}}
    normal = normal3d(intersection(ray))
    pos = hit_point(ray)
    dir = reflection3d(direction(ray), normal)
    return Beam(Ray(pos, dir, wavelength(ray)))
end

@inline function _beamsplitter_transmitted_beam(bs::AbstractBeamsplitter, ::Beam{T, R},
        ray::R) where {T <: Real, R <: PolarizedRay{T}}
    J = SPBasis(transmittance(bs), 0, 0, transmittance(bs))
    pos = hit_point(ray)
    dir = direction(ray)
    E0 = _calculate_global_E0(bs, ray, dir, J)
    return Beam(PolarizedRay(pos, dir, wavelength(ray), E0))
end

@inline function _beamsplitter_reflected_beam(bs::AbstractBeamsplitter, ::Beam{T, R},
        ray::R) where {T <: Real, R <: PolarizedRay{T}}
    J = SPBasis(-reflectance(bs), 0, 0, reflectance(bs))
    normal = normal3d(intersection(ray))
    pos = hit_point(ray)
    in_dir = direction(ray)
    out_dir = reflection3d(in_dir, normal)
    E0 = _calculate_global_E0(bs, ray, out_dir, J)
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
    refractive_index!(child, refractive_index_out(t))
    direction!(child, dir)
    return nothing
end

function transmit!(child::PolarizedRay, ray::PolarizedRay, t::Transition)
    dir, _ = refraction3d(direction(ray), normal3d(t), refractive_index_in(t), refractive_index_out(t))
    normal = normal3d(t)
    P = _calculate_global_E0(direction(ray), dir, normal3d(intersection(ray)), SPBasis(1, 0, 0, 1))
    cos_in = abs(dot(direction(ray), normal))
    cos_out = abs(dot(dir, normal))
    scale = sqrt(refractive_index_in(t) * cos_in / (refractive_index_out(t) * cos_out))
    polarization!(child, scale * _transverse(P * polarization(child), dir))
    refractive_index!(child, refractive_index_out(t))
    direction!(child, dir)
    return nothing
end

"""
    reflect!(child::AbstractRay, t::Transition)

Sets the medium of the first ray `child` of the reflected beam of a coating: the medium before the
boundary of the transition `t`.
"""
reflect!(child::AbstractRay, t::Transition) = refractive_index!(child, refractive_index_in(t))

"""
    interact3d(system::AbstractSystem, bs::ThinBeamsplitter, beam::Beam, ray)

Splits the `ray` into a transmitted and a reflected child beam. The media on both sides of the splitter
follow from the [`transition`](@ref) at the hit: the transmitted beam is refracted into the medium
behind the splitter, e.g. the substrate that it is coated on, the reflected beam stays in the medium
before it. Beyond the critical angle of these media, the `ray` is totally reflected and not split.
"""
interact3d(system::AbstractSystem, bs::ThinBeamsplitter, beam::Beam{T, R}, ray::R) where {T <: Real, R <: AbstractRay{T}} =
    interact3d(system, bs, transition(system, ray), beam, ray)

function interact3d(::AbstractSystem, bs::ThinBeamsplitter, t::Transition, beam::Beam{T, R},
        ray::R) where {T <: Real, R <: AbstractRay{T}}
    is_totally_reflected(ray, t) && return cross_boundary(t, beam, ray)
    transmitted = _beamsplitter_transmitted_beam(bs, beam, ray)
    reflected = _beamsplitter_reflected_beam(bs, beam, ray)
    transmit!(first(rays(transmitted)), ray, t)
    reflect!(first(rays(reflected)), t)
    # Push transmitted and reflected beams to system
    children!(beam, [transmitted, reflected])
    # Stop for beam spawning
    return nothing
end

@inline function _beamsplitter_transmitted_beam(
        bs::AbstractBeamsplitter, gauss::GaussianBeamlet, ray_id::Int)
    # Transmitted Gaussian (no phase flip)
    chief = _beamsplitter_transmitted_beam(bs, gauss.chief, rays(gauss.chief)[ray_id])
    waist = _beamsplitter_transmitted_beam(bs, gauss.waist, rays(gauss.waist)[ray_id])
    divergence = _beamsplitter_transmitted_beam(
        bs, gauss.divergence, rays(gauss.divergence)[ray_id])
    λ = wavelength(gauss)
    w0 = gauss_parameters(gauss, length(gauss))[4]
    E0 = transmittance(bs) * electric_field(gauss) * (beam_waist(gauss) / w0)
    return GaussianBeamlet(chief, waist, divergence, λ, w0, E0)
end

@inline function _beamsplitter_reflected_beam(
        bs::AbstractBeamsplitter, gauss::GaussianBeamlet, ray_id::Int)
    # Reflected Gaussian (no phase flip)
    chief = _beamsplitter_reflected_beam(bs, gauss.chief, rays(gauss.chief)[ray_id])
    waist = _beamsplitter_reflected_beam(bs, gauss.waist, rays(gauss.waist)[ray_id])
    divergence = _beamsplitter_reflected_beam(
        bs, gauss.divergence, rays(gauss.divergence)[ray_id])
    λ = wavelength(gauss)
    w0 = gauss_parameters(gauss, length(gauss))[4]
    E0 = reflectance(bs) * electric_field(gauss) * (beam_waist(gauss) / w0)
    return GaussianBeamlet(chief, waist, divergence, λ, w0, E0)
end

@inline function _beamsplitter_transmitted_beam(
        bs::AbstractBeamsplitter, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    c = _beamsplitter_transmitted_beam(bs, agb.c, rays(agb.c)[ray_id])
    wxp = _beamsplitter_transmitted_beam(bs, agb.wxp, rays(agb.wxp)[ray_id])
    wxm = _beamsplitter_transmitted_beam(bs, agb.wxm, rays(agb.wxm)[ray_id])
    wyp = _beamsplitter_transmitted_beam(bs, agb.wyp, rays(agb.wyp)[ray_id])
    wym = _beamsplitter_transmitted_beam(bs, agb.wym, rays(agb.wym)[ray_id])
    dxp = _beamsplitter_transmitted_beam(bs, agb.dxp, rays(agb.dxp)[ray_id])
    dxm = _beamsplitter_transmitted_beam(bs, agb.dxm, rays(agb.dxm)[ray_id])
    dyp = _beamsplitter_transmitted_beam(bs, agb.dyp, rays(agb.dyp)[ray_id])
    dym = _beamsplitter_transmitted_beam(bs, agb.dym, rays(agb.dym)[ray_id])
    return AstigmaticGaussianBeamlet(c, wxp, wxm, wyp, wym, dxp, dxm, dyp, dym)
end

@inline function _beamsplitter_reflected_beam(
        bs::AbstractBeamsplitter, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    c = _beamsplitter_reflected_beam(bs, agb.c, rays(agb.c)[ray_id])
    wxp = _beamsplitter_reflected_beam(bs, agb.wxp, rays(agb.wxp)[ray_id])
    wxm = _beamsplitter_reflected_beam(bs, agb.wxm, rays(agb.wxm)[ray_id])
    wyp = _beamsplitter_reflected_beam(bs, agb.wyp, rays(agb.wyp)[ray_id])
    wym = _beamsplitter_reflected_beam(bs, agb.wym, rays(agb.wym)[ray_id])
    dxp = _beamsplitter_reflected_beam(bs, agb.dxp, rays(agb.dxp)[ray_id])
    dxm = _beamsplitter_reflected_beam(bs, agb.dxm, rays(agb.dxm)[ray_id])
    dyp = _beamsplitter_reflected_beam(bs, agb.dyp, rays(agb.dyp)[ray_id])
    dym = _beamsplitter_reflected_beam(bs, agb.dym, rays(agb.dym)[ray_id])
    return AstigmaticGaussianBeamlet(c, wxp, wxm, wyp, wym, dxp, dxm, dyp, dym)
end

"""
    interact3d(::AbstractSystem, bs::ThinBeamsplitter, beamlet, ray_id::Int)

Models the interaction between a [`ThinBeamsplitter`](@ref) and a [`GaussianBeamlet`](@ref) or an
[`AstigmaticGaussianBeamlet`](@ref): it splits the `beamlet` into a transmitted and a reflected child.
Each component ray is refracted into the medium behind the splitter or stays in the medium before it,
see [`transition`](@ref). Beyond the critical angle of these media, the `beamlet` is totally reflected
and not split.

# Reflection phase jump

The reflection phase jump is modeled here as θᵣ = π for simplicity. This is since in practice it will have only a relative effect on the signal at the detector for interferometric setups.
The phase jump is applied to the reflected portion of any incoming beam that faces the [`ThinBeamsplitter`](@ref) normal vector, which assumes that the splitter has an unambigous normal, i.e. a 2D mesh.
This is intended to model the effect of the Fresnel equations without full polarization calculus.
"""
function interact3d(system::AbstractSystem, bs::ThinBeamsplitter, beamlet::GaussianBeamlets, ray_id::Int)
    ray = rays(_chief_beam(beamlet))[ray_id]
    is_totally_reflected(ray, transition(system, ray)) && return cross_boundary(system, beamlet, ray_id)
    t = _beamsplitter_transmitted_beam(bs, beamlet, ray_id)
    r = _beamsplitter_reflected_beam(bs, beamlet, ray_id)
    # Each component ray crosses the boundary at its own hit. The component beams differ in their ray
    # type, hence `map` over the tuples, which is unrolled
    map(_component_beams(beamlet), _component_beams(t), _component_beams(r)) do b, tb, rb
        ray_b = rays(b)[ray_id]
        tr = transition(system, ray_b)
        transmit!(first(rays(tb)), ray_b, tr)
        reflect!(first(rays(rb)), tr)
    end
    # Add conditional phase flip to the reflected beam
    dot(direction(ray), normal3d(intersection(ray))) < 0 && phase_shift!(r, π)
    children!(beamlet, [t, r])
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
