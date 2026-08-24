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

@inline function _beamsplitter_transmitted_ray(
        bs::AbstractBeamsplitter, trans::Transition, int::AbstractIntersection{T}, ray::Ray{T}) where {T <: Real}
    λ = wavelength(ray)
    n1 = refractive_index(trans.medium_in, λ)
    n2 = refractive_index(trans.medium_out, λ)
    normal = dot(normal3d(int), direction(ray)) < 0 ? normal3d(int) : -normal3d(int)
    new_dir, _ = refraction3d(direction(ray), normal, n1, n2)
    pos = position(int)
    return Ray{T}(pos, new_dir, λ, n2)
end

@inline function _beamsplitter_reflected_ray(
        bs::AbstractBeamsplitter, trans::Transition, int::AbstractIntersection{T}, ray::Ray{T}) where {T <: Real}
    λ = wavelength(ray)
    n1 = refractive_index(trans.medium_in, λ)
    normal = dot(normal3d(int), direction(ray)) < 0 ? normal3d(int) : -normal3d(int)
    new_dir = reflection3d(direction(ray), normal)
    pos = position(int)
    return Ray{T}(pos, new_dir, λ, n1)
end

@inline function _beamsplitter_transmitted_ray(
        bs::AbstractBeamsplitter, trans::Transition, int::AbstractIntersection{T}, ray::PolarizedRay{T}) where {T <: Real}
    λ = wavelength(ray)
    n1 = refractive_index(trans.medium_in, λ)
    n2 = refractive_index(trans.medium_out, λ)
    normal = dot(normal3d(int), direction(ray)) < 0 ? normal3d(int) : -normal3d(int)
    new_dir, _ = refraction3d(direction(ray), normal, n1, n2)
    pos = position(int)
    J = SPBasis(transmittance(bs), 0, 0, transmittance(bs))
    P = _calculate_global_E0(direction(ray), new_dir, normal, J)
    new_E0 = P * polarization(ray)
    return PolarizedRay{T}(pos, new_dir, λ, n2, new_E0)
end

@inline function _beamsplitter_reflected_ray(
        bs::AbstractBeamsplitter, trans::Transition, int::AbstractIntersection{T}, ray::PolarizedRay{T}) where {T <: Real}
    λ = wavelength(ray)
    n1 = refractive_index(trans.medium_in, λ)
    normal = dot(normal3d(int), direction(ray)) < 0 ? normal3d(int) : -normal3d(int)
    new_dir = reflection3d(direction(ray), normal)
    pos = position(int)
    J = SPBasis(-reflectance(bs), 0, 0, reflectance(bs))
    P = _calculate_global_E0(direction(ray), new_dir, normal, J)
    new_E0 = P * polarization(ray)
    return PolarizedRay{T}(pos, new_dir, λ, n1, new_E0)
end

function interact3d(
        system::AbstractSystem,
        bs::ThinBeamsplitter,
        trans::Transition,
        int::AbstractIntersection{T},
        beam::Beam{T},
        ray::AbstractRay{T}) where {T <: Real}
    t_ray = _beamsplitter_transmitted_ray(bs, trans, int, ray)
    r_ray = _beamsplitter_reflected_ray(bs, trans, int, ray)
    children!(beam, [Beam(t_ray), Beam(r_ray)])
    return nothing
end

function interact3d(
        system::AbstractSystem,
        bs::ThinBeamsplitter,
        beam::Beam{T},
        ray::AbstractRay{T}) where {T <: Real}
    int = isempty(beam.segments) ? intersect3d(shape(bs), ray) : intersection(last(beam.segments))
    isnothing(int) && return nothing
    ambient = ambient_medium(system)
    trans = resolve_transition(ambient, ambient, ray, normal3d(int))
    return interact3d(system, bs, trans, int, beam, ray)
end

function interact3d(
        system::AbstractSystem,
        bs::ThinBeamsplitter,
        trans::Transition,
        int::AbstractIntersection,
        gauss::GaussianBeamlet{R},
        ray_id::Int) where {R}
    chief_ray = rays(gauss.chief)[ray_id]
    waist_ray = rays(gauss.waist)[ray_id]
    div_ray = rays(gauss.divergence)[ray_id]

    int_c = ray_id <= length(gauss.chief.segments) ? (isnothing(intersection(gauss.chief.segments[ray_id])) ? int : intersection(gauss.chief.segments[ray_id])) : int
    int_w = ray_id <= length(gauss.waist.segments) ? (isnothing(intersection(gauss.waist.segments[ray_id])) ? int : intersection(gauss.waist.segments[ray_id])) : int
    int_d = ray_id <= length(gauss.divergence.segments) ? (isnothing(intersection(gauss.divergence.segments[ray_id])) ? int : intersection(gauss.divergence.segments[ray_id])) : int

    t_c = _beamsplitter_transmitted_ray(bs, trans, int_c, chief_ray)
    t_w = _beamsplitter_transmitted_ray(bs, trans, int_w, waist_ray)
    t_d = _beamsplitter_transmitted_ray(bs, trans, int_d, div_ray)

    r_c = _beamsplitter_reflected_ray(bs, trans, int_c, chief_ray)
    r_w = _beamsplitter_reflected_ray(bs, trans, int_w, waist_ray)
    r_d = _beamsplitter_reflected_ray(bs, trans, int_d, div_ray)

    λ = wavelength(gauss)
    w0 = gauss_parameters(gauss, length(gauss))[4]
    E0_t = transmittance(bs) * electric_field(gauss) * (beam_waist(gauss) / w0)
    E0_r = reflectance(bs) * electric_field(gauss) * (beam_waist(gauss) / w0)

    # Reflection phase flip if facing normal
    df = dot(direction(chief_ray), normal3d(int_c))
    ϕ = df < 0 ? π : 0
    E0_r *= exp(im * ϕ)

    t_beam = GaussianBeamlet(Beam(t_c), Beam(t_w), Beam(t_d), λ, w0, E0_t)
    r_beam = GaussianBeamlet(Beam(r_c), Beam(r_w), Beam(r_d), λ, w0, E0_r)

    children!(gauss, [t_beam, r_beam])
    return nothing
end

function interact3d(
        system::AbstractSystem,
        bs::ThinBeamsplitter,
        gauss::GaussianBeamlet{R},
        ray_id::Int) where {R}
    chief_ray = rays(gauss.chief)[ray_id]
    int = ray_id <= length(gauss.chief.segments) ? intersection(gauss.chief.segments[ray_id]) : intersect3d(shape(bs), chief_ray)
    isnothing(int) && return nothing
    ambient = ambient_medium(system)
    trans = resolve_transition(ambient, ambient, chief_ray, normal3d(int))
    return interact3d(system, bs, trans, int, gauss, ray_id)
end

function interact3d(
        system::AbstractSystem,
        bs::ThinBeamsplitter,
        trans::Transition,
        int::AbstractIntersection,
        agb::AstigmaticGaussianBeamlet{R},
        ray_id::Int) where {R}
    c_ray = rays(agb.c)[ray_id]
    int_c = ray_id <= length(agb.c.segments) ? (isnothing(intersection(agb.c.segments[ray_id])) ? int : intersection(agb.c.segments[ray_id])) : int
    
    t_c = _beamsplitter_transmitted_ray(bs, trans, int_c, c_ray)
    r_c = _beamsplitter_reflected_ray(bs, trans, int_c, c_ray)
    
    df = dot(direction(c_ray), normal3d(int_c))
    ϕ = df < 0 ? π : 0
    if ϕ != 0
        r_c = PolarizedRay(position(r_c), direction(r_c), wavelength(r_c), refractive_index(r_c), polarization(r_c) * exp(im * ϕ))
    end
    
    t_aux = Ray{R}[]
    r_aux = Ray{R}[]
    for b in _aux_beams(agb)
        ray_a = rays(b)[ray_id]
        int_a = ray_id <= length(b.segments) ? (isnothing(intersection(b.segments[ray_id])) ? int : intersection(b.segments[ray_id])) : int
        push!(t_aux, _beamsplitter_transmitted_ray(bs, trans, int_a, ray_a))
        push!(r_aux, _beamsplitter_reflected_ray(bs, trans, int_a, ray_a))
    end
    
    t_agb = AstigmaticGaussianBeamlet(
        Beam(t_c),
        Beam(t_aux[1]), Beam(t_aux[2]), Beam(t_aux[3]), Beam(t_aux[4]),
        Beam(t_aux[5]), Beam(t_aux[6]), Beam(t_aux[7]), Beam(t_aux[8])
    )
    r_agb = AstigmaticGaussianBeamlet(
        Beam(r_c),
        Beam(r_aux[1]), Beam(r_aux[2]), Beam(r_aux[3]), Beam(r_aux[4]),
        Beam(r_aux[5]), Beam(r_aux[6]), Beam(r_aux[7]), Beam(r_aux[8])
    )
    
    children!(agb, [t_agb, r_agb])
    return nothing
end

function interact3d(
        system::AbstractSystem,
        bs::ThinBeamsplitter,
        agb::AstigmaticGaussianBeamlet{R},
        ray_id::Int) where {R}
    c_ray = rays(agb.c)[ray_id]
    int = ray_id <= length(agb.c.segments) ? intersection(agb.c.segments[ray_id]) : intersect3d(shape(bs), c_ray)
    isnothing(int) && return nothing
    ambient = ambient_medium(system)
    trans = resolve_transition(ambient, ambient, c_ray, normal3d(int))
    return interact3d(system, bs, trans, int, agb, ray_id)
end

