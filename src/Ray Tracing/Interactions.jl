#=
Implements photodetector, efield calculation during solve_system!
=#
abstract type AbstractDetector{T, S <: AbstractShape{T}} <: AbstractObject{T, S} end

"""
    Photodetector{T, S <: AbstractShape{T}} <: AbstractDetector{T, S}

Represents a **flat** rectangular or quadratic surface in R³ that is the active surface of a photodetector.
The active surface is discretized in the local R² x-y-coordinate system.
Field contributions Eᵢ are added by the corresponding [`interact3d`](@ref) method.

# Fields

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
    Currently, only the [`GaussianBeamlet`](@ref) is supported.
"""
mutable struct Photodetector{T, S <: AbstractShape{T}} <: AbstractDetector{T, S}
    const shape::S
    x::LinRange{T, Int}
    y::LinRange{T, Int}
    field::Matrix{Complex{T}}
end

function Photodetector(width::T, n::Int) where {T}
    shape = QuadraticFlatMesh(width)
    sz = maximum(vertices(shape))
    x = y = LinRange(-sz, sz, n)
    field = zeros(Complex{T}, n, n)
    return Photodetector{T, typeof(shape)}(shape, x, y, field)
end

function interact3d(::AbstractSystem, ::Photodetector, ::B, ::Ray) where {B <: AbstractBeam}
    @warn "Photodetection for $B not implemented"
    return nothing
end

"""
    interact3d(::AbstractSystem, pd::Photodetector, gauss::GaussianBeamlet, ray_id::Int)

Implements the [`Photodetector`](@ref) interaction with a [`GaussianBeamlet`](@ref).
On hit, the scalar E-field of the `gauss` is added to the current PD field matrix.
Tilt and tip between beam and PD surface are considered via projection factors.
"""
function interact3d(
        ::AbstractSystem, pd::Photodetector, gauss::GaussianBeamlet, ray_id::Int)
    # Select final ray of chief beam
    ray = gauss.chief.rays[ray_id]
    l0 = length(gauss, opl = true) - length(ray, opl = true)
    p0 = position(ray)
    d0 = direction(ray)
    # Preallocate transforms
    T = transpose(orientation(shape(pd)))
    p = position(shape(pd))
    # E-field projection scalar reduction factor
    ray_int = intersection(ray)
    isnothing(ray_int) && return nothing

    proj = abs(dot(d0, normal3d(ray_int)))

    # Add current E-field contribution
    Threads.@threads for j in eachindex(pd.y) # FIXME row column major order?
        y = pd.y[j]
        @inbounds for i in eachindex(pd.x)
            x = pd.x[i]
            # Transform point p on PD into world coordinates
            p1 = Point3(
                T[1, 1] * x + T[1, 3] * y + p[1],
                T[2, 1] * x + T[2, 3] * y + p[2],
                T[3, 1] * x + T[3, 3] * y + p[3]
            )
            # Find projection of p1 onto Gaussian optical axis, i.e. local r and z
            l1 = dot(p1 - p0, d0)
            p2 = p0 + l1 * d0
            r = norm(p1 - p2)
            z = l0 + l1
            # Add field contribution, projection factor accounts for beam spot stretching
            pd.field[i, j] += electric_field(gauss, r, z) * sqrt(proj)
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
reset_photodetector!(pd::Photodetector{T, S}) where {T, S} = (pd.field .= zero(Complex{T}))

"""
    photodetector_resolution!(pd::Photodetector, n::Int)

Sets the resolution of `pd` to `n` × `n`. Note that this resets the current `pd.field`.
"""
function photodetector_resolution!(pd::Photodetector{T, S}, n::Int) where {T, S}
    pd.x = LinRange(pd.x.start, pd.x.stop, n)
    pd.y = LinRange(pd.y.start, pd.y.stop, n)
    pd.field = zeros(Complex{T}, n, n)
    return nothing
end

#=
Implements thin beamsplitter, beam spawning
=#
abstract type AbstractBeamSplitter{T, S <: AbstractShape{T}} <: AbstractObject{T, S} end

"""Models a generic beamsplitter"""
struct BeamSplitter{T <: Real, S <: AbstractShape{T}} <: AbstractBeamSplitter{T, S}
    shape::S
    reflectance::T
    transmittance::T
end

reflectance(bs::BeamSplitter) = bs.reflectance
transmittance(bs::BeamSplitter) = bs.transmittance

"""
    ThinBeamSplitter(width::T, reflectance::Real=0.5) where {T}

Creates a zero-thickness, lossless, non-polarizing quadratic rectangle beamsplitter where

- `width`: is the edge length
- `reflectance`: determines how much light is **reflected**, i.e. 0.7 for a 70:30 splitter

# Additional information

!!! info "Reflectance"
    The input value for the `reflectance` R is normed such that R² + T² = 1, where T is the `transmittance`.
    The transmittance is calculated via T = √(1 - R²).

!!! warning "Reflection phase jump"
    Note that the reflection phase jump θᵣ is implemented by the individual [`interact3d`](@ref)-methods. Refer to them for more information.
"""
function ThinBeamSplitter(width::T, reflectance::Real = 0.5) where {T}
    if reflectance ≥ 1 || reflectance ≈ 0
        error("Splitting ratio ∈ (0, 1)!")
    end
    shape = QuadraticFlatMesh(width)
    Reflected = sqrt(reflectance)
    Transmitted = sqrt(1 - Reflected^2)
    return BeamSplitter(shape, Reflected, Transmitted)
end

Base.isvalid(bs::AbstractBeamSplitter) = reflectance(bs)^2 + transmittance(bs)^2 ≈ 1

@inline function _beamsplitter_transmitted_beam(
        ::AbstractBeamSplitter, ::Beam{T, R}, ray::R) where {T <: Real, R <: Ray{T}}
    pos = position(ray) + length(ray) * direction(ray)
    dir = direction(ray)
    return Beam(Ray(pos, dir, wavelength(ray)))
end

@inline function _beamsplitter_reflected_beam(
        ::AbstractBeamSplitter, ::Beam{T, R}, ray::R) where {T <: Real, R <: Ray{T}}
    normal = normal3d(intersection(ray))
    pos = position(ray) + length(ray) * direction(ray)
    dir = reflection3d(direction(ray), normal)
    return Beam(Ray(pos, dir, wavelength(ray)))
end

@inline function _beamsplitter_transmitted_beam(bs::AbstractBeamSplitter, ::Beam{T, R},
        ray::R) where {T <: Real, R <: PolarizedRay{T}}
    J = @SArray [transmittance(bs) 0 0; 0 transmittance(bs) 0; 0 0 1]
    pos = position(ray) + length(ray) * direction(ray)
    dir = direction(ray)
    E0 = _calculate_global_E0(dir, dir, J, polarization(ray))
    return Beam(PolarizedRay(pos, dir, wavelength(ray), E0))
end

@inline function _beamsplitter_reflected_beam(bs::AbstractBeamSplitter, ::Beam{T, R},
    ray::R) where {T <: Real, R <: PolarizedRay{T}}
    J = @SArray [-reflectance(bs) 0 0; 0 reflectance(bs) 0; 0 0 1]
    normal = normal3d(intersection(ray))
    pos = position(ray) + length(ray) * direction(ray)
    in_dir = direction(ray)
    out_dir = reflection3d(in_dir, normal)
    E0 = _calculate_global_E0(in_dir, out_dir, J, polarization(ray))
    return Beam(PolarizedRay(pos, out_dir, wavelength(ray), E0))
end

function interact3d(::AbstractSystem, bs::BeamSplitter, beam::Beam{T, R},
        ray::R) where {T <: Real, R <: AbstractRay{T}}
    # Push transmitted and reflected beams to system
    children!(beam,
        [_beamsplitter_transmitted_beam(bs, beam, ray), _beamsplitter_reflected_beam(bs, beam, ray)])
    # Stop for beam spawning
    return nothing
end

@inline function _beamsplitter_transmitted_beam(bs::AbstractBeamSplitter, gauss::GaussianBeamlet, ray_id::Int)
    # Transmitted Gaussian (no phase flip)
    chief = _beamsplitter_transmitted_beam(bs, gauss.chief, rays(gauss.chief)[ray_id])
    waist = _beamsplitter_transmitted_beam(bs, gauss.waist, rays(gauss.waist)[ray_id])
    divergence = _beamsplitter_transmitted_beam(bs, gauss.divergence, rays(gauss.divergence)[ray_id])
    λ = wavelength(gauss)
    w0 = gauss_parameters(gauss, length(gauss))[4]
    E0 = transmittance(bs) * electric_field(gauss) * (beam_waist(gauss) / w0)
    return GaussianBeamlet(chief, waist, divergence, λ, w0, E0)
end

@inline function _beamsplitter_reflected_beam(bs::AbstractBeamSplitter, gauss::GaussianBeamlet, ray_id::Int)
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
    interact3d(::AbstractSystem, bs::BeamSplitter, gauss::GaussianBeamlet, ray_id::Int)

Models the interaction between a [`BeamSplitter`](@ref) and a [`GaussianBeamlet`](@ref).
For more information refer to [`ThinBeamSplitter`](@ref).

# Reflection phase jump

The reflection phase jump is modeled here as θᵣ = π for simplicity. This is since in practice it will have only a relative effect on the signal at the detector for interferometric setups.
The phase jump is applied to the reflected portion of any incoming beam that faces the `BeamSplitter` normal vector, which assumes that the splitter has an unambigous normal, i.e. a 2D mesh.
This is intended to model the effect of the Fresnel equations without full polarization calculus.
"""
function interact3d(::AbstractSystem, bs::BeamSplitter, gauss::GaussianBeamlet, ray_id::Int)
    # Phase flip
    ray = gauss.chief.rays[end]
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
