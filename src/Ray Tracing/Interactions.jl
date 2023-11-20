
"""
    AbstractReflectiveOptic <: AbstractObject

A generic type to represent `AbstractObject`s which reflect incoming rays. The main function of `interact3d` should be akin to [`reflection3d`](@ref).
"""
abstract type AbstractReflectiveOptic <: AbstractObject end

function interact3d(::AbstractSystem,
        ::AbstractReflectiveOptic,
        ::Beam{R},
        ray::Ray{R}) where {R <: Real}
    normal = intersection(ray).n
    npos = position(ray) + length(ray) * direction(ray)
    ndir = reflection3d(direction(ray), normal)
    return BeamInteraction{R}(nothing, Ray{R}(uuid4(), npos, ndir, nothing, ray.parameters))
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
Subtypes of `AbstractRefractiveOptic` should implement all supertype reqs as well as:

# Fields
- `n::Function`: a function which returns the refractive index for a wavelength λ
"""
abstract type AbstractRefractiveOptic <: AbstractObject end

refractive_index(object::AbstractRefractiveOptic) = object.n

function interact3d(::AbstractSystem,
        object::AbstractRefractiveOptic,
        ::Beam{R},
        ray::Ray{R}) where {R}
    # Check dir. of ray and surface normal
    normal = intersection(ray).n
    λ = wavelength(ray)
    if dot(direction(ray), normal) < 0
        # "Outside prism"
        n1 = 1.0
        n2 = refractive_index(object)(λ)
    else
        # "Inside prism"
        n1 = refractive_index(object)(λ)
        n2 = 1.0
        normal *= -1
    end
    # Calculate new dir. and pos.
    ndir = refraction3d(direction(ray), normal, n1, n2)
    npos = position(ray) + length(ray) * direction(ray)
    # Hint is the current object ID
    return BeamInteraction{R}(id(object),
        Ray{R}(uuid4(), npos, ndir, nothing, Parameters(λ, n2)))
end

struct Lens{S <: AbstractShape, T <: Function} <: AbstractRefractiveOptic
    id::UUID
    shape::S
    n::T
end

SphericalLens(r1, r2, l, d, n) = SphericalLens(r1, r2, l, d, λ -> n)

"""
    SphericalLens(r1, r2, l, d=1inch, n=λ->1.5)

Creates a spherical lens based on:

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

If `l` is set to zero, an ideal [`ThinLens`](@ref) will be created. However, note that the actual lens thickness will be different from zero.
"""
function SphericalLens(r1::Real, r2::Real, l::Real, d::Real = 1inch, n::Function = λ -> 1.5)
    # # Test for thin lens
    if iszero(l)
        shape = ThinLensSDF(r1, r2, d)
        # goto to avoid overwrite
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

ThinLens(R1::Real, R2::Real, d::Real, n::Real) = ThinLens(R1, R2, d, x -> n)
function ThinLens(R1::Real, R2::Real, d::Real, n::Function)
    shape = ThinLensSDF(R1, R2, d)
    return Lens(uuid4(), shape, n)
end

struct Prism{S <: AbstractShape, T <: Function} <: AbstractRefractiveOptic
    id::UUID
    shape::S
    n::T
end

#=
Implements photodetector, efield calculation during solve_system!
=#
abstract type AbstractDetector <: AbstractObject end

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
    l = length(gauss)
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

@inline function _beamsplitter_transmitted_beam(ray::Ray)
    pos = position(ray) + length(ray) * direction(ray)
    dir = direction(ray)
    return Beam(Ray(pos, dir, wavelength(ray)))
end

@inline function _beamsplitter_reflected_beam(ray::Ray)
    normal = normal3d(intersection(ray))
    pos = position(ray) + length(ray) * direction(ray)
    dir = reflection3d(direction(ray), normal)
    return Beam(Ray(pos, dir, wavelength(ray)))
end

function interact3d(::AbstractSystem, ::BeamSplitter, beam::Beam{R}, ray::Ray{R}) where {R}
    # Push transmitted and reflected beams to system
    children!(beam,
        [_beamsplitter_transmitted_beam(ray), _beamsplitter_reflected_beam(ray)])
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
    chief = _beamsplitter_transmitted_beam(rays(gauss.chief)[ray_id])
    waist = _beamsplitter_transmitted_beam(rays(gauss.waist)[ray_id])
    divergence = _beamsplitter_transmitted_beam(rays(gauss.divergence)[ray_id])
    λ = wavelength(gauss)
    w0 = gauss_parameters(gauss, length(gauss))[4]
    E0 = transmittance(bs) * beam_amplitude(gauss)
    t = GaussianBeamlet(chief, waist, divergence, λ, w0, E0)
    # Reflected gauss
    chief = _beamsplitter_reflected_beam(rays(gauss.chief)[ray_id])
    waist = _beamsplitter_reflected_beam(rays(gauss.waist)[ray_id])
    divergence = _beamsplitter_reflected_beam(rays(gauss.divergence)[ray_id])
    λ = wavelength(gauss)
    w0 = gauss_parameters(gauss, length(gauss))[4]
    E0 = reflectance(bs) * beam_amplitude(gauss)
    r = GaussianBeamlet(chief, waist, divergence, λ, w0, E0)
    children!(gauss, [t, r])
    return nothing
end
