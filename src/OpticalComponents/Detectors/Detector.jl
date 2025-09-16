#=
TODO
1. implement intensity for ray, pray and gbeamlet
2. implement generic autolims
3. think about solutions for post-solve kinematic changes of Detector
=#

abstract type AbstractDetectorHit end
abstract type AbstractRayHit{T} <: AbstractDetectorHit end
abstract type AbstractBeamletHit{T} <: AbstractDetectorHit end

position(hit::AbstractRayHit) = position(hit.ray) + length(hit.ray) * direction(hit.ray)
direction(hit::AbstractRayHit) = direction(hit.ray)
optical_path_length(hit::AbstractRayHit) = hit.opl

struct RayHit{T} <: AbstractRayHit{T}
    ray::Ray{T}
    opl::T
end

projection_factor(hit::RayHit) = abs(dot(direction(hit), normal3d(intersection(hit.ray))))

struct PolarizedRayHit{T} <: AbstractRayHit{T}
    ray::PolarizedRay{T}
    opl::T
end

polarization(hit::PolarizedRayHit) = polarization(hit.ray)

struct GaussianBeamletHit{T} <: AbstractBeamletHit{T}
    gauss::GaussianBeamlet{T}
    id::Int 
end

mutable struct Detector{T, S <: AbstractShape{T}} <: AbstractDetector{T, S}
    const shape::S
    hits::NullableVector{<:AbstractDetectorHit}
    stop::Bool
end

function Detector(len::Real, cont=false)
    shape = QuadraticFlatMesh(len)
    return Detector(shape, nothing, cont)
end

empty!(d::Detector) = hits!(d, nothing)

hits(d::Detector) = d.hits
hits!(d::Detector, new) = (d.hits = new)

"Allows continued tracing if set to false"
stop(d::Detector) = d.stop

function push!(d::Detector, new::AbstractDetectorHit)
    # set data type on first push
    if isnothing(hits(d))
        hits!(d, [new])
        return nothing
    end
    # if new<:AbstractData does not match, throws push! error 
    push!(hits(d), new)
end

function interact3d(::AbstractSystem, d::Detector, beam::Beam{T, R}, ray::R) where {T <: Real, R <: Ray{T}}
    # Optical path length and wavenumber
    opl = optical_path_length(beam)
    # Push hit data into detector, determine stop
    hit = RayHit(ray, opl)
    push!(d, hit)
    if stop(d)
        # Stop solver (hard target)
        return nothing
    else
        # Continue tracing with hit pos. as starting
        return BeamInteraction{T, R}(
            nothing, 
            Ray{T}(
                position(hit),
                direction(hit),
                nothing,
                wavelength(ray),
                refractive_index(ray)
            )
        )
    end
end

function interact3d(::AbstractSystem, d::Detector, beam::Beam{T, R}, ray::R) where {T <: Real, R <: PolarizedRay{T}}
    # Optical path length and wavenumber
    opl = optical_path_length(beam)
    # Push hit data into detector, determine stop
    hit = PolarizedRayHit(ray, opl)
    push!(d, hit)
    if stop(d)
        # Stop solver (hard target)
        return nothing
    else
        # Continue tracing with hit pos. as starting
        return BeamInteraction{T, R}(
            nothing,
            PolarizedRay{T}(
                position(hit),
                direction(hit),
                nothing,
                wavelength(ray),
                refractive_index(ray),
                polarization(ray)
            )
        )
    end
end

function interact3d(::AbstractSystem, d::Detector, g::GaussianBeamlet{R}, id::Int) where {R}
    push!(d, GaussianBeamletHit(g, id))
    if stop(d)
        # Stop solver (hard target)
        return nothing
    else
        # Continue tracing
        throw(ErrorException("Continued tracing for GaussianBeamlet not yet implemented."))
    end
end