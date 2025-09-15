abstract type AbstractDetectorHit end

struct RayHit{T} <: AbstractDetectorHit
    hit::Point3{T}
    dir::Point3{T}
    proj::T
    opl::T
    k::T
end

struct PolarizedRayHit{T} <: AbstractDetectorHit
    hit::Point3{T}
    dir::Point3{T}
    E0::Point3{Complex{T}}
    opl::T
    k::T
end

struct GaussianBeamletHit{T} <: AbstractDetectorHit
    gauss::GaussianBeamlet{T}
    id::Int 
end

mutable struct Detector{T, S <: AbstractShape{T}} <: AbstractDetector{T, S}
    const shape::S
    hits::NullableVector{<:AbstractDetectorHit}
    stop::Bool
end

function Detector(len::Real, cont=false)
    shape = BMO.QuadraticFlatMesh(len)
    return Detector(shape, nothing, cont)
end 

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
    # Global hit pos
    pos = position(ray) + length(ray) * direction(ray)
    # Global hit dir
    dir = direction(ray)
    # Optical path length and wavenumber
    opl = optical_path_length(beam)
    wvn = 2π/wavelength(ray)
    # Projection factor
    proj = abs(dot(dir, normal3d(intersection(ray))))
    hit = RayHit(pos, dir, proj, opl, wvn)
    # Push hit data into detector, determine stop
    push!(d, hit)
    if stop(d)
        # Stop solver (hard target)
        return nothing
    else
        # Continue tracing
        return BeamInteraction{T, R}(nothing, Ray{T}(pos, dir, nothing, wavelength(ray), refractive_index(ray)))
    end
end

function interact3d(::AbstractSystem, d::Detector, b::Beam{T, R}, r::R) where {T <: Real, R <: PolarizedRay{T}}
    # Global hit pos
    pos = position(ray) + length(ray) * direction(ray)
    # Global hit dir
    dir = direction(ray)
    # Optical path length and wavenumber
    opl = optical_path_length(beam)
    wvn = 2π/wavelength(ray)
    # Push hit data into detector, determine stop
    hit = PolarizedRayHit(pos, dir, polarization(r), opl, wvn)
    push!(d, hit)
    if stop(d)
        # Stop solver (hard target)
        return nothing
    else
        # Continue tracing
        return BeamInteraction{T, R}(nothing, PolarizedRay{T}(pos, dir, nothing, wavelength(ray), refractive_index(ray), polarization(r)))
    end
end

function interact3d(::AbstractSystem, d::Detector, g::GaussianBeamlet{R}, id::Int) where {R}
    push!(d, GaussianBeamletHit(g, id))
    if stop(d)
        # Stop solver (hard target)
        return nothing
    else
        # Continue tracing
        throw(ErrorException("Not yet implemented."))
    end
end