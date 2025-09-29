"""
    AbstractDetectorHit

Abstract supertype for all [`Detector`](@ref) hit records.

A hit of this type encapsulates the interaction of a ray or beamlet with the [`Detector`](@ref) 
surface. Instead of directly accumulating values (e.g. optical power), the hit 
object stores all relevant information about the incident field for a posteriori
evaluation (such as plotting, power integration, or polarization analysis).

Concrete subtypes include:

- [`RayHit`](@ref) — unpolarized geometric ray
- [`PolarizedRayHit`](@ref) — ray carrying a polarization state
- [`GaussianBeamletHit`](@ref) — single Gaussian beamlet
"""
abstract type AbstractDetectorHit end

"""
    AbstractRayHit{T} <: AbstractDetectorHit

Abstract supertype for detector hits produced by [`AbstractRay`](@ref)s. 
Provides a common interface for extracting positional and optical path 
information from rays stored in a detector hit.
Currently the following concrete types are implemented:

- [`RayHit`](@ref)
- [`PolarizedRayHit`](@ref)

# Implementation reqs.

Subtypes of `AbstractRayHit` must implement the following:

## Fields

- `ray`: stores the `AbstractRay` that has intersected the detector
- `opl`: stores the [`optical_path_length`](@ref) of the parent beam (incl. the ray)

## Functions

The interface provides the following functions for the fields above:

- `position`: returns the `ray` position
- `direction`: returns the `ray` direction
- `length`: returns the `ray` length
- `optical_path_length`: returns the `opl`
- `wavenumber`: returns the `ray` wavenumber
- `hit_point`: returns the R³ point of intersection
- `projection_factor`: returns the scalar projection between the surface normal and ray dir.
"""
abstract type AbstractRayHit{T} <: AbstractDetectorHit end

position(hit::AbstractRayHit) = position(hit.ray)
direction(hit::AbstractRayHit) = direction(hit.ray)

length(hit::AbstractRayHit) = length(hit.ray)
optical_path_length(hit::AbstractRayHit) = hit.opl
wavenumber(hit::AbstractRayHit) = wavenumber(hit.ray)

hit_point(hit::AbstractRayHit) = position(hit) + length(hit.ray) * direction(hit)

projection_factor(hit::AbstractRayHit) = abs(dot(direction(hit), normal3d(intersection(hit.ray))))

"Stores a [`Ray`](@ref) hit"
struct RayHit{T} <: AbstractRayHit{T}
    ray::Ray{T}
    opl::T
end

"Stores a [`PolarizedRay`](@ref) hit"
struct PolarizedRayHit{T} <: AbstractRayHit{T}
    ray::PolarizedRay{T}
    opl::T
end

polarization(hit::PolarizedRayHit) = polarization(hit.ray)

"""
    AbstractBeamletHit

Stores beamlet hits. Currently implemented:

- [`GaussianBeamletHit`](@ref)
"""
abstract type AbstractBeamletHit{T} <: AbstractDetectorHit end

"""
    GaussianBeamletHit{T} <: AbstractBeamletHit{T}

Stores a [`GaussianBeamlet`], where `l0` represents the length of the parent beam
up until the current beam section, identified by the `id` index.
"""
struct GaussianBeamletHit{T} <: AbstractBeamletHit{T}
    gauss::GaussianBeamlet{T}
    l0::T
    id::Int 
end

position(hit::GaussianBeamletHit) = position(hit.gauss.chief.rays[hit.id])
direction(hit::GaussianBeamletHit) = direction(hit.gauss.chief.rays[hit.id])

hit_point(hit::GaussianBeamletHit) = position(hit) + length(hit.gauss.chief.rays[hit.id]) * direction(hit)

projection_factor(hit::GaussianBeamletHit) = abs(dot(direction(hit), normal3d(intersection(hit.gauss.chief.rays[hit.id]))))

"""
    Detector <: AbstractDetector

Represents a **flat** rectangular or quadratic, infinitesimally thin surface in R³.
The detector surface is a detection screen that captures incoming ray or beamlet data.
The active surface is discretized in the local R² x-y-coordinate system.
If configured, beams or beamlets can continue tracing after hitting the detector.

# Hits 

Hits are represented via the [`AbstractDetectorHit`](@ref) interface. An empty detector is able to 
detect any kind of incoming hit, but as soon as the initial hit type has been determined, all following
hits must share the same type, i.e. no cross-interaction between hit types is allowed.

# Functions

The following functions allow a posteriori evaluation of hit contributions via e.g. `f(detector)`.
Refer to the respective function documentation.

- [`spot_diagram`](@ref)
- [`electric_field`](@ref)
- [`intensity`](@ref)

# Additional information

In general, the detection surface is represented by a flat [`Mesh`](@ref) that has been rotated such that
the surface normals point towards the negative y-axis for the initial positioning. This allows for the definition
of a **left-handed** (x, z) surface coordinate system, where incoming beams intersect against the detector surface normal.   

!!! warning "Reset behavior"
    The `Detector` must be reset between each call of [`solve_system!`](@ref) in order to
    overwrite previous results using the [`empty!`](@ref) function.
    Otherwise, the current result will be added onto the previous result.

# Fields

- `shape`
  geometry of the active surface, must represent 2D-`field` in `x` any `y` dimensions,
  normal vector direction must adhere to definition above
- `hits`
  a [`NullableVector`](@ref) of hit data, resetable via `empty!`
- `stop`
  a boolean value that allows for continued tracing after "passing through" the detector
- `lock`
  locks the `Detector` for multithreading-safe `push!`ing to the hits vector
"""
mutable struct Detector{T, S <: AbstractShape{T}} <: AbstractDetector{T, S}
    const shape::S
    hits::NullableVector{<:AbstractDetectorHit}
    stop::Bool
    lock::ReentrantLock
end

"""
    Detector(edge_length, stop)

Spawns a quadratic [`Detector`](@ref) surface that is aligned with the neg. y-axis.
The detector edge length can be configured via `edge_length`.
Additionally, continued tracing can be configured via the `stop` flag where

- `true` indicates continued tracing
- `false` stops the incoming beams as with any hard target 
"""
function Detector(edge_length::Real, stop::Bool=true)
    shape = QuadraticFlatMesh(edge_length)
    # rotate surface normal along neg. y-axis
    zrotate3d!(shape, π)
    return Detector(shape, nothing, stop, ReentrantLock())
end

empty!(d::Detector) = hits!(d, nothing)

hits(d::Detector) = d.hits
hits!(d::Detector, new) = (d.hits = new)

"Allows continued tracing if set to false"
stop(d::Detector) = d.stop

function push!(d::Detector, new::AbstractDetectorHit)
    # Lock d to avoid race conditions
    lock(d.lock) do 
        # set data type on first push
        if isnothing(hits(d))
            hits!(d, [new])
            return nothing
        end
        # if new<:AbstractData does not match, throws push! error 
        push!(hits(d), new)
    end
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
    l0 = length(g) - length(g.chief.rays[id])
    push!(d, GaussianBeamletHit(g, l0, id))
    if stop(d)
        # Stop solver (hard target)
        return nothing
    else
        # Continue tracing # FIXME
        throw(ErrorException("Continued tracing for GaussianBeamlet not yet implemented."))
    end
end