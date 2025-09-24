#=
TODO
1. implement intensity for 
    1. ray ✓
    2. pray
    3. gbeamlet ✓
2. implement autolims
    1. rays ✓
    2. gbeamlet ✓
3. fix GB hit getter syntax
4. think about solutions for post-solve kinematic changes of Detector
5. add docs to fcts
    1. calc_local_pos ✓
    2. calc_local_lims ✓
    3. electric_field ✓
    4. intensity ✓
    5. ellipse ✓
    6. Detector 
6. replace old detectors and fix docs/testcases
    1. Photodetector
    2. Spotdetector
    3. PSFDetector
=#

abstract type AbstractDetectorHit end
abstract type AbstractRayHit{T} <: AbstractDetectorHit end
abstract type AbstractBeamletHit{T} <: AbstractDetectorHit end

position(hit::AbstractRayHit) = position(hit.ray)
direction(hit::AbstractRayHit) = direction(hit.ray)

length(hit::AbstractRayHit) = length(hit.ray)
optical_path_length(hit::AbstractRayHit) = hit.opl
wavenumber(hit::AbstractRayHit) = wavenumber(hit.ray)

hit_point(hit::AbstractRayHit) = position(hit) + length(hit.ray) * direction(hit)

projection_factor(hit::AbstractRayHit) = abs(dot(direction(hit), normal3d(intersection(hit.ray))))

struct RayHit{T} <: AbstractRayHit{T}
    ray::Ray{T}
    opl::T
end

struct PolarizedRayHit{T} <: AbstractRayHit{T}
    ray::PolarizedRay{T}
    opl::T
end

polarization(hit::PolarizedRayHit) = polarization(hit.ray)

struct GaussianBeamletHit{T} <: AbstractBeamletHit{T}
    gauss::GaussianBeamlet{T}
    id::Int 
end

position(hit::GaussianBeamletHit) = position(hit.gauss.chief.rays[hit.id])
direction(hit::GaussianBeamletHit) = direction(hit.gauss.chief.rays[hit.id])

length(hit::GaussianBeamletHit) = length(hit.gauss) - length(hit.gauss.chief.rays[hit.id])

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
    The `Photodetector` must be reset between each call of [`solve_system!`](@ref) in order to
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

"""
mutable struct Detector{T, S <: AbstractShape{T}} <: AbstractDetector{T, S}
    const shape::S
    hits::NullableVector{<:AbstractDetectorHit}
    stop::Bool
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
    return Detector(shape, nothing, stop)
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
        # Continue tracing # FIXME
        throw(ErrorException("Continued tracing for GaussianBeamlet not yet implemented."))
    end
end