"""
    AbstractDetector <: AbstractObject

A generic representation of a detector that evaluates [`AbstractBeam`](@ref) data during and/or after interaction.
Refer to e.g. [`Detector`](@ref) for more information.

# Implementation reqs.

Subtypes of `AbstractDetector` should implement all supertype requirements as well as:

## Functions

- `interact3d`: see e.g. [`Detector`](@ref) for reference
- `empty!`: resets data stored in the detector, see below

# Additional information

The information provided below applies to the standard functional implementation of this type and may be overwritten
by specialized subtypes.

## Data mutability

In order to model field superposition effects, the concrete implementation of an `AbstractDetector` should be a **mutable struct**
with a `const`ant `shape` field. This is necessary, since (sub-)beams will interact sequentially with the detector during [`solve_system!`](@ref).
Only if the data can be accumulated sequentially, multiple beam interactions can be captured for a complex system, e.g. an interferometer.

## Data reset

Since e.g. E-field data is supposed to be accumulated by mutability of the detector data, the burden of resetting the data for a new solver call
is placed on the user. This function should be called `empty!`.
"""
abstract type AbstractDetector{T, S <: AbstractShape{T}} <: AbstractObject{T, S} end

"""
    empty!(detector)

Resets the field data of the `detector`. Must be implemented for each concrete subtype of [`AbstractDetector`](@ref).
"""
function empty!(::D) where D <: AbstractDetector
    @warn "Detector reset logic for $D not implemented"
    return nothing
end

"""
    AbstractDetectorHit

Abstract supertype for all [`Detector`](@ref) hit records.

A hit of this type encapsulates the interaction of a ray or beamlet with the [`Detector`](@ref) 
surface. Instead of directly accumulating values (e.g. optical power), the hit 
object stores all relevant information about the incident field for a posteriori
evaluation (such as plotting, power integration, or polarization analysis).

!!! info
    Note that new subtypes of `AbstractDetectorHit` must be manually added to the `Detector`
    definition to be eligible.

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

This type must provide all fields required to satisfy the getter functions listed below.

## Functions

The interface provides the following functions for the fields above:

- `position`: returns the `ray` position
- `direction`: returns the `ray` direction
- `normal3d`: returns the surface normal at the point of intersection

- `optical_path_length`: returns the `opl`
- `wavelength`: returns the `ray` wavelength
- `wavenumber`: returns the `ray` wavenumber

- `hit_point`: returns the R³ point of intersection

- `projection_factor`: returns the scalar projection between the surface normal and ray dir.
"""
abstract type AbstractRayHit{T} <: AbstractDetectorHit end

position(hit::AbstractRayHit) = hit.pos
direction(hit::AbstractRayHit) = hit.dir
normal3d(hit::AbstractRayHit) = hit.nml

optical_path_length(hit::AbstractRayHit) = hit.opl
wavelength(hit::AbstractRayHit) = hit.λ
wavenumber(hit::AbstractRayHit) = wavenumber(wavelength(hit))

hit_point(h::AbstractRayHit) = h.hit

projection_factor(hit::AbstractRayHit) = abs(dot(direction(hit), normal3d(hit)))

"Stores a [`Ray`](@ref) hit"
struct RayHit{T} <: AbstractRayHit{T}
    pos::Point3{T}
    dir::Point3{T}
    nml::Point3{T}
    hit::Point3{T}
    opl::T
    λ::T
    function RayHit(ray::Ray{T}, opl) where T
        # ray data
        pos = position(ray)
        dir = direction(ray)
        nml = normal3d(intersection(ray))
        hit = pos + length(ray) * dir
        return new{T}(pos, dir, nml, hit, opl, wavelength(ray))
    end
end

"Stores a [`PolarizedRay`](@ref) hit"
struct PolarizedRayHit{T} <: AbstractRayHit{T}
    pos::Point3{T}
    dir::Point3{T}
    nml::Point3{T}
    hit::Point3{T}
    opl::T
    λ::T
    E0::Point3{Complex{T}}
    function PolarizedRayHit(ray::PolarizedRay{T}, opl) where T
        # ray data
        pos = position(ray)
        dir = direction(ray)
        nml = normal3d(intersection(ray))
        hit = pos + length(ray) * dir
        return new{T}(pos, dir, nml, hit, opl, wavelength(ray), polarization(ray))
    end
end

polarization(hit::PolarizedRayHit) = hit.E0

"""
    AbstractBeamletHit

Stores beamlet hits. Currently implemented:

- [`GaussianBeamletHit`](@ref)
"""
abstract type AbstractBeamletHit{T} <: AbstractDetectorHit end

include("Hits/GaussianBeamletHit.jl")

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

!!! warning "Moving after solving"
    Do not move the detector before calculating all parameters of interest for the current system configuration.
    Since the detector stores pointers to the current system and beam states, silent errors might occur.

# Fields

- `shape`:
  geometry of the active surface, must represent 2D-`field` in `x` any `y` dimensions,
  normal vector direction must adhere to definition above
- `hits`:
  a union of `Nothing` and all implemented [`AbstractDetectorHit`](@ref)s, resettable via `empty!`
  (note that only one type is allowed at any time)
- `stop`:
  a boolean value that allows for continued tracing after "passing through" the detector
- `lock`:
  locks the `Detector` for multithreading-safe `push!`ing to the hits vector
"""
mutable struct Detector{T, S <: AbstractShape{T}} <: AbstractDetector{T, S}
    const shape::S
    # direct reference to avoid UnionAny from AbstractDetectorHit
    hits::Union{
        Nothing,
        Vector{RayHit{T}},
        Vector{PolarizedRayHit{T}},
        Vector{GaussianBeamletHit{T}}
    }
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
    push!(d, GaussianBeamletHit(g, id))
    if stop(d)
        # Stop solver (hard target)
        return nothing
    else
        # Continue tracing # FIXME
        throw(ErrorException("Continued tracing for GaussianBeamlet not yet implemented."))
    end
end

# include eval and helper functions
include("DetectorUtils.jl")
include("SpotDiagram.jl")
include("FieldCalculation.jl")