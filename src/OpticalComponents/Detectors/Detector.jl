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
abstract type AbstractDetector{T} <: AbstractObject{T} end

"""
    empty!(detector)

Resets the field data of the `detector`. Must be implemented for each concrete subtype of [`AbstractDetector`](@ref).
"""
function empty!(::D) where {D <: AbstractDetector}
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

length(hit::AbstractRayHit) = length(hit.intersection)
optical_path_length(hit::AbstractRayHit) = hit.opl
wavenumber(hit::AbstractRayHit) = wavenumber(hit.ray)

hit_point(hit::AbstractRayHit) = position(hit.intersection)

function projection_factor(hit::AbstractRayHit)
    abs(dot(direction(hit), normal3d(hit.intersection)))
end

"Stores a [`Ray`](@ref) hit"
struct RayHit{T, R <: Ray{T}} <: AbstractRayHit{T}
    ray::R
    intersection::Intersection{T}
    opl::T
end

RayHit(ray::Ray{T}, int::Intersection{T}, opl::Real) where {T} = RayHit{T, typeof(ray)}(ray, int, T(opl))
RayHit(ray::Ray{T}, opl::Real) where {T} = RayHit{T, typeof(ray)}(ray, Intersection(zero(T), position(ray), Point3{T}(0,0,1)), T(opl))
RayHit(ray::Ray{T}) where {T} = RayHit{T, typeof(ray)}(ray, isnothing(intersection(ray)) ? Intersection(zero(T), position(ray), Point3{T}(0,0,1)) : intersection(ray), accumulated_opl(ray))

"Stores a [`PolarizedRay`](@ref) hit"
struct PolarizedRayHit{T, R <: PolarizedRay{T}} <: AbstractRayHit{T}
    ray::R
    intersection::Intersection{T}
    opl::T
end

PolarizedRayHit(ray::PolarizedRay{T}, int::Intersection{T}, opl::Real) where {T} = PolarizedRayHit{T, typeof(ray)}(ray, int, T(opl))
PolarizedRayHit(ray::PolarizedRay{T}, opl::Real) where {T} = PolarizedRayHit{T, typeof(ray)}(ray, Intersection(zero(T), position(ray), Point3{T}(0,0,1)), T(opl))
PolarizedRayHit(ray::PolarizedRay{T}) where {T} = PolarizedRayHit{T, typeof(ray)}(ray, isnothing(intersection(ray)) ? Intersection(zero(T), position(ray), Point3{T}(0,0,1)) : intersection(ray), accumulated_opl(ray))

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
    # Cache
    p0::Point3{T}
    d0::Point3{T}
    sqrt_proj::T
    w_max::T
end

"""
    AstigmaticGaussianBeamletHit{T} <: AbstractBeamletHit{T}

Stores an [`AstigmaticGaussianBeamlet`], where `l0` represents the length of the parent beam
up until the current beam section, identified by the `id` index.
"""
struct AstigmaticGaussianBeamletHit{T} <: AbstractBeamletHit{T}
    agb::AstigmaticGaussianBeamlet{T}
    l0::T
    id::Int
    # Cache
    p0::Point3{T}
    d0::Point3{T}
    h1::Point3{Complex{T}}
    u1::Point3{Complex{T}}
    h2::Point3{Complex{T}}
    u2::Point3{Complex{T}}
    area_ref::Complex{T}
    k0::T
    Δl::T
    n_eff::T
    E_ref_amp::Complex{T}
    sqrt_proj::T
    w_max::T
end

position(hit::GaussianBeamletHit) = position(rays(hit.gauss.chief)[hit.id])
direction(hit::GaussianBeamletHit) = direction(rays(hit.gauss.chief)[hit.id])

position(hit::AstigmaticGaussianBeamletHit) = position(rays(hit.agb.c)[hit.id])
direction(hit::AstigmaticGaussianBeamletHit) = direction(rays(hit.agb.c)[hit.id])

function hit_point(hit::GaussianBeamletHit)
    seg = hit.gauss.chief.segments[hit.id]
    return position(seg) + length(seg) * direction(seg)
end
function hit_point(hit::AstigmaticGaussianBeamletHit)
    seg = hit.agb.c.segments[hit.id]
    return position(seg) + length(seg) * direction(seg)
end

function projection_factor(hit::GaussianBeamletHit)
    seg = hit.gauss.chief.segments[hit.id]
    abs(dot(direction(hit), normal3d(seg)))
end
function projection_factor(hit::AstigmaticGaussianBeamletHit)
    seg = hit.agb.c.segments[hit.id]
    abs(dot(direction(hit), normal3d(seg)))
end


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
mutable struct Detector{T, S <: AbstractShape{T}} <: AbstractDetector{T}
    const shape::S
    # direct reference to avoid UnionAny from AbstractDetectorHit
    hits::Union{
        Nothing,
        Vector{<:RayHit{T}},
        Vector{<:PolarizedRayHit{T}},
        Vector{GaussianBeamletHit{T}},
        Vector{AstigmaticGaussianBeamletHit{T}}
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
function Detector(edge_length::Real, stop::Bool = true)
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

function interact3d(::AbstractSystem, d::Detector, beam::Beam{T},
        ray::Ray{T}) where {T <: Real}
    int = isnothing(intersection(ray)) ? intersect3d(shape(d), ray) : intersection(ray)
    opl = optical_path_length(beam)
    hit = RayHit(ray, int, opl)
    push!(d, hit)
    if stop(d)
        # Stop solver (hard target)
        return nothing
    else
        # Continue tracing with hit pos. as starting
        return BeamInteraction{T, typeof(ray)}(
            nothing,
            Ray(
                position(hit),
                direction(hit),
                wavelength(ray),
                refractive_index(ray)
            )
        )
    end
end

function interact3d(::AbstractSystem, d::Detector, beam::Beam{T},
        ray::PolarizedRay{T}) where {T <: Real}
    int = isnothing(intersection(ray)) ? intersect3d(shape(d), ray) : intersection(ray)
    opl = optical_path_length(beam)
    hit = PolarizedRayHit(ray, int, opl)
    push!(d, hit)
    if stop(d)
        # Stop solver (hard target)
        return nothing
    else
        # Continue tracing with hit pos. as starting
        return BeamInteraction{T, typeof(ray)}(
            nothing,
            PolarizedRay(
                position(hit),
                direction(hit),
                wavelength(ray),
                refractive_index(ray),
                polarization(ray)
            )
        )
    end
end

function interact3d(::AbstractSystem, d::Detector, g::GaussianBeamlet{R}, id::Int) where {R}
    # rays(g.chief) and g.chief.segments are the same vector (a beam is never empty),
    # so `id` is either a valid index into both or genuinely out of range.
    seg = g.chief.segments[id]
    l0 = length(g) - length(seg)
    p0 = position(seg)
    d0 = direction(seg)
    nml = isnothing(intersection(seg)) ? Point3{R}(0, 0, 1) : normal3d(intersection(seg))
    sqrt_proj = sqrt(abs(dot(d0, nml)))

    # Calculate actual beam radius at detector for auto-limits
    w_at_detector, _, _, _ = gauss_parameters(g, length(g))
    w_max = w_at_detector

    push!(d, GaussianBeamletHit(g, l0, id, p0, d0, sqrt_proj, w_max))
    if stop(d)
        # Stop solver (hard target)
        return nothing
    else
        # Continue tracing # FIXME
        throw(ErrorException("Continued tracing for GaussianBeamlet not yet implemented."))
    end
end

function interact3d(system::AbstractSystem, d::Detector,
        agb::AstigmaticGaussianBeamlet{R}, id::Int) where {R}
    # rays(agb.c) and agb.c.segments are the same vector (a beam is never empty),
    # so `id` is either a valid index into both or genuinely out of range.
    chief = agb.c.segments[id]
    l0 = length(agb) - length(chief)

    p0 = position(chief)
    d0 = direction(chief)
    k0 = 2π / wavelength(chief)
    nml = isnothing(intersection(chief)) ? Point3{R}(0, 0, 1) : normal3d(intersection(chief))
    sqrt_proj = sqrt(abs(dot(d0, nml)))

    # Parabasal parameters at segment start (p0)
    h1, u1, h2, u2, _ = parabasal_ray_parameters(agb, p0, id)

    # Reference normalization at z=0
    p0n, in_ = point_on_beam(agb, 0.0)
    dirn = direction(rays(agb.c)[in_])
    h1n, _, h2n, _, _ = parabasal_ray_parameters(agb, p0n, in_)
    area_ref = _pseudo_cross2d(h1n, h2n, dirn)

    # Extract complex reference amplitude
    E_vec = polarization(rays(agb.c)[in_])
    max_idx = argmax(abs.(E_vec))
    E_ref_amp = Complex{R}(norm(E_vec) * cis(angle(E_vec[max_idx])))

    # OPL correction (Δl)
    p_parent = agb.parent
    l_parent = isnothing(p_parent) ? 0.0 : length(p_parent)
    opl_parent = isnothing(p_parent) ? 0.0 : optical_path_length(p_parent)

    Δl = opl_parent - l_parent
    z_sum = l_parent
    for j in 1:(id - 1)
        seg_j = agb.c.segments[j]
        Δl += optical_path_length(seg_j) - length(seg_j)
        z_sum += length(seg_j)
    end
    # Note: the (n-1)*z term in parabasal_field depends on (z - z_sum).
    # Propagate semi-axes to detector for w_max (auto-limits)
    H1 = h1 + length(chief) * u1
    H2 = h2 + length(chief) * u2
    w_max = max(norm(H1), norm(H2))
    n_eff = refractive_index(agb, id)

    push!(d,
        AstigmaticGaussianBeamletHit(
            agb, l0, id, p0, d0, h1, u1, h2, u2, area_ref, k0, Δl,
            n_eff, E_ref_amp, sqrt_proj, w_max
        ))
    if stop(d)
        # Stop solver (hard target)
        return nothing
    else
        # Continue tracing # FIXME
        throw(ErrorException("Continued tracing for AstigmaticGaussianBeamlet not yet implemented."))
    end
end


# include eval and helper functions
include("DetectorUtils.jl")
include("SpotDiagram.jl")
include("FieldCalculation.jl")
