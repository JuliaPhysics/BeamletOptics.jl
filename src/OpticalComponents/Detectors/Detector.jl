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
- [`ray_hits`](@ref)
- [`electric_field`](@ref)
- [`intensity`](@ref)

# Additional information

In general, the detection surface is represented by a flat [`Mesh`](@ref) that has been rotated such that
the surface normals point towards the negative y-axis for the initial positioning. This allows for the definition
of a **left-handed** (x, z) surface coordinate system, where incoming beams intersect against the detector surface normal.

!!! warning "Reset behavior"
    A `Detector` accumulates hits across calls of [`solve_system!`](@ref), on purpose: all beams of a beam group
    and several sources solved one after another superpose on it. To discard the results of previous solves,
    e.g. after changing a component in a loop, call [`initialize!`](@ref)`(system)` for all detectors of a system
    (nested groups included), `empty!(detector)` for a single one, or solve with `initialize = true`.
    Otherwise, the current result is added onto the previous result.

!!! note "Thread safety"
    Beam groups are solved multithreaded. Pushing hits is guarded by the `lock` of the detector, so
    all beams can hit the same detector concurrently. Each task of the solve of a beam group collects its
    hits and stores them at its end, see [`with_hit_buffers`](@ref). The order of the stored hits is not
    deterministic.

!!! warning "Reading hits"
    A hit stores the values of the ray or beamlet at the moment of the solve, in world coordinates.
    Solving again or moving a source or another object does not change the stored hits. Moving the
    detector itself empties it, since its hits would no longer lie on its surface. The hits of
    consecutive solves are kept side by side, and it is up to the user when to evaluate them:
    [`spot_diagram`](@ref) shows all of them, whereas [`electric_field`](@ref), [`intensity`](@ref) and
    [`optical_power`](@ref) superpose all of them coherently, as for beams that exist at the same time.

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
        Vector{RayHit{T}},
        Vector{PolarizedRayHit{T}},
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

- `false` indicates continued tracing
- `true` stops the incoming beams as with any hard target
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

"""
    length(d::Detector) -> Int

The number of hits stored in the detector `d`, `0` before the first solve and after
`empty!(d)`.
"""
length(d::Detector) = isnothing(hits(d)) ? 0 : length(hits(d))

"Allows continued tracing if set to false"
stop(d::Detector) = d.stop

"""
    local_axes(pd::Detector) -> (x̂, ẑ)

The axes of the local (x, z) coordinates of `pd` in global coordinates, with
`x̂ × ẑ = direction(pd)`. `x̂` is flipped since the detector mesh is rotated by `π`.
"""
local_axes(pd::Detector) = (Point3(-orientation(pd)[:, 1]), Point3(orientation(pd)[:, 3]))

"""
    push!(d::Detector, hit::AbstractDetectorHit)

Stores the `hit` in the detector `d`, guarded by the `lock` of `d`, and returns `d`. Within
[`with_hit_buffers`](@ref), e.g. in the solve of a beam group, the `hit` is collected in the buffer
of the current task instead and reaches `d` when the buffer is flushed.
"""
function push!(d::Detector, hit::AbstractDetectorHit)
    # Buffers of the current task, otherwise the detector
    sink = get(task_local_storage(), hit_buffers_key, LockedHits())::Union{LockedHits, HitBuffers}
    push!(sink, d => hit)
    return d
end

"""
    append!(d::Detector, new::Vector{<:AbstractDetectorHit})

Stores all hits of `new` in the detector `d` at once, guarded by the `lock` of `d`. The vector
`new` may become the vector of hits of `d`, i.e. it must not be used afterwards.
"""
function Base.append!(d::Detector, new::Vector{<:AbstractDetectorHit})
    lock(d.lock) do
        # set data type on first push
        if isnothing(hits(d))
            hits!(d, new)
            return nothing
        end
        # if the hit type does not match, throws a conversion error
        append!(hits(d), new)
        return nothing
    end
    return d
end

function interact3d(::AbstractSystem, d::Detector, beam::Beam{T, R},
        ray::R) where {T <: Real, R <: Ray{T}}
    # Optical path length and wavenumber
    opl = optical_path_length(beam)
    # Push hit data into detector, determine stop
    hit = RayHit(ray, opl)
    push!(d, hit)
    if stop(d)
        # Stop solver (hard target)
        return nothing
    else
        # Continue tracing with hit point as starting position
        return BeamInteraction{T, R}(
            nothing,
            Ray{T}(
                hit_point(hit),
                direction(hit),
                nothing,
                wavelength(ray),
                refractive_index(ray)
            )
        )
    end
end

function interact3d(::AbstractSystem, d::Detector, beam::Beam{T, R},
        ray::R) where {T <: Real, R <: PolarizedRay{T}}
    # Optical path length and wavenumber
    opl = optical_path_length(beam)
    # Push hit data into detector, determine stop
    hit = PolarizedRayHit(ray, opl)
    push!(d, hit)
    if stop(d)
        # Stop solver (hard target)
        return nothing
    else
        # Continue tracing with hit point as starting position
        return BeamInteraction{T, R}(
            nothing,
            PolarizedRay{T}(
                hit_point(hit),
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

function interact3d(system::AbstractSystem, d::Detector,
        agb::AstigmaticGaussianBeamlet{R}, id::Int) where {R}
    push!(d, AstigmaticGaussianBeamletHit(agb, id))
    if stop(d)
        # Stop solver (hard target)
        return nothing
    else
        # Continue tracing # FIXME
        throw(ErrorException("Continued tracing for AstigmaticGaussianBeamlet not yet implemented."))
    end
end
