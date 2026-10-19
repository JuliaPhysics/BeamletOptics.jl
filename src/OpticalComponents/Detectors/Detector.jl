"""
    AbstractDetector <: AbstractObject

A generic representation of a detector that evaluates [`AbstractBeam`](@ref) data during and/or after interaction.
Refer to e.g. [`Detector`](@ref) for more information.

# Implementation reqs.

Subtypes of `AbstractDetector` should implement all supertype requirements as well as:

## Functions

- `interact3d`: see e.g. [`Detector`](@ref) for reference
- `empty!`: resets data stored in the detector, see below. [`initialize!`](@ref) calls it for all detectors of a system

# Additional information

The information provided below applies to the standard functional implementation of this type and may be overwritten
by specialized subtypes.

## Data mutability

In order to model field superposition effects, the concrete implementation of an `AbstractDetector` should be a **mutable struct**
with a `const`ant `shape` field. This is necessary, since (sub-)beams will interact sequentially with the detector during [`solve_system!`](@ref).
Only if the data can be accumulated sequentially, multiple beam interactions can be captured for a complex system, e.g. an interferometer.

## Data reset

Since e.g. E-field data is supposed to be accumulated by mutability of the detector data, the burden of resetting the data for a new solver call
is placed on the user. This function should be called `empty!`. A single detector is reset with `empty!(detector)`, all
detectors of a [`System`](@ref) or [`ObjectGroup`](@ref) with [`initialize!`](@ref), which calls `empty!` on every `AbstractDetector`.
The solver does not reset detectors unless asked to (keyword `initialize` of [`solve_system!`](@ref)), since the beams of a beam group
and several sources solved one after another are meant to superpose.

Every kinematic function empties a detector, since the stored data belong to the pose of the detector during the solve, see
[`translate3d!`](@ref)`(::Movable, ::AbstractDetector, offset)`. A subtype gets this from its `empty!` method without further code.
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

initialize!(d::AbstractDetector) = (empty!(d); nothing)

"""
    translate3d!(::Movable, detector::AbstractDetector, offset)

Empties the `detector` (see `empty!`) and moves it by `offset`. Its hits lie on the surface
where it was during the solve, hence they are discarded when it moves. All kinematic functions
are derived from this one and [`rotate3d!`](@ref)`(::Movable, ::AbstractDetector, R)`, so each of
them empties the detector, also when it moves as a member of an [`ObjectGroup`](@ref).
"""
function translate3d!(::Movable, d::AbstractDetector, offset)
    empty!(d)
    translate3d!(shape_trait_of(d), d, offset)
    return nothing
end

"""
    rotate3d!(::Movable, detector::AbstractDetector, R::AbstractMatrix)

Empties the `detector` (see `empty!`) and rotates it by `R` about its `position`, see
[`translate3d!`](@ref)`(::Movable, ::AbstractDetector, offset)`.
"""
function rotate3d!(::Movable, d::AbstractDetector, R::AbstractMatrix)
    empty!(d)
    rotate3d!(shape_trait_of(d), d, R)
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
Currently the following concrete types are implemented:

- [`RayHit`](@ref)
- [`PolarizedRayHit`](@ref)

A hit stores the values of the ray at the moment it hit the detector, not the ray itself: a
new solve of the beam or a move of any object does not change it.

# Implementation reqs.

Subtypes of `AbstractRayHit` must implement the following functions:

- `hit_point`: the R³ point of intersection
- `direction`: the direction of the ray
- `optical_path_length`: the optical path length of the beam up to the hit point (incl. its parents)
- `wavelength`: the vacuum wavelength of the ray
- `wavenumber`: the vacuum wavenumber `2π / λ`
- `refractive_index`: the refractive index in front of the detector
- `projection_factor`: `|d · n|` of the ray direction `d` and the surface normal `n`; `1` if the
  ray ended without an intersection
"""
abstract type AbstractRayHit{T} <: AbstractDetectorHit end

hit_point(hit::AbstractRayHit) = hit.p
direction(hit::AbstractRayHit) = hit.d
optical_path_length(hit::AbstractRayHit) = hit.opl
wavelength(hit::AbstractRayHit) = hit.λ
wavenumber(hit::AbstractRayHit) = wavenumber(wavelength(hit))
refractive_index(hit::AbstractRayHit) = hit.n
projection_factor(hit::AbstractRayHit) = hit.proj

"""
    RayHit(ray::Ray, opl)

Stores the hit of a [`Ray`](@ref) as values: hit point `p`, direction `d`, optical path length `opl`
of the beam up to the hit point, wavelength `λ`, refractive index `n` and projection factor `proj`.
"""
struct RayHit{T} <: AbstractRayHit{T}
    p::Point3{T}
    d::Point3{T}
    opl::T
    λ::T
    n::T
    proj::T
end

RayHit(ray::Ray{T}, opl::Real) where {T} =
    RayHit{T}(hit_point(ray), direction(ray), opl, wavelength(ray), refractive_index(ray), projection_factor(ray))

"""
    PolarizedRayHit(ray::PolarizedRay, opl)

Stores the hit of a [`PolarizedRay`](@ref) as values: those of a [`RayHit`](@ref) and the field
vector `E` of the ray in global coordinates.
"""
struct PolarizedRayHit{T} <: AbstractRayHit{T}
    p::Point3{T}
    d::Point3{T}
    opl::T
    λ::T
    n::T
    proj::T
    E::Point3{Complex{T}}
end

PolarizedRayHit(ray::PolarizedRay{T}, opl::Real) where {T} =
    PolarizedRayHit{T}(hit_point(ray), direction(ray), opl, wavelength(ray), refractive_index(ray),
        projection_factor(ray), polarization(ray))

polarization(hit::PolarizedRayHit) = hit.E

"""
    AbstractBeamletHit

Stores beamlet hits. Currently implemented:

- [`GaussianBeamletHit`](@ref)
"""
abstract type AbstractBeamletHit{T} <: AbstractDetectorHit end

"""
    GaussianBeamletHit{T} <: AbstractBeamletHit{T}

Stores the segment of a [`GaussianBeamlet`](@ref) that hit the detector as values, without a
reference to the beamlet: `p0` and `d0` are start and direction of the segment, `l_seg` its
length, `l0` the length of the beam (including its parents) up to `p0`, `Δl0` the optical minus
the geometric path and `ψ0` the Gouy phase that the beam has collected up to there. `segment`
holds the waist and divergence rays of the segment (see [`BeamletOptics.GaussianSegment`](@ref)),
from which the field follows as from the beamlet itself, see [`gauss_parameters`](@ref).
`E_amp` is the field amplitude times the beam radius, `E0 w0` of the source.
"""
struct GaussianBeamletHit{T} <: AbstractBeamletHit{T}
    l0::T
    Δl0::T
    ψ0::T
    p0::Point3{T}
    d0::Point3{T}
    l_seg::T
    segment::GaussianSegment{T}
    E_amp::Complex{T}
    sqrt_proj::T
    w_max::T
end

"""
    AstigmaticGaussianBeamletHit{T} <: AbstractBeamletHit{T}

Stores the segment of an [`AstigmaticGaussianBeamlet`](@ref) that hit the detector as values,
without a reference to the beamlet: `p0` and `d0` are start and direction of the segment,
`l_seg` its length, `l0` the length of the beam (including its parents) up to `p0`, `λ` the
wavelength and `E_pol` the polarization of the chief ray of the segment.
"""
struct AstigmaticGaussianBeamletHit{T} <: AbstractBeamletHit{T}
    l0::T
    # Cache
    p0::Point3{T}
    d0::Point3{T}
    h1::Point3{Complex{T}}
    u1::Point3{Complex{T}}
    h2::Point3{Complex{T}}
    u2::Point3{Complex{T}}
    # Amplitude and Gouy factor of the complex beam area a at the local distance l behind
    # p0, as gouy0 / (√(1 − ρ1 l) √(1 − ρ2 l)): gouy0 is the change of the amplitude with the
    # beam area and the Gouy phase up to p0, followed from the start of the root beam (see
    # `_area_from_root`), ρ1 and ρ2 are the inverse roots of a(l) (see `_area_inverse_roots`)
    gouy0::Complex{T}
    ρ1::Complex{T}
    ρ2::Complex{T}
    k0::T
    Δl::T
    n_eff::T
    E_ref_amp::Complex{T}
    sqrt_proj::T
    w_max::T
    l_seg::T
    λ::T
    E_pol::Point3{Complex{T}}
end

direction(hit::AbstractBeamletHit) = hit.d0
wavelength(hit::GaussianBeamletHit) = hit.segment.λ
wavelength(hit::AstigmaticGaussianBeamletHit) = hit.λ
refractive_index(hit::GaussianBeamletHit) = hit.segment.n
refractive_index(hit::AstigmaticGaussianBeamletHit) = hit.n_eff

hit_point(hit::AbstractBeamletHit) = hit.p0 + hit.l_seg * hit.d0

projection_factor(hit::AbstractBeamletHit) = hit.sqrt_proj^2

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

"Allows continued tracing if set to false"
stop(d::Detector) = d.stop

"""
    push!(d::Detector, hit::AbstractDetectorHit)

Stores the `hit` in the detector `d`, guarded by the `lock` of `d`. Within [`with_hit_buffers`](@ref),
e.g. in the solve of a beam group, the `hit` is collected in the buffer of the current task instead
and reaches `d` when the buffer is flushed.
"""
push!(d::Detector, new::AbstractDetectorHit) = push!(current_hit_sink(), d => new)

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

#=
Hit buffers

The tasks that solve the beams of a group push hits concurrently. Under the lock of a detector each
push waits for the others, e.g. 8 threads were slower than 4 for a detector only. Hence each task
collects its hits and appends them once per detector at its end.
=#

"""
The hits of the current task per [`Detector`](@ref), see [`with_hit_buffers`](@ref). Pushing
`detector => hit` collects the `hit` in the buffer of the `detector`. The buffer of the last
detector is kept at hand, since consecutive hits mostly go to the same detector.
"""
mutable struct HitBuffers
    const buffers::IdDict{Detector, Vector}
    last_detector::Union{Nothing, Detector}
    last_buffer::Vector
    HitBuffers() = new(IdDict{Detector, Vector}(), nothing, [])
end

# Key of the buffers of the current task in the task-local storage
const HIT_BUFFERS_KEY = :BeamletOptics_hit_buffers

"Stores each pushed `detector => hit` directly in its detector, guarded by the `lock` of the detector."
struct LockedHits end

# Where the current task pushes its hits: its buffers within `with_hit_buffers`, otherwise the detectors
current_hit_sink() = get(task_local_storage(), HIT_BUFFERS_KEY, LockedHits())::Union{LockedHits, HitBuffers}

"""
    with_hit_buffers(f)

Calls `f()` such that the hits that the current task pushes to a [`Detector`](@ref) are collected
in a buffer of the task, and appends them to their detectors when `f` returns or throws, one
`append!` per detector. The solve of a beam group calls it in each of its tasks, which then do not
wait for the lock of a detector at every hit. The hits are stored when `with_hit_buffers` returns,
not while `f` runs.
"""
function with_hit_buffers(f)
    buffers = HitBuffers()
    try
        return task_local_storage(f, HIT_BUFFERS_KEY, buffers)
    finally
        for (d, new) in buffers.buffers
            append!(d, new)
        end
    end
end

function push!(::LockedHits, (d, new)::Pair{<:Detector, <:AbstractDetectorHit})
    lock(d.lock) do
        # set data type on first push
        if isnothing(hits(d))
            hits!(d, [new])
            return nothing
        end
        # if the hit type does not match, throws a conversion error
        push!(hits(d), new)
        return nothing
    end
end

function push!(b::HitBuffers, (d, new)::Pair{<:Detector, H}) where {H <: AbstractDetectorHit}
    if b.last_detector !== d
        b.last_buffer = get!(() -> H[], b.buffers, d)
        b.last_detector = d
    end
    push!(b.last_buffer::Vector{H}, new)
    return nothing
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

"""
    GaussianBeamletHit(g::GaussianBeamlet, id)

Hit record of the traced segment `id` of `g`: the data that [`beamlet_hit_field`](@ref)
needs to evaluate the field of that segment at any point near it. [`Detector`](@ref)s store
one per beamlet that hits them; it can also be built for any traced segment, e.g. to sample
a beamlet on a plane of one's own. `l0` is the length of the beam (including its parents)
up to the start of the segment. The projection factor is `1` if the segment ends without an
intersection. The hit keeps no reference to `g`: a new solve of `g` does not change it.
"""
function GaussianBeamletHit(g::GaussianBeamlet{R}, id::Integer) where {R}
    ray = g.chief.rays[id]
    l0 = _length_before(g.chief, id)
    p0 = position(ray)
    d0 = direction(ray)
    segment = GaussianSegment(g, id)
    l_seg = length(ray)
    # Beam radius at the end of the segment, for the auto-limits of the field plots
    w_max = first(gauss_parameters(segment, isfinite(l_seg) ? l_seg : zero(R)))
    return GaussianBeamletHit{R}(l0, _path_excess_before(g, id), _gouy_before(g, id), p0, d0, l_seg,
        segment, electric_field(g) * beam_waist(g), sqrt(projection_factor(ray)), w_max)
end

"""
    AstigmaticGaussianBeamletHit(agb::AstigmaticGaussianBeamlet, id)

Hit record of the traced segment `id` of `agb`: the data that [`beamlet_hit_field`](@ref)
needs to evaluate the field of that segment at any point near it, including the optical
path from the source (parents included), the Gouy factor and the amplitude. The amplitude is
that of the polarization of the chief ray of the segment, times the change with the beam
area within the segments since the start of the root beam; the Gouy phase is followed from
there as well, i.e. a child beamlet continues its parent.
[`Detector`](@ref)s store one per beamlet that hits them; it can also be built for any
traced segment, e.g. to sample a beamlet on a plane of one's own. The projection factor is
`1` if the segment ends without an intersection. The hit keeps no reference to `agb`: a new
solve of `agb` does not change it.
"""
function AstigmaticGaussianBeamletHit(agb::AstigmaticGaussianBeamlet{R}, id::Integer) where {R}
    # Pre-calculate cache
    chief = rays(agb.c)[id]
    p0 = position(chief)
    d0 = direction(chief)
    k0 = 2π / wavelength(chief)
    sqrt_proj = sqrt(projection_factor(chief))

    # Parabasal parameters at segment start (p0)
    h1, u1, h2, u2, _ = parabasal_ray_parameters(agb, p0, id)

    # Change of the amplitude with the beam area and Gouy phase up to p0, from the start of the
    # root beam, which a child continues through its parents
    start = _area_from_root(agb, id, zero(R))
    # Gouy factor at p0 and inverse roots of the area along the hit segment, constant per
    # hit, so that the field loop over the detector pixels does not recompute them
    A = _area_coefficients(h1, u1, h2, u2, d0)
    gouy0 = start.scale * cis(-start.arg / 2)
    ρ1, ρ2 = _area_inverse_roots(A)

    # Complex amplitude of the segment: its chief ray carries what the surfaces on the way did
    # to the field (Fresnel coefficients, polarizers, beamsplitters)
    E_vec = polarization(chief)
    max_idx = argmax(abs.(E_vec))
    E_ref_amp = Complex{R}(norm(E_vec) * cis(angle(E_vec[max_idx])))

    # OPL correction (Δl) and geometric length up to the start of the segment (l0)
    p_parent = agb.parent
    l_parent = isnothing(p_parent) ? 0.0 : length(p_parent)
    opl_parent = isnothing(p_parent) ? 0.0 : optical_path_length(p_parent)

    Δl = opl_parent - l_parent
    l0 = l_parent
    for j in 1:(id - 1)
        ray_j = rays(agb.c)[j]
        Δl += optical_path_length(ray_j) - length(ray_j)
        l0 += length(ray_j)
    end
    # Note: the (n-1)*z term in parabasal_field depends on (z - l0).
    # Propagate semi-axes to the end of the segment for w_max (auto-limits)
    l_seg = isfinite(length(chief)) ? length(chief) : zero(R)
    H1 = h1 + l_seg * u1
    H2 = h2 + l_seg * u2
    w_max = max(norm(H1), norm(H2))
    n_eff = refractive_index(agb, id)

    return AstigmaticGaussianBeamletHit{R}(
        R(l0), p0, d0, h1, u1, h2, u2, gouy0, ρ1, ρ2, k0, R(Δl),
        n_eff, E_ref_amp, sqrt_proj, w_max, length(chief), wavelength(chief), E_vec)
end

# Geometric length of `beam` (including its parents) up to the start of ray `id`
function _length_before(beam::Beam{T}, id::Integer) where {T}
    p = AbstractTrees.parent(beam)
    l = isnothing(p) ? zero(T) : length(p)
    for j in 1:(id - 1)
        l += length(rays(beam)[j])
    end
    return l
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

# include eval and helper functions
include("DetectorUtils.jl")
include("SpotDiagram.jl")
include("FieldCalculation.jl")
