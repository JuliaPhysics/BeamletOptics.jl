"""
    AbstractDetectorHit

Abstract supertype for all [`Detector`](@ref) hit records.

A hit of this type encapsulates the interaction of a ray or beamlet with the [`Detector`](@ref)
surface. Instead of directly accumulating values (e.g. optical power), the hit
object stores all relevant information about the incident field for a posteriori
evaluation (such as plotting, power integration, or polarization analysis).

# Implementation reqs.

## Functions

Subtypes of `AbstractDetectorHit` must implement the following functions:

- `electric_field(hit, p)`: the field of the hit at a global point `p`, see
  [`electric_field(::BeamletOptics.AbstractDetectorHit, ::AbstractVector)`](@ref)
- `direction(hit)`: the direction of propagation at the hit
- `hit_point(hit)`: the global point where the ray or chief ray crossed the detector, for
  [`ray_hits`](@ref)
- `refractive_index(hit)`: the refractive index in front of the detector
- `intensity(hits, p, n̂)`: the irradiance of a vector of such hits at a global point `p`
  through a plane with the normal `n̂`, see [`intensity`](@ref)
- `calc_local_lims(pd, hits, crop_factor::Real; center::AbstractCenterAlgorithm = Centroid())`:
  the window of the default [`RegularGrid`](@ref) around the hits
- `calc_local_lims(pd, hits, ::Nothing; center::AbstractCenterAlgorithm = Centroid())`: the
  same window with the default `crop_factor` of the hit type
- `calc_local_pos(pd, hits)`: the local (x, z) positions of the hits for [`spot_diagram`](@ref)

The readouts of a [`Detector`](@ref) ([`electric_field`](@ref), [`intensity`](@ref),
[`optical_power`](@ref)) need all of them except `hit_point` and `calc_local_pos`.

!!! info
    Note that new subtypes of `AbstractDetectorHit` must be manually added to the `Detector`
    definition to be eligible.

Concrete subtypes include:

- [`RayHit`](@ref) — unpolarized geometric ray
- [`PolarizedRayHit`](@ref) — ray carrying a polarization state
- [`GaussianBeamletHit`](@ref) — single Gaussian beamlet
- [`AstigmaticGaussianBeamletHit`](@ref) — single astigmatic Gaussian beamlet with a vector field
"""
abstract type AbstractDetectorHit end

hit_point(hit::AbstractDetectorHit) = error(lazy"hit_point not defined for $(typeof(hit))")

"""
    AbstractRayHit{T} <: AbstractDetectorHit

Abstract supertype for detector hits produced by [`AbstractRay`](@ref)s.
Currently the following concrete types are implemented:

- [`RayHit`](@ref)
- [`PolarizedRayHit`](@ref)

A hit stores the values of the ray at the moment it hit the detector, not the ray itself: a
new solve of the beam or a move of any object does not change it.

# Implementation reqs.

## Fields

The default getters read these fields:

- `p`: the R³ point of intersection
- `dir`: the direction of the ray
- `opl`: the optical path length of the beam up to the hit point (incl. its parents)
- `λ`: the vacuum wavelength of the ray
- `n`: the refractive index in front of the detector
- `proj`: the projection factor, see below

## Functions

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
direction(hit::AbstractRayHit) = hit.dir
optical_path_length(hit::AbstractRayHit) = hit.opl
wavelength(hit::AbstractRayHit) = hit.λ
wavenumber(hit::AbstractRayHit) = wavenumber(wavelength(hit))
refractive_index(hit::AbstractRayHit) = hit.n
projection_factor(hit::AbstractRayHit) = hit.proj

"""
    RayHit(ray::Ray, opl)

Stores the hit of a [`Ray`](@ref) as values: hit point `p`, direction `dir`, optical path length `opl`
of the beam up to the hit point, wavelength `λ`, refractive index `n` and projection factor `proj`.
"""
struct RayHit{T} <: AbstractRayHit{T}
    p::Point3{T}
    dir::Point3{T}
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
vector `E0` of the ray in global coordinates.
"""
struct PolarizedRayHit{T} <: AbstractRayHit{T}
    p::Point3{T}
    dir::Point3{T}
    opl::T
    λ::T
    n::T
    proj::T
    E0::Point3{Complex{T}}
end

PolarizedRayHit(ray::PolarizedRay{T}, opl::Real) where {T} =
    PolarizedRayHit{T}(hit_point(ray), direction(ray), opl, wavelength(ray), refractive_index(ray),
        projection_factor(ray), polarization(ray))

polarization(hit::PolarizedRayHit) = hit.E0

"""
    AbstractBeamletHit

Stores beamlet hits. Currently implemented:

- [`GaussianBeamletHit`](@ref)
- [`AstigmaticGaussianBeamletHit`](@ref)
"""
abstract type AbstractBeamletHit{T} <: AbstractDetectorHit end

"""
    GaussianBeamletHit{T} <: AbstractBeamletHit{T}

Stores the segment of a [`GaussianBeamlet`](@ref) that hit the detector as values, without a
reference to the beamlet: `p0` and `dir` are start and direction of the segment, `l_seg` its
length, `l0` the length of the beam (including its parents) up to `p0`, `Δl0` the optical minus
the geometric path and `ψ0` the Gouy phase that the beam has collected up to there. `segment`
holds the waist and divergence rays of the segment (see [`BeamletOptics.GaussianSegment`](@ref)),
from which the field follows as from the beamlet itself, see [`gauss_parameters`](@ref).
`E_amp` is the reference amplitude times the beam radius, `E0 w0` of the source, which is normalized to the
power of the beamlet: the physical field in a medium of index `n` is smaller by `√n`, see
[`electric_field(::BeamletOptics.AbstractDetectorHit, ::AbstractVector)`](@ref).
"""
struct GaussianBeamletHit{T} <: AbstractBeamletHit{T}
    l0::T
    Δl0::T
    ψ0::T
    p0::Point3{T}
    dir::Point3{T}
    l_seg::T
    segment::GaussianSegment{T}
    E_amp::Complex{T}
    w_max::T
end

"""
    AstigmaticGaussianBeamletHit{T} <: AbstractBeamletHit{T}

Stores the segment of an [`AstigmaticGaussianBeamlet`](@ref) that hit the detector as values,
without a reference to the beamlet: `p0` and `dir` are start and direction of the segment,
`l_seg` its length, `l0` the length of the beam (including its parents) up to `p0`, `Δl0` the
optical minus the geometric path up to there, `n` the refractive index of the segment, `λ` the
wavelength, `k0` the vacuum wavenumber and `E0` the polarization of the chief ray of the
segment. `h1`, `u1`, `h2` and `u2` are the parabasal ray parameters at `p0`. The amplitude and
Gouy factor of the complex beam area `a` at the local distance `l` behind `p0` is
`gouy0 / (√(1 − ρ1 l) √(1 − ρ2 l))`: `gouy0` is the change of the amplitude with the beam area
and the Gouy phase up to `p0`, followed from the start of the root beam (see
`_area_from_root`), and `ρ1` and `ρ2` are the inverse roots of `a(l)` (see
`area_inverse_roots`).
"""
struct AstigmaticGaussianBeamletHit{T} <: AbstractBeamletHit{T}
    l0::T
    p0::Point3{T}
    dir::Point3{T}
    h1::Point3{Complex{T}}
    u1::Point3{Complex{T}}
    h2::Point3{Complex{T}}
    u2::Point3{Complex{T}}
    gouy0::Complex{T}
    ρ1::Complex{T}
    ρ2::Complex{T}
    k0::T
    Δl0::T
    n::T
    w_max::T
    l_seg::T
    λ::T
    E0::Point3{Complex{T}}
end

direction(hit::AbstractBeamletHit) = hit.dir
wavelength(hit::GaussianBeamletHit) = hit.segment.λ
wavelength(hit::AstigmaticGaussianBeamletHit) = hit.λ
refractive_index(hit::GaussianBeamletHit) = hit.segment.n
refractive_index(hit::AstigmaticGaussianBeamletHit) = hit.n

hit_point(hit::AbstractBeamletHit) = hit.p0 + hit.l_seg * hit.dir

"""
    GaussianBeamletHit(g::GaussianBeamlet, id)

Hit record of the traced segment `id` of `g`: the data that [`electric_field`](@ref)`(hit, p)`
needs to evaluate the field of that segment at any point near it. [`Detector`](@ref)s store
one per beamlet that hits them; it can also be built for any traced segment, e.g. to sample
a beamlet on a plane of one's own. `l0` is the length of the beam (including its parents)
up to the start of the segment.
"""
function GaussianBeamletHit(g::GaussianBeamlet{R}, id::Integer) where {R}
    ray = rays(g)[id]
    l0 = length_before(g, id)
    p0 = position(ray)
    dir = direction(ray)
    segment = GaussianSegment(g, id)
    l_seg = length(ray)
    # Beam radius at the segment end for the plot limits
    w_max = first(gauss_parameters(segment, isfinite(l_seg) ? l_seg : zero(R)))
    return GaussianBeamletHit{R}(l0, path_excess_before(g, id), _gouy_before(g, id), p0, dir, l_seg,
        segment, electric_field(g) * beam_waist(g), w_max)
end

"""
    AstigmaticGaussianBeamletHit(agb::AstigmaticGaussianBeamlet, id)

Hit record of the traced segment `id` of `agb`: the data that [`electric_field`](@ref)`(hit, p)`
needs to evaluate the field of that segment at any point near it, including the optical
path from the source (parents included), the Gouy factor and the amplitude. The amplitude is
that of the polarization of the chief ray of the segment, times the change with the beam
area within the segments since the start of the root beam; the Gouy phase is followed from
there as well, i.e. a child beamlet continues its parent.
[`Detector`](@ref)s store one per beamlet that hits them; it can also be built for any
traced segment, e.g. to sample a beamlet on a plane of one's own. The hit keeps no reference
to `agb`: a new solve of `agb` does not change it.
"""
function AstigmaticGaussianBeamletHit(agb::AstigmaticGaussianBeamlet{R}, id::Integer) where {R}
    # Chief ray of the segment
    chief = rays(agb)[id]
    p0 = position(chief)
    dir = direction(chief)
    k0 = wavenumber(chief)

    # Parabasal parameters at segment start (p0)
    h1, u1, h2, u2, _ = parabasal_ray_parameters(agb, p0, id)

    # Area change and Gouy phase up to p0 from the root beam
    start = _area_from_root(agb, id, zero(R))
    # Gouy factor at p0 and inverse area roots along the segment
    A = area_coefficients(h1, u1, h2, u2, dir)
    gouy0 = start.scale * cis(-start.arg / 2)
    ρ1, ρ2 = area_inverse_roots(A)

    # Amplitude from the chief ray polarization
    E0 = polarization(chief)

    # Geometric length and path excess up to p0
    l0 = length_before(agb, id)
    Δl0 = path_excess_before(agb, id)

    # Semi-axes at the segment end for the plot limits
    l_seg = isfinite(length(chief)) ? length(chief) : zero(R)
    H1 = h1 + l_seg * u1
    H2 = h2 + l_seg * u2
    w_max = max(norm(H1), norm(H2))
    n = refractive_index(agb, id)

    return AstigmaticGaussianBeamletHit{R}(
        R(l0), p0, dir, h1, u1, h2, u2, gouy0, ρ1, ρ2, k0, R(Δl0),
        n, w_max, length(chief), wavelength(chief), E0)
end
