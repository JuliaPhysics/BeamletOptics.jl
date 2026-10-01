# Couples BeamletOptics (BMO) to OpticsBase's single exchange format, `PlaneField`
# (tangential E and H sampled on a plane, see `OpticsBase.jl`). Two directions:
#
# - `OpticsBase.PlaneField(detector::BMO.Detector; size, spacing, origin, axes)`: samples
#   the coherent sum of a BMO `Detector`'s beamlet hits into a `PlaneField`.
# - `BeamletOptics.WavefrontBeamletDecomposition(f::OpticsBase.PlaneField; kwargs...)`:
#   the return path, decomposing `forward(f)` into `AstigmaticGaussianBeamlet`s on `f`'s
#   plane, travelling along `f`'s normal `n`.
#
# The detector port convention (plane origin/axes from a BMO `Detector`) and the
# equivalence between BMO's own field summation and an externally sampled field were
# verified against OpticsBase's pre-rewrite BMO glue (`OpticsBaseBeamletOpticsExt` and
# its integration tests, commit 00509b3 of the OpticsBase repository) rather than
# re-derived from scratch; see the docstrings below for the resulting conventions.
module BeamletOpticsOpticsBaseExt

using BeamletOptics
import BeamletOptics: Detector, GaussianBeamletHit, AstigmaticGaussianBeamletHit,
                       position, direction, orientation, wavelength, refractive_index,
                       beamlet_hit_field, beamlet_hit_polarization,
                       WavefrontBeamletDecomposition, GaussianModeDecomposition,
                       translate_to3d!, _with_progress,
                       _tick!
const BMO = BeamletOptics

import OpticsBase
using OpticsBase: PlaneField, VACUUM_IMPEDANCE, forward, reference_phase
using LinearAlgebra: cross, dot, norm
using StaticArrays: SVector, SMatrix

const _PREFIX = "BeamletOptics/OpticsBase extension"

# -----------------------------------------------------------------------------------------
# BMO Detector -> OpticsBase.PlaneField

"""
    OpticsBase.PlaneField(detector::BeamletOptics.Detector;
        size::NTuple{2,Int}, spacing::NTuple{2,Real},
        origin = BeamletOptics.position(detector), axes = <default detector frame, see below>)

Samples the coherent field of the beamlets hitting a BeamletOptics (BMO) `Detector` into
an [`OpticsBase.PlaneField`](@ref). Call it after `solve_system!` and before moving or
emptying the detector. Requires `using OpticsBase` (package extension).

# Arguments

- `detector`: a BMO `Detector` whose hits are all [`GaussianBeamletHit`](@ref
  BeamletOptics.GaussianBeamletHit) or all [`AstigmaticGaussianBeamletHit`](@ref
  BeamletOptics.AstigmaticGaussianBeamletHit) (i.e. Gaussian beamlet sources; pure-ray
  hits have no scaled field in BMO and throw an `ArgumentError`, as does an empty
  detector). All hits must share the vacuum wavelength and the refractive index of the
  medium in front of the detector.
- `size` (keyword, required): number of samples `(nx, ny)` along the plane's `u`, `v`
  axes.
- `spacing` (keyword, required): sample spacings `(Δu, Δv)` in \\[m\\], as in
  [`OpticsBase.PlaneField`](@ref).
- `origin` (keyword): center of the plane in global coordinates in \\[m\\]. Default: the
  detector position.
- `axes` (keyword): `3×3` matrix with columns `u`, `v`, `n` (see
  [`OpticsBase.PlaneField`](@ref)). Default: `n = −orientation(detector)[:, 2]` (BMO's
  detector normal points against the beam; `n` points downstream, i.e. the direction the
  light was travelling in when it hit the detector), `u = −orientation(detector)[:, 1]`
  (BMO's local detector x axis) and `v = n × u` (BMO's local −z axis), matching the
  default port of BMO's detector hit readers.

# Conventions

- `E` is the physical field: the coherent sum of the beamlet fields at the plane's
  sample points, projected onto `u`, `v`. Unlike `electric_field(detector)`, it carries
  no `√|cos θ|` projection factor; for oblique beams `power(f)` is the full beam power
  crossing the plane, not the detector-area power density integral. `H` is the coherent sum, **over hits**, of `(n_medium/Z₀) dᵢ × Eᵢ`, with
  `dᵢ` the propagation direction of the chief/stigmatic ray of hit `i` and `Eᵢ` that
  hit's own 3D field contribution at the sample point (not the paraxial `(n/Z₀) n × E`
  built from the *total* field, which would be exact only for normal incidence); this
  matches the surface-equivalence-theorem content of `PlaneField` even for beams that hit
  the plane off-normal. `R = Inf` (no reference sphere removed: BMO's beamlet field
  already carries the full spatial phase without needing one).
- `Eᵢ` needs a 3D direction for each hit's otherwise-scalar contribution:
  [`AstigmaticGaussianBeamletHit`](@ref BeamletOptics.AstigmaticGaussianBeamletHit)s use
  their chief ray's own polarization ([`beamlet_hit_polarization`](@ref
  BeamletOptics.beamlet_hit_polarization)). Stigmatic
  [`GaussianBeamletHit`](@ref BeamletOptics.GaussianBeamletHit)s carry no polarization in
  BMO (the model is scalar); by convention, matching `OpticsBase`'s own rule for a scalar
  field (`PlaneField(E::AbstractMatrix, ...)` puts a scalar field entirely into `Eu`),
  such a hit's contribution is polarized along the plane's `u` axis, made transverse to
  the hit direction (along `v` if `u` is parallel to it). At normal incidence this is
  exactly `u`.
- `n` (medium refractive index) is that of the hits; `λ` is their shared wavelength.
- Time convention `exp(−iωt)`, as in BMO and `OpticsBase`.

Throws an `ArgumentError` if the detector has no hits, if the hits are not beamlet hits,
or if the hits do not share wavelength and refractive index.
"""
function OpticsBase.PlaneField(detector::Detector;
        size::NTuple{2, Int}, spacing::NTuple{2, <:Real},
        origin = position(detector), axes = _default_plane_axes(detector),
        progress::Bool = true)
    hs = _beamlet_hits(detector)
    λ, n_medium = _common_wavelength_index(hs)
    return _sample_plane_field(hs, λ, n_medium, size, spacing, origin, axes, progress)
end

function _default_plane_axes(detector::Detector)
    R = orientation(detector)
    n = -SVector{3}(R[:, 2])
    u = -SVector{3}(R[:, 1])
    v = cross(n, u)
    return SMatrix{3, 3}(hcat(u, v, n))
end

function _beamlet_hits(detector::Detector)
    hs = BeamletOptics.hits(detector)
    (hs === nothing || isempty(hs)) &&
        throw(ArgumentError("$_PREFIX: the detector has no hits; run `solve_system!` first"))
    return hs
end
_beamlet_hits_error(hs) = throw(ArgumentError(
    "$_PREFIX: PlaneField needs Gaussian beamlet hits (GaussianBeamletHit or AstigmaticGaussianBeamletHit), got $(eltype(hs))"))

function _common_wavelength_index(hs::AbstractVector{<:Union{GaussianBeamletHit, AstigmaticGaussianBeamletHit}})
    λ = wavelength(first(hs))
    n_medium = refractive_index(first(hs))
    tol = sqrt(eps(float(typeof(λ))))
    for (i, h) in enumerate(hs)
        isapprox(wavelength(h), λ; rtol = tol) ||
            throw(ArgumentError("$_PREFIX: all hits must share the wavelength, got $(wavelength(h)) m for hit $i and $λ m for hit 1"))
        isapprox(refractive_index(h), n_medium; rtol = tol) ||
            throw(ArgumentError("$_PREFIX: all hits must share the refractive index, got $(refractive_index(h)) for hit $i and $n_medium for hit 1"))
    end
    return λ, n_medium
end
_common_wavelength_index(hs::AbstractVector) = _beamlet_hits_error(hs)

# Per-hit 3D field direction: the chief ray's own polarization for astigmatic hits. A
# stigmatic hit has no polarization in BMO; it gets the plane's `u` axis made transverse
# to the hit direction (`v` if `u` is parallel to it), so that E ⟂ d and H = (n/Z₀) d × E
# form a consistent plane wave.
function _hit_field_direction(hit::GaussianBeamletHit, u, v)
    d = direction(hit)
    e = norm(cross(d, u)) > 1e-6 ? u : v
    p = e - dot(d, e) * d
    return p / norm(p)
end
_hit_field_direction(hit::AstigmaticGaussianBeamletHit, u, v) = beamlet_hit_polarization(hit)

function _sample_plane_field(hs, λ, n_medium, (nx, ny), (Δu, Δv), origin, axes, progress)
    T = typeof(float(λ))
    u = SVector{3, T}(axes[:, 1])
    v = SVector{3, T}(axes[:, 2])
    r0 = SVector{3, T}(origin)

    us = ((0:(nx - 1)) .- nx ÷ 2) .* Δu
    vs = ((0:(ny - 1)) .- ny ÷ 2) .* Δv

    E = zeros(Complex{T}, nx, ny, 2)
    H = zeros(Complex{T}, nx, ny, 2)

    hit_dirs = [SVector{3, T}(direction(h)) for h in hs]
    hit_pols = [SVector{3, Complex{T}}(_hit_field_direction(h, u, v)) for h in hs]

    Y = T(n_medium / VACUUM_IMPEDANCE)

    _with_progress(progress, ny, "PlaneField: ") do prog
        Threads.@threads for j in eachindex(vs)
            v_grid = vs[j]
            @inbounds for i in eachindex(us)
                p = r0 + us[i] * u + v_grid * v
                E3 = zero(SVector{3, Complex{T}})
                C3 = zero(SVector{3, Complex{T}})    # Σᵢ dᵢ × Eᵢ
                for (hit, d, pol) in zip(hs, hit_dirs, hit_pols)
                    Ei = beamlet_hit_field(hit, p) .* pol
                    E3 += Ei
                    C3 += cross(d, Ei)
                end
                H3 = Y .* C3
                # `dot` conjugates its first argument: the real axis goes first
                E[i, j, 1] = dot(u, E3)
                E[i, j, 2] = dot(v, E3)
                H[i, j, 1] = dot(u, H3)
                H[i, j, 2] = dot(v, H3)
            end
            _tick!(prog)
        end
    end
    return PlaneField(E, H, (T(Δu), T(Δv)), r0, SMatrix{3, 3, T}(axes), T(λ); n = T(n_medium))
end

# -----------------------------------------------------------------------------------------
# OpticsBase.PlaneField -> BeamletOptics.WavefrontBeamletDecomposition

# Physical field of the part of `f` travelling along +n. A field without forward light
# (e.g. a purely backward wave, whose forward part is rounding noise) cannot be decomposed.
function _forward_field(f::OpticsBase.PlaneField)
    E = forward(f).E
    norm(E) > sqrt(eps(real(eltype(E)))) * norm(f.E) || throw(ArgumentError(
        "$_PREFIX: the PlaneField has no forward-travelling light to decompose"))
    return E .* reference_phase(f)
end

"""
    BeamletOptics.WavefrontBeamletDecomposition(f::OpticsBase.PlaneField; kwargs...)

Decomposes the forward-travelling part of a [`PlaneField`](@ref OpticsBase.PlaneField)
`f` into an [`AstigmaticBeamGroup`](@ref BeamletOptics.AstigmaticBeamGroup) of
`AstigmaticGaussianBeamlet`s, placed on `f`'s plane and travelling along `f`'s normal
`n`. This is the return path for `OpticsBase.PlaneField(detector; ...)`: passing the
result back into this constructor and tracing the resulting beam group reproduces the
original field (round trip), up to beamlet-grid discretization error.

# Arguments

- `f`: the field to decompose. Only `forward(f)` (the part travelling along `+n`) is
  used; the backward part is dropped (a `PlaneField` has no ready single propagation
  direction for its backward part in this constructor, which needs one direction per
  call).
- `kwargs...`: forwarded to [`WavefrontBeamletDecomposition`](@ref
  BeamletOptics.WavefrontBeamletDecomposition)`(x, y, Eu, Ev, dir, λ; kwargs...)`
  (`threshold`, `overlap`, `randomize_axes`, `rng`); `basis` is fixed to `f`'s own `(u,
  v)` and cannot be overridden.

# Conventions

- The physical field used is `forward(f).E .* reference_phase(forward(f))`: the stored
  array with the reference-sphere phase (if any) and the local-plane-wave forward/backward
  split both applied, in \\[V/m\\] peak amplitude, matching the units
  `AstigmaticGaussianBeamlet`'s `E0` expects.
- The beam group's sampling grid is `f`'s own `(u, v)` grid ([`coordinates`](@ref
  OpticsBase.coordinates)), placed at `f.origin` with local axes `f.axes[:, 1:2]`
  (`translate_to3d!` after construction; [`WavefrontBeamletDecomposition`](@ref
  BeamletOptics.WavefrontBeamletDecomposition) itself always centers the undecomposed
  group on the global origin).
- `f.n` (the medium refractive index) and `f.R` are not carried over:
  `AstigmaticGaussianBeamlet` traces in vacuum from the plane (BMO assigns the medium
  index by tracing through the scene), and `f`'s reference-sphere phase is folded into
  the sampled amplitude/phase before decomposition, so no curvature information is lost.

Throws an `ArgumentError` if `f` has no forward-travelling light (e.g. a purely backward
wave), and whatever `WavefrontBeamletDecomposition(x, y, Eu, Ev, dir, λ)` throws.
"""
function BeamletOptics.WavefrontBeamletDecomposition(f::OpticsBase.PlaneField; kwargs...)
    haskey(kwargs, :basis) &&
        throw(ArgumentError("$_PREFIX: `basis` is fixed to the PlaneField's own (u, v) axes and cannot be overridden"))
    phys = _forward_field(f)
    Eu = phys[:, :, 1]
    Ev = phys[:, :, 2]

    x = collect(OpticsBase.coordinates(f, 1))
    y = collect(OpticsBase.coordinates(f, 2))
    n = SVector{3}(f.axes[:, 3])
    u = SVector{3}(f.axes[:, 1])
    v = SVector{3}(f.axes[:, 2])

    group = WavefrontBeamletDecomposition(x, y, Eu, Ev, n, f.λ; basis = (u, v), kwargs...)
    translate_to3d!(group, SVector{3}(f.origin))
    return group
end

# -----------------------------------------------------------------------------------------
# OpticsBase.PlaneField -> BeamletOptics.GaussianModeDecomposition

"""
    BeamletOptics.GaussianModeDecomposition(f::OpticsBase.PlaneField)

Fits a single `AstigmaticGaussianBeamlet` to the forward-travelling part of a
[`PlaneField`](@ref OpticsBase.PlaneField), e.g. the output of a single-mode fiber. The
beamlet starts on `f`'s plane at the centroid of the field; see
[`GaussianModeDecomposition`](@ref BeamletOptics.GaussianModeDecomposition) for the fit.
Use it for fields close to one Gaussian mode that are too small for
[`WavefrontBeamletDecomposition`](@ref BeamletOptics.WavefrontBeamletDecomposition)
(which needs beamlets much larger than λ).

`optical_power(beamlet) / OpticsBase.power(f)` is the fraction of the power captured by
the Gaussian mode. The field used is `forward(f).E .* reference_phase(f)`. `f` must lie in
vacuum or air (`f.n == 1`), where BMO starts beams, and must contain forward-travelling
light; an `ArgumentError` is thrown otherwise.
"""
function BeamletOptics.GaussianModeDecomposition(f::OpticsBase.PlaneField)
    f.n == 1 || throw(ArgumentError(
        "$_PREFIX: the PlaneField must lie in vacuum or air (n = 1), got n = $(f.n)"))
    E = _forward_field(f)
    u, v, n = (SVector{3}(f.axes[:, i]) for i in 1:3)
    return GaussianModeDecomposition(collect(OpticsBase.coordinates(f, 1)),
        collect(OpticsBase.coordinates(f, 2)), E[:, :, 1], E[:, :, 2], n, f.λ;
        basis = (u, v), origin = SVector{3}(f.origin))
end

end # module
