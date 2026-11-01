# BeamletOptics <-> OpticsBase.PlaneField
module BeamletOpticsOpticsBaseExt

using BeamletOptics
import BeamletOptics: Detector, AbstractDetectorHit, GaussianBeamletHit,
                       AstigmaticGaussianBeamletHit, GaussianBeamlet, AstigmaticGaussianBeamlet,
                       position, direction, wavelength, refractive_index, electric_field,
                       isparallel3d, local_axes,
                       WavefrontBeamletDecomposition, GaussianModeDecomposition,
                       hits, _sample

import OpticsBase
using OpticsBase: PlaneField, VACUUM_IMPEDANCE, forward
using LinearAlgebra: cross, dot, norm
using StaticArrays: SVector, SMatrix

# -----------------------------------------------------------------------------------------
# BMO Detector -> OpticsBase.PlaneField

"""
    OpticsBase.PlaneField(detector::BeamletOptics.Detector;
        size::NTuple{2,Int}, spacing::NTuple{2,Real},
        origin = BeamletOptics.position(detector), axes = <default detector frame, see below>,
        progress = true)

Samples the coherent field of the beamlets hitting a BeamletOptics (BMO) `Detector` into
an [`OpticsBase.PlaneField`](https://StackEnjoyer.github.io/OpticsBase.jl/dev/). Call it after `solve_system!`. The hits keep the state of
their solve: solving again or moving other objects does not change them, moving the detector
empties it. Requires `using OpticsBase` (package extension).

# Arguments

- `detector`: a BMO `Detector` whose hits are all [`GaussianBeamletHit`](@ref
  BeamletOptics.GaussianBeamletHit) or all [`AstigmaticGaussianBeamletHit`](@ref
  BeamletOptics.AstigmaticGaussianBeamletHit) (i.e. Gaussian beamlet sources). All hits must
  share the vacuum wavelength and the refractive index of the medium in front of the
  detector.
- `size` (keyword, required): number of samples `(nx, ny)` along the plane's `u`, `v`
  axes.
- `spacing` (keyword, required): sample spacings `(Δu, Δv)` in \\[m\\], as in
  [`OpticsBase.PlaneField`](https://StackEnjoyer.github.io/OpticsBase.jl/dev/).
- `origin` (keyword): center of the plane in global coordinates in \\[m\\]. Default: the
  detector position.
- `axes` (keyword): `3×3` matrix with columns `u`, `v`, `n` (see
  [`OpticsBase.PlaneField`](https://StackEnjoyer.github.io/OpticsBase.jl/dev/)). Default:
  `n = −direction(detector)` (downstream, the direction of the light that hit the detector),
  `u` the local detector x axis of BMO and `v = n × u`, the default port of BMO's detector
  readouts.
- `progress` (keyword): show a progress bar. Default: `true`.

# Conventions

- `E` is the physical field: the coherent sum of the beamlet fields at the plane's
  sample points, projected onto `u`, `v`. Like `electric_field(detector)`, it carries
  no projection factor; for oblique beams `power(f)` is the full beam power crossing the
  plane, as `optical_power(detector)`.
- `H` is the coherent sum over hits of `(n_medium/Z₀) dᵢ × Eᵢ`, with `dᵢ` the propagation
  direction of the chief/stigmatic ray of hit `i` and `Eᵢ` that hit's own 3D field at the
  sample point.
- `R = Inf`: no reference sphere is removed.
- `Eᵢ` is the 3D field of each hit:
  [`AstigmaticGaussianBeamletHit`](@ref BeamletOptics.AstigmaticGaussianBeamletHit)s
  carry the polarization of their chief ray, `Eᵢ` is their field vector
  ([`electric_field`](@ref BeamletOptics.electric_field)`(hit, p)`). A stigmatic
  [`GaussianBeamletHit`](@ref BeamletOptics.GaussianBeamletHit) (scalar model) is polarized
  along the plane's `u` axis made transverse to the hit direction (along `v` if `u` is
  parallel to it). At normal incidence this is exactly `u`.
- `n` (medium refractive index) is that of the hits; `λ` is their shared wavelength.
- Time convention `exp(−iωt)`, as in BMO and `OpticsBase`.

Throws an `ErrorException` if the detector has no hits, since a `PlaneField` takes `λ` and `n`
from them, and an `ArgumentError` if the hits are not beamlet hits (e.g. ray hits) or do not
share wavelength and refractive index.
"""
function OpticsBase.PlaneField(detector::Detector;
        size::NTuple{2, Int}, spacing::NTuple{2, <:Real},
        origin = position(detector), axes = default_plane_axes(detector),
        progress::Bool = true)
    hs = hits(detector)
    isnothing(hs) && throw(ErrorException("No hits available on detector."))
    return PlaneField(hs; size, spacing, origin, axes, progress)
end

"""
    OpticsBase.PlaneField(beam::BeamletOptics.AstigmaticGaussianBeamlet, id::Integer;
        size::NTuple{2,Int}, spacing::NTuple{2,Real}, origin, axes, progress = false)
    OpticsBase.PlaneField(beam::BeamletOptics.GaussianBeamlet, id::Integer; kwargs...)

Samples the field of segment `id` of a traced beamlet `beam` on a plane, as
`PlaneField(detector; ...)` does for the beamlets that hit a detector.

`size` and `spacing` are as for `PlaneField(detector; ...)`; `origin` (in \\[m\\]) and
`axes` (columns `u`, `v`, `n`) place the plane and have no default. The field is that of
the segment's straight continuation through the plane (no further interaction), including
the optical path from the source (parents included), and follows the conventions of
`PlaneField(detector; ...)`: no projection factor, `H` from the beamlet direction, a scalar
`GaussianBeamlet` polarized along `u`, `R = Inf`, `n` the refractive index of the segment.

Both beamlet models can be sampled on any of their segments, also within a medium.

!!! warning
    The scalar `GaussianBeamlet` takes no Fresnel losses at surfaces: within a medium its
    field carries the full power of the beamlet, while that of an `AstigmaticGaussianBeamlet`
    is reduced by the transmission of each surface it passed.
"""
OpticsBase.PlaneField(beam::AstigmaticGaussianBeamlet, id::Integer; kwargs...) =
    PlaneField([AstigmaticGaussianBeamletHit(beam, id)]; kwargs...)

OpticsBase.PlaneField(beam::GaussianBeamlet, id::Integer; kwargs...) =
    PlaneField([GaussianBeamletHit(beam, id)]; kwargs...)

"""
    OpticsBase.PlaneField(hits::AbstractVector{<:Union{GaussianBeamletHit, AstigmaticGaussianBeamletHit}};
        size::NTuple{2,Int}, spacing::NTuple{2,Real}, origin, axes, progress = false)

Samples the coherent field of a vector of beamlet hits on a plane. Build the hits with
`BeamletOptics.GaussianBeamletHit(beam, id)` or `BeamletOptics.AstigmaticGaussianBeamletHit(beam, id)`
(public, not exported).

# Arguments

- `hits`: all [`GaussianBeamletHit`](@ref BeamletOptics.GaussianBeamletHit)s or all
  [`AstigmaticGaussianBeamletHit`](@ref BeamletOptics.AstigmaticGaussianBeamletHit)s, with
  a shared vacuum wavelength and refractive index.
- `size` (keyword, required): number of samples `(nx, ny)` along `u`, `v`.
- `spacing` (keyword, required): sample spacings `(Δu, Δv)` in \\[m\\].
- `origin` (keyword, required): center of the plane in global coordinates in \\[m\\].
- `axes` (keyword, required): `3×3` matrix with columns `u`, `v`, `n`.
- `progress` (keyword): show a progress bar. Default: `false`.

The field follows the conventions of `PlaneField(detector; ...)`.

Throws an `ArgumentError` if `hits` is empty, if the hits are not beamlet hits, or if they
do not share wavelength and refractive index.
"""
function OpticsBase.PlaneField(
        hits::AbstractVector{<:Union{GaussianBeamletHit, AstigmaticGaussianBeamletHit}};
        size::NTuple{2, Int}, spacing::NTuple{2, <:Real}, origin, axes,
        progress::Bool = false)
    isempty(hits) && throw(ArgumentError("PlaneField: the vector of hits is empty"))
    λ, n_medium = common_wavelength_index(hits)
    T = typeof(float(λ))
    u = SVector{3, T}(axes[:, 1])
    v = SVector{3, T}(axes[:, 2])
    r0 = SVector{3, T}(origin)

    # Sampling grid
    (nx, ny), (Δu, Δv) = size, spacing
    us = ((0:(nx - 1)) .- nx ÷ 2) .* Δu
    vs = ((0:(ny - 1)) .- ny ÷ 2) .* Δv
    points = [SVector(x, y) for x in us, y in vs]

    hit_dirs = [SVector{3, T}(direction(h)) for h in hits]
    hit_pols = [hit_field_direction(h, u, v) for h in hits]
    Y = T(n_medium / VACUUM_IMPEDANCE)

    values = _sample(hits, r0, u, v, points, progress) do hs, p
        E3 = zero(SVector{3, Complex{T}})
        C3 = zero(SVector{3, Complex{T}})    # Σᵢ dᵢ × Eᵢ
        for (hit, d, pol) in zip(hs, hit_dirs, hit_pols)
            Ei = hit_field(hit, p, pol)
            E3 += Ei
            C3 += cross(d, Ei)
        end
        H3 = Y .* C3
        # Real axis first
        return SVector(dot(u, E3), dot(v, E3), dot(u, H3), dot(v, H3))
    end
    E = Array{Complex{T}}(undef, nx, ny, 2)
    H = similar(E)
    for k in CartesianIndices(values)
        E[k, 1], E[k, 2], H[k, 1], H[k, 2] = values[k]
    end
    return PlaneField(E, H, (T(Δu), T(Δv)), r0, SMatrix{3, 3, T}(axes), T(λ); n = T(n_medium))
end

OpticsBase.PlaneField(hits::AbstractVector{<:AbstractDetectorHit}; kwargs...) =
    throw(ArgumentError("PlaneField: needs Gaussian beamlet hits (GaussianBeamletHit or AstigmaticGaussianBeamletHit), got $(eltype(hits))"))

function default_plane_axes(detector::Detector)
    u = SVector{3}(first(local_axes(detector)))
    n = -SVector{3}(direction(detector))
    v = cross(n, u)
    return SMatrix{3, 3}(hcat(u, v, n))
end

function common_wavelength_index(
        hits::AbstractVector{<:Union{GaussianBeamletHit, AstigmaticGaussianBeamletHit}})
    λ = wavelength(first(hits))
    n_medium = refractive_index(first(hits))
    # Relative tolerance of the shared λ and n
    tol = sqrt(eps(float(typeof(λ))))
    for (i, h) in enumerate(hits)
        isapprox(wavelength(h), λ; rtol = tol) ||
            throw(ArgumentError("PlaneField: all hits must share the wavelength, got $(wavelength(h)) m for hit $i and $λ m for hit 1"))
        isapprox(refractive_index(h), n_medium; rtol = tol) ||
            throw(ArgumentError("PlaneField: all hits must share the refractive index, got $(refractive_index(h)) for hit $i and $n_medium for hit 1"))
    end
    return λ, n_medium
end

# Field direction of a hit on the plane with the axes `u` and `v`
function hit_field_direction(hit::GaussianBeamletHit, u::SVector{3, T}, v) where {T}
    d = direction(hit)
    # Degeneracy guard
    e = isparallel3d(d, u; atol = 1e-6) ? v : u
    p = e - dot(d, e) * d
    return SVector{3, Complex{T}}(p / norm(p))
end
hit_field_direction(::AstigmaticGaussianBeamletHit, u, v) = nothing

# 3D field of a hit at `p`, with its direction `pol` from `hit_field_direction`
hit_field(hit::GaussianBeamletHit, p, pol::SVector{3}) = electric_field(hit, p) .* pol
hit_field(hit::AstigmaticGaussianBeamletHit, p, ::Nothing) = electric_field(hit, p)

# -----------------------------------------------------------------------------------------
# OpticsBase.PlaneField -> BeamletOptics.WavefrontBeamletDecomposition

# Forward part of `f`, which must lie in vacuum (n = 1)
function forward_field(f::OpticsBase.PlaneField)
    f.n == 1 || throw(ArgumentError(
        "PlaneField: the field must lie in vacuum or air (n = 1) to be decomposed, got n = $(f.n)"))
    E = forward(f).E
    norm(E) > sqrt(eps(real(eltype(E)))) * norm(f.E) || throw(ArgumentError(
        "PlaneField: the field has no forward-travelling light to decompose"))
    return E
end

# Plane axes `(u, v, n)` of `f`
plane_frame(f::OpticsBase.PlaneField) = (SVector{3}(f.axes[:, 1]), SVector{3}(f.axes[:, 2]),
    SVector{3}(f.axes[:, 3]))

"""
    BeamletOptics.WavefrontBeamletDecomposition(f::OpticsBase.PlaneField; kwargs...)

Decomposes the forward-travelling part of a [`PlaneField`](https://StackEnjoyer.github.io/OpticsBase.jl/dev/)
`f` into an [`AstigmaticBeamGroup`](@ref BeamletOptics.AstigmaticBeamGroup) of
`AstigmaticGaussianBeamlet`s, placed on `f`'s plane and travelling along `f`'s normal
`n`. This is the return path for `OpticsBase.PlaneField(detector; ...)`.

# Arguments

- `f`: the field to decompose. Only `forward(f)` (the part travelling along `+n`) is
  used; the backward part is dropped.
- `kwargs...`: forwarded to [`WavefrontBeamletDecomposition`](@ref
  BeamletOptics.WavefrontBeamletDecomposition)`(x, y, Eu, Ev, pos, dir, λ; kwargs...)`
  (`threshold`, `overlap`, `randomize_axes`, `rng`); `basis` is fixed to `f`'s own `(u,
  v)` and `R` to `f.R`, neither can be overridden.

# Conventions

- The physical field used is `forward(f).E .* reference_phase(f)`, in \\[V/m\\] peak
  amplitude, matching the units `AstigmaticGaussianBeamlet`'s `E0` expects. The samples and
  the reference sphere are handed on separately (`R = f.R`).
- The beam group is placed at `f.origin` (`pos = f.origin`) with the local axes
  `f.axes[:, 1:2]`; its sampling grid is `f`'s own `(u, v)` grid (`OpticsBase.coordinates`).
- `f` must lie in vacuum or air (`f.n == 1`).

!!! warning
    `forward` splits each sample as a local plane wave along the reference direction (`n` for
    `R = Inf`); a beam at the angle `θ` to it loses an amplitude of the order
    `(1 − cos θ)/2`. Keep the plane close to normal to the beam.

Throws an `ArgumentError` if `f.n != 1`, if `f` has no forward-travelling light (e.g. a
purely backward wave), or if `kwargs` contains `basis` or `R`, and whatever
`WavefrontBeamletDecomposition(x, y, Eu, Ev, pos, dir, λ)` throws.
"""
function BeamletOptics.WavefrontBeamletDecomposition(f::OpticsBase.PlaneField; kwargs...)
    # Keywords fixed by the PlaneField
    for key in (:basis, :R)
        haskey(kwargs, key) && throw(ArgumentError(
            "WavefrontBeamletDecomposition: `$key` is taken from the PlaneField and cannot be overridden"))
    end
    E = forward_field(f)
    u, v, n = plane_frame(f)
    return WavefrontBeamletDecomposition(collect(OpticsBase.coordinates(f, 1)),
        collect(OpticsBase.coordinates(f, 2)), E[:, :, 1], E[:, :, 2], SVector{3}(f.origin), n,
        f.λ; basis = (u, v), R = f.R, kwargs...)
end

# -----------------------------------------------------------------------------------------
# OpticsBase.PlaneField -> BeamletOptics.GaussianModeDecomposition

"""
    BeamletOptics.GaussianModeDecomposition(f::OpticsBase.PlaneField)

Fits a single `AstigmaticGaussianBeamlet` to the forward-travelling part of a
[`PlaneField`](https://StackEnjoyer.github.io/OpticsBase.jl/dev/), e.g. the output of a single-mode fiber. The
beamlet starts on `f`'s plane at the centroid of the field; see
[`GaussianModeDecomposition`](@ref BeamletOptics.GaussianModeDecomposition) for the fit.
Use it for fields close to one Gaussian mode that are too small for
[`WavefrontBeamletDecomposition`](@ref BeamletOptics.WavefrontBeamletDecomposition)
(which needs beamlets much larger than λ).

`optical_power(beamlet) / OpticsBase.power(f)` is the fraction of the power captured by
the Gaussian mode. The field used is `forward(f).E .* reference_phase(f)`, with the
reference sphere passed on as `R = f.R`. `f` must lie in vacuum or air (`f.n == 1`) and must
contain forward-travelling light; an `ArgumentError` is thrown otherwise.
"""
function BeamletOptics.GaussianModeDecomposition(f::OpticsBase.PlaneField)
    E = forward_field(f)
    u, v, n = plane_frame(f)
    return GaussianModeDecomposition(collect(OpticsBase.coordinates(f, 1)),
        collect(OpticsBase.coordinates(f, 2)), E[:, :, 1], E[:, :, 2], SVector{3}(f.origin), n,
        f.λ; basis = (u, v), R = f.R)
end

end # module
