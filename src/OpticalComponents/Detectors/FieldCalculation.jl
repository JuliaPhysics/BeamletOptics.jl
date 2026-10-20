"""
    electric_field(pd::Detector; kwargs...)

Samples the electric field of the hits of the [`Detector`](@ref) `pd` on a regular `n×n` grid in
its local (x, z) plane, in a window around the hits. The field at a point is the coherent sum of
the fields of all hits, see [`electric_field(::BeamletOptics.AbstractDetectorHit, ::AbstractVector)`](@ref).
It is not projected onto the detector: the power per detector area follows from
[`intensity`](@ref), which takes the angle of incidence of each hit into account.
Note that the `pd` local coordinates are given in a (x, z) basis where the normal vector forms a
left-handed system. To sample on a grid of one's own, use
[`electric_field(::Detector, ::AbstractVector, ::AbstractVector)`](@ref).

!!! note "Resetting detectors"
    Be sure to call `empty!(pd)` (or [`initialize!`](@ref)`(system)` for all detectors) before each new measurement if reusing the same detector.

# Keyword Arguments

The following generic kwargs can be used for all hit types:

- `n::Int=100`
  Number of sample points per axis.
- `crop_factor::Real`
  Scales the width of the sampling window returned by
  [`calc_local_lims`](@ref); values >1 expand, <1 shrink. The default is `1` for ray hits
  (the box around the hit points) and `3` for beamlet hits (three beam radii around each
  hit, which holds the power of a Gaussian beam except for a part of `10⁻⁸`).
- `x_min, x_max, z_min, z_max`
  Manually override the sampling bounds in the local x or z directions.
  If left as `Inf`, the bounds from `calc_local_lims` are used.
- `x0_shift::Real=0, z0_shift::Real=0`
  Applies a constant offset to the entire x or z coordinate arrays,
  useful for recentring or testing alignment.
- `progress::Bool=true`
  Shows a progress bar once the calculation has run for `get_progress_threshold()`
  seconds (default 5 s). It is only drawn if `stderr` is a terminal.

## Ray specific keyword arguments

- `center::AbstractCenterAlgorithm=Centroid()`
  How the sampling window is centred.
  `Centroid()` uses the projection‑weighted centroid,
  `MinMax()` uses the geometric mid‑point of the bounding box.

!!! note "Scaling"
    The returned values for ray hits correspond to the E-field of the point spread function.
    The values are raw/unscaled and not equal to a Strehl ratio: a [`Ray`](@ref) has the
    amplitude `1`, a [`PolarizedRay`](@ref) its field vector `E0`, which carries the
    Fresnel/Jones amplitude factors but not the ray-tube area or the pupil sampling density.

# Returns

A tuple `(xs, zs, E)` where
- `xs::LinRange{T}` and `zs::LinRange{T}` are the sampled coordinates
  in the detector's local x and z axes,
- `E` is the field at `(xs[i], zs[j])`: a `Matrix{Complex{T}}` for [`Ray`](@ref)s and
  [`GaussianBeamlet`](@ref)s, a `Matrix{Point3{Complex{T}}}` of 3D vectors in global coordinates
  for [`PolarizedRay`](@ref)s and [`AstigmaticGaussianBeamlet`](@ref)s, which carry a polarization.

!!! note "Unpolarized light"
    An unpolarized PSF is **not** obtained from a single coherent trace. Instead, trace twice
    with orthogonal input polarizations (e.g. `E0 = [1,0,0]` and `E0 = [0,0,1]`), evaluate
    [`intensity`](@ref)`(pd, xs, zs)` on the same grid after each trace and add the results.
"""
function electric_field(pd::Detector; progress::Bool = true, kwargs...)
    xs, zs = _sampling_grid(pd; kwargs...)
    return xs, zs, electric_field(pd, xs, zs; progress)
end

"""
    electric_field(pd::Detector, xs::AbstractVector, zs::AbstractVector; progress = true)

Samples the electric field of the hits of the [`Detector`](@ref) `pd` on the grid `xs × zs` of
its local coordinates in \\[m\\]: `E[i, j]` is the field at `(xs[i], zs[j])`. The type of `E` and
the meaning of its values are those of [`electric_field`](@ref)`(pd; kwargs...)`, which picks
the grid automatically.
"""
electric_field(pd::Detector, xs::AbstractVector, zs::AbstractVector; progress::Bool = true) =
    _sample(electric_field, pd, xs, zs, progress)

"""
    electric_field(hits::AbstractVector{<:AbstractDetectorHit}, p::AbstractVector)

Coherent sum of the fields of the detector `hits` at the global point `p` \\[m\\], see
[`electric_field(::BeamletOptics.AbstractDetectorHit, ::AbstractVector)`](@ref).
"""
electric_field(hits::AbstractVector{<:AbstractDetectorHit}, p::AbstractVector) =
    sum(hit -> electric_field(hit, p), hits)

"""
    electric_field(hit::AbstractDetectorHit, p::AbstractVector)

Electric field of a single detector `hit` at the global point `p` \\[m\\], not projected onto
the detector: a complex scalar for [`RayHit`](@ref)s and [`GaussianBeamletHit`](@ref)s, a complex
3D vector in global coordinates for [`PolarizedRayHit`](@ref)s and
[`AstigmaticGaussianBeamletHit`](@ref)s. Beamlet hits give the physical field in \\[V/m\\]. A
`RayHit` has the amplitude `1` and a `PolarizedRayHit` the field vector of its ray, see the
note on scaling in [`electric_field`](@ref)`(pd; kwargs...)`.

A new hit type implements this method, `direction(hit)` and `refractive_index(hit)`; the
readouts of a [`Detector`](@ref) on a grid ([`electric_field`](@ref), [`intensity`](@ref),
[`optical_power`](@ref)) then work for it.
"""
electric_field(::AbstractDetectorHit, ::AbstractVector)

electric_field(hit::RayHit, p::AbstractVector) = _ray_phasor(hit, p)

electric_field(hit::PolarizedRayHit, p::AbstractVector) = polarization(hit) * _ray_phasor(hit, p)

# The scalar beamlet model keeps |E|²/(2Z₀) at the beam power in any medium, i.e. its amplitude
# is √n larger than the physical field
electric_field(hit::GaussianBeamletHit, p::AbstractVector) =
    beamlet_hit_field(hit, p) / sqrt(refractive_index(hit))

electric_field(hit::AstigmaticGaussianBeamletHit, p::AbstractVector) =
    hit.E_pol * _astigmatic_phasor(hit, p)

# Phase factor of a ray hit at `p`, from its optical path and the distance along the ray
_ray_phasor(hit::AbstractRayHit, p) =
    cis(wavenumber(hit) * (optical_path_length(hit) + dot(p - hit_point(hit), direction(hit))))

"""
    intensity(pd::Detector; kwargs...)

Samples the irradiance on the [`Detector`](@ref) `pd`, the power per detector area in \\[W/m²\\],
on the grid of [`electric_field`](@ref)`(pd; kwargs...)`, with the same keyword arguments.
Returns a tuple `(xs, zs, I)`. To sample on a grid of one's own, use
[`intensity(::Detector, ::AbstractVector, ::AbstractVector)`](@ref).

`I` is the flux of the coherent sum of all hits through the detector, the Poynting vector along
its normal `n̂`, `I = ½ Re(E × H*) ⋅ n̂`. Each hit adds `(n/Z₀) s × E` to `H`, with its direction
`s`, the refractive index `n` in front of the detector and the vacuum impedance `Z₀`. A single
hit under the angle `θ` to the normal gives `I = n |E|² cos θ / (2Z₀)`, so that
[`optical_power`](@ref) is the power of the beam for any angle of incidence. Hits under
different angles interfere with the mean of their `cos θ`. Hits without a polarization
([`Ray`](@ref)s, [`GaussianBeamlet`](@ref)s) have no vector `H` and use the scalar form of the
same flux, `½ Re(E* F)` with `F = Σ (n/Z₀) cos θ E` over the hits. Hits on the back side of the
detector count like hits on its front side.
"""
function intensity(pd::Detector; progress::Bool = true, kwargs...)
    xs, zs = _sampling_grid(pd; kwargs...)
    return xs, zs, intensity(pd, xs, zs; progress)
end

"""
    intensity(pd::Detector, xs::AbstractVector, zs::AbstractVector; progress = true)

Samples the irradiance on the [`Detector`](@ref) `pd` in \\[W/m²\\] on the grid `xs × zs` of its
local coordinates in \\[m\\]: `I[i, j]` is the irradiance at `(xs[i], zs[j])`, see
[`intensity`](@ref)`(pd; kwargs...)`.
"""
function intensity(pd::Detector, xs::AbstractVector, zs::AbstractVector; progress::Bool = true)
    # downstream normal: the detector normal points against the incoming beams
    n̂ = Point3(-orientation(pd)[:, 2])
    return _sample((hs, p) -> intensity(hs, p, n̂), pd, xs, zs, progress)
end

"""
    intensity(hits::AbstractVector{<:AbstractDetectorHit}, p::AbstractVector, n̂::AbstractVector)

Irradiance \\[W/m²\\] of the coherent sum of the detector `hits` at the global point `p` through
a plane with the unit normal `n̂`, see [`intensity`](@ref)`(pd; kwargs...)`.
"""
intensity(::AbstractVector{<:AbstractDetectorHit}, ::AbstractVector, ::AbstractVector)

function intensity(hits::AbstractVector{<:Union{PolarizedRayHit, AstigmaticGaussianBeamletHit}},
        p::AbstractVector, n̂::AbstractVector)
    E = H = zero(Point3{Complex{eltype(p)}})
    for hit in hits
        Ej = electric_field(hit, p)
        # a hit on the back side counts like one on the front side
        s = sign(dot(direction(hit), n̂)) * direction(hit)
        E += Ej
        H += refractive_index(hit) * cross(s, Ej)
    end
    return real(dot(n̂, cross(E, conj(H)))) / (2 * Z_vacuum)
end

function intensity(hits::AbstractVector{<:Union{RayHit, GaussianBeamletHit}},
        p::AbstractVector, n̂::AbstractVector)
    E = F = zero(Complex{eltype(p)})
    for hit in hits
        Ej = electric_field(hit, p)
        E += Ej
        F += refractive_index(hit) * abs(dot(direction(hit), n̂)) * Ej
    end
    return real(conj(E) * F) / (2 * Z_vacuum)
end

"""
    optical_power(pd::Detector; kwargs...)

Calculates the total optical power on `pd` in [W] by integration of its [`intensity`](@ref)
over the sampling window. For more information on keyword argument options, refer to the
[`electric_field`](@ref) docs.
"""
function optical_power(pd::Detector; kwargs...)
    x, y, I = intensity(pd; kwargs...)
    return trapz((x, y), I)
end

# The local grid of `electric_field(pd; kwargs...)`: `n` points per axis in the window of
# `calc_local_lims`, or within the given bounds
function _sampling_grid(pd::Detector; n::Int = 100, x_min = Inf, x_max = Inf, z_min = Inf,
        z_max = Inf, x0_shift::Real = 0, z0_shift::Real = 0, kwargs...)
    _x_min, _x_max, _z_min, _z_max = calc_local_lims(pd, _detector_hits(pd); kwargs...)
    if x_min != Inf && x_max != Inf
        _x_min = x_min
        _x_max = x_max
    end
    if z_min != Inf && z_max != Inf
        _z_min = z_min
        _z_max = z_max
    end
    xs = LinRange(_x_min, _x_max, n) .+ x0_shift
    zs = LinRange(_z_min, _z_max, n) .+ z0_shift
    return xs, zs
end

# Samples `f(hits, p)` at the global points `p` of the local detector coordinates `xs × zs`,
# threaded over `zs`. All readouts of a detector on a grid share this loop.
function _sample(f, pd::Detector, xs::AbstractVector, zs::AbstractVector, progress::Bool)
    # PD local coordinate axes (left-handed coords.), flipped x due to the rotated pd mesh
    local_x = Point3(-orientation(pd)[:, 1])
    local_z = Point3(orientation(pd)[:, 3])
    return _sample(f, _detector_hits(pd), position(pd), local_x, local_z, xs, zs, progress)
end

# function barrier: `hits` has a concrete type here
function _sample(f, hits::AbstractVector, origin, local_x, local_z, xs, zs, progress::Bool)
    point(i, j) = origin + xs[i] * local_x + zs[j] * local_z
    values = Matrix{typeof(f(hits, point(1, 1)))}(undef, length(xs), length(zs))
    _with_progress(progress, length(zs), "Detector field: ") do prog
        Threads.@threads for j in eachindex(zs)
            @inbounds for i in eachindex(xs)
                values[i, j] = f(hits, point(i, j))
            end
            _tick!(prog)
        end
    end
    return values
end

function _detector_hits(pd::Detector)
    hs = hits(pd)
    isnothing(hs) && throw(ErrorException("No hits available on detector."))
    return hs
end

"""
    beamlet_hit_field(hit::AbstractBeamletHit, p::AbstractArray)

Complex scalar field of a single beamlet `hit` at the global point `p` \\[m\\]. For an
[`AstigmaticGaussianBeamletHit`](@ref) it is the field along the polarization of the hit, see
[`beamlet_hit_polarization`](@ref) for the 3D field vector; for a [`GaussianBeamletHit`](@ref) it
is the field of the scalar model, whose amplitude is `√n` too large in a medium of index `n`.
The fields that the readouts of a [`Detector`](@ref) add are those of
[`electric_field(::BeamletOptics.AbstractDetectorHit, ::AbstractVector)`](@ref). It is kept as
its own function so that other code that needs the *individual* hit contributions at a point,
e.g. the `BeamletOpticsOpticsBaseExt` package extension, does not have to duplicate the beamlet
field formulas.
"""
function beamlet_hit_field(hit::GaussianBeamletHit{G}, p::AbstractArray) where {G}
    v = p - hit.p0
    l1 = _pseudo_dot(v, hit.d0)

    # Transverse distance r
    r = norm(v - l1 * hit.d0)

    # Beam radius, curvature and Gouy phase from the waist and divergence rays of the segment
    w, curv, ψ, _ = gauss_parameters(hit.segment, l1)
    ψ += hit.ψ0

    # Optical path up to the point, l0 + l1 along the beam. The curvature term is a geometric
    # length like l1, hence it enters with the wavenumber n k0 of the medium as well.
    k0 = wavenumber(hit.segment.λ)
    n = hit.segment.n
    opl = hit.l0 + l1 + hit.Δl0 + (n - 1) * l1
    return hit.E_amp / w * exp(-r^2 / w^2) * cis(k0 * opl + ψ + n * k0 * r^2 * curv / 2)
end

beamlet_hit_field(hit::AstigmaticGaussianBeamletHit, p::AbstractArray) =
    hit.E_ref_amp * _astigmatic_phasor(hit, p)

# Field of an astigmatic beamlet hit at `p` per unit amplitude: Gouy factor and phase
function _astigmatic_phasor(hit::AstigmaticGaussianBeamletHit{G}, p::AbstractArray) where {G}
    v = p - hit.p0
    l1 = _pseudo_dot(v, hit.d0)
    r_vec = v - l1 * hit.d0

    h1_z = hit.h1 + l1 * hit.u1
    h2_z = hit.h2 + l1 * hit.u2
    area_z = _pseudo_cross2d(h1_z, h2_z, hit.d0)
    if abs(area_z) < 1e-25
        area_z = Complex{G}(1e-25, 1e-25)
    end
    # Gouy factor √(a_ref / a(l1)) with the continuous argument of the area: each factor
    # of a(l1) stays off the negative real axis (see `_area_inverse_roots`), so the
    # principal square roots follow the Gouy phase through a focus
    gouy = hit.gouy0 / (sqrt(1 - hit.ρ1 * l1) * sqrt(1 - hit.ρ2 * l1))

    ξ1 = _pseudo_cross2d(h1_z, r_vec, hit.d0)
    ξ2 = _pseudo_cross2d(h2_z, r_vec, hit.d0)
    w = (ξ1 * _pseudo_dot(hit.u2, r_vec) -
         ξ2 * _pseudo_dot(hit.u1, r_vec)) / (2 * area_z)

    phase_corr = (hit.n_eff - 1) * l1
    z_total = hit.l0 + l1
    # the transverse term w enters with the wavenumber n k0 of the medium, like l1
    return gouy * cis(hit.k0 * (z_total + hit.n_eff * w + hit.Δl + phase_corr))
end

"""
    beamlet_hit_polarization(hit::AstigmaticGaussianBeamletHit)

Complex 3D vector (global frame) that turns the scalar
[`beamlet_hit_field`](@ref)`(hit, p)` of an astigmatic beamlet hit into its 3D field
vector at `p`: the polarization of the hit segment's chief ray divided by the complex
amplitude `E_ref_amp` that the scalar field already contains, i.e. a unit vector with the
phase of its largest component removed. It is the zero vector if the segment carries no
light, e.g. behind a crossed polarizer.

Stigmatic [`GaussianBeamletHit`](@ref)s have no method: the underlying [`GaussianBeamlet`](@ref)
model is scalar and carries no polarization in BMO, so there is no principled 3D
direction to return. Callers that need a vector for such a hit (e.g. to build a
tangential field on some plane) must supply their own convention, such as the local `u`
axis of that plane (this matches `OpticsBase`'s rule for a purely scalar field: it is
placed entirely along `u`).
"""
function beamlet_hit_polarization(hit::AstigmaticGaussianBeamletHit)
    # `beamlet_hit_field` already carries the complex amplitude `E_ref_amp` (magnitude of the
    # polarization of the segment and phase of its largest component), so dividing it out
    # keeps amplitude and phase from being counted twice.
    E = hit.E_pol
    iszero(hit.E_ref_amp) && return zero(E)
    return E / hit.E_ref_amp
end
