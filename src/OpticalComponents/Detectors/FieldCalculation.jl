"""
    electric_field(pd::Detector, s::AbstractDetectorSampling = RegularGrid(); progress = true)

Samples the electric field of the hits of the [`Detector`](@ref) `pd` at the points of the
sampling `s` in its local (x, z) plane, by default a [`RegularGrid`](@ref) in a window around
the hits. The field at the local point `(x, z)` is the coherent sum of the fields `Eⱼ` of all
hits `j`,

```math
\\mathbf{E}(x, z) = \\sum_j \\mathbf{E}_j(\\mathbf{p}), \\qquad \\mathbf{p} = \\mathbf{p}_d + x\\,\\hat{\\mathbf{x}} + z\\,\\hat{\\mathbf{z}},
```

with the position `p_d` of the detector and its local axes `x̂`, `ẑ`, which form a left-handed
system with the detector normal. `Eⱼ` is defined per hit type in
[`electric_field(::BeamletOptics.AbstractDetectorHit, ::AbstractVector)`](@ref). The field is
not projected onto the detector: it is the field that is there, for any angle of incidence.
The power per detector area follows from [`intensity`](@ref), which takes the angle of incidence
of each hit into account. To sample on a grid of one's own, pass a [`RectilinearGrid`](@ref);
other samplings are subtypes of [`BeamletOptics.AbstractDetectorSampling`](@ref).

!!! warning "Coherent sum"
    All hits are summed as one monochromatic, mutually coherent field, also hits of different
    wavelengths or of sources that are mutually incoherent. For those, evaluate the
    [`intensity`](@ref) of each wavelength or source on its own and add the intensities.

!!! note "Resetting detectors"
    Be sure to call `empty!(pd)` (or [`initialize!`](@ref)`(system)` for all detectors) before each new measurement if reusing the same detector.

# Keyword arguments

- `progress::Bool = true`: shows a progress bar once the calculation has run for
  `get_progress_threshold()` seconds (default 5 s). It is only drawn if `stderr` is a terminal.

The grid, its window and its resolution are set by the sampling, e.g.
`electric_field(pd, RegularGrid(n = 200, crop_factor = 5))`, see [`RegularGrid`](@ref).

!!! note "Scaling"
    The returned values for ray hits correspond to the E-field of the point spread function.
    The values are raw/unscaled and not equal to a Strehl ratio: a [`Ray`](@ref) has the
    amplitude `1`, a [`PolarizedRay`](@ref) its field vector `E0`, which carries the
    Fresnel/Jones amplitude factors but not the ray-tube area or the pupil sampling density.

# Returns

The coordinates of the sampling followed by the field at its points, for the grids
`(xs, zs, E)` with `E[i, j]` at `(xs[i], zs[j])`, see [`sample_points`](@ref). The values are
complex scalars for [`Ray`](@ref)s and [`GaussianBeamlet`](@ref)s, `Point3{Complex{T}}` 3D vectors
in global coordinates for [`PolarizedRay`](@ref)s and [`AstigmaticGaussianBeamlet`](@ref)s,
which carry a polarization.

!!! note "Unpolarized light"
    An unpolarized PSF is **not** obtained from a single coherent trace. Instead, trace twice
    with orthogonal input polarizations (e.g. `E0 = [1,0,0]` and `E0 = [0,0,1]`), evaluate
    [`intensity`](@ref)`(pd, s)` with the same [`RectilinearGrid`](@ref) `s` after each trace
    and add the results.
"""
function electric_field(pd::Detector, s::AbstractDetectorSampling = RegularGrid();
        progress::Bool = true)
    grid = sample_points(pd, s)
    return (grid.coordinates..., electric_field(pd, grid.points; progress))
end

"""
    electric_field(pd::Detector, points::AbstractArray{<:AbstractVector}; progress = true)

The field of the hits of the [`Detector`](@ref) `pd` at the local detector coordinates
`points` in \\[m\\], an array of vectors `(x, z)` (e.g. `Point2`s), see
[`electric_field`](@ref)`(pd, s)`. Returns the values only, in the shape of `points`.
"""
electric_field(pd::Detector, points::AbstractArray{<:AbstractVector}; progress::Bool = true) =
    _sample(electric_field, pd, points, progress)

"""
    electric_field(hits::AbstractVector{<:AbstractDetectorHit}, p::AbstractVector)

Coherent sum of the fields of the detector `hits` at the global point `p` \\[m\\], see
[`electric_field(::BeamletOptics.AbstractDetectorHit, ::AbstractVector)`](@ref).
"""
electric_field(hits::AbstractVector{<:AbstractDetectorHit}, p::AbstractVector) =
    sum(hit -> electric_field(hit, p), hits)

"""
    electric_field(hit::AbstractDetectorHit, p::AbstractVector)

Electric field `Eⱼ(p)` of a single detector `hit` at the global point `p` \\[m\\], not projected
onto the detector: a complex scalar for [`RayHit`](@ref)s and [`GaussianBeamletHit`](@ref)s, a
complex 3D vector in global coordinates for [`PolarizedRayHit`](@ref)s and
[`AstigmaticGaussianBeamletHit`](@ref)s. It is a phasor with the time dependence `exp(-iωt)`, so
that a plane wave in a medium of index `n` carries the intensity `n|E|²/(2Z₀)`. With the vacuum
wavenumber `k₀ = 2π/λ`, the direction `s` of the hit, the refractive index `n` in front of the
detector and the optical path length `L` from the source to the hit point `p_hit`:

- [`RayHit`](@ref): a plane wave of amplitude `1`,
  ```math
  E_j(\\mathbf{p}) = \\exp\\!\\big(i k_0 \\,[L + n\\, (\\mathbf{p} - \\mathbf{p}_\\mathrm{hit}) \\cdot \\mathbf{s}]\\big)
  ```
- [`PolarizedRayHit`](@ref): the same plane wave times the field vector `E₀` of the ray,
  ``\\mathbf{E}_j(\\mathbf{p}) = \\mathbf{E}_0 \\exp(i k_0 [L + n\\, (\\mathbf{p} - \\mathbf{p}_\\mathrm{hit}) \\cdot \\mathbf{s}])``
- [`GaussianBeamletHit`](@ref): the TEM₀₀ field of the traced segment, at the distance `l` along
  its chief ray and the distance `r` from it,
  ```math
  E_j(\\mathbf{p}) = \\frac{E_0 w_0}{\\sqrt{n}\\, w(l)} \\exp\\!\\left(-\\frac{r^2}{w(l)^2}\\right)
      \\exp\\!\\left(i\\left[k_0 L(l) + \\psi(l) + \\frac{n k_0 r^2}{2 R(l)}\\right]\\right)
  ```
  with the field `E₀` and waist `w₀` of the beamlet, the beam radius `w`, the radius of curvature
  `R` and the Gouy phase `ψ` of the segment (see [`gauss_parameters`](@ref)), the index `n` of
  the segment and the optical path `L(l)` along the chief ray.
- [`AstigmaticGaussianBeamletHit`](@ref): the field vector of the chief ray of the segment, `E_pol`,
  times the general astigmatic Gaussian beam of its parabasal rays,
  ``\\mathbf{E}_j(\\mathbf{p}) = \\mathbf{E}_\\mathrm{pol}\\, g(l) \\exp(i k_0 [L(l) + n\\, W(\\mathbf{p})])``,
  with the Gouy factor `g`, which also holds the change of the amplitude with the beam area, and
  the complex transverse phase `W`, whose imaginary part is the Gaussian envelope, see
  [Astigmatic polarized beamlets](@ref).

Beamlet hits give the physical field in \\[V/m\\]. Ray hits are raw and unscaled: `|E|²` does not
carry the power of the ray, see the note on scaling in [`electric_field`](@ref)`(pd, s)`.

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

# Phase factor of a ray hit at `p`, from its optical path and the optical distance along the ray
_ray_phasor(hit::AbstractRayHit, p) = cis(wavenumber(hit) *
    (optical_path_length(hit) + refractive_index(hit) * dot(p - hit_point(hit), direction(hit))))

"""
    intensity(pd::Detector, s::AbstractDetectorSampling = RegularGrid(); progress = true)

Samples the irradiance on the [`Detector`](@ref) `pd`, the power per detector area in \\[W/m²\\],
at the points of the sampling `s`, by default a [`RegularGrid`](@ref) around the hits. Returns
the coordinates of the sampling followed by the irradiance at its points, for the grids
`(xs, zs, I)`; `progress` as for [`electric_field`](@ref)`(pd, s)`.

`I` is the flux of the coherent sum of all hits through the detector: the time-averaged Poynting
vector along the detector normal `n̂`, which points along the incoming beams. Hits with a
polarization ([`PolarizedRay`](@ref)s, [`AstigmaticGaussianBeamlet`](@ref)s) give

```math
I = \\frac{1}{2} \\operatorname{Re}\\big[(\\mathbf{E} \\times \\mathbf{H}^*) \\cdot \\hat{\\mathbf{n}}\\big],
\\qquad \\mathbf{E} = \\sum_j \\mathbf{E}_j, \\qquad
\\mathbf{H} = \\sum_j \\frac{n_j}{Z_0}\\, \\sigma_j \\mathbf{s}_j \\times \\mathbf{E}_j,
```

with the field `Eⱼ` of hit `j` (see [`electric_field`](@ref)`(hit, p)`), its direction `sⱼ`, the
refractive index `nⱼ` in front of the detector, the vacuum impedance `Z₀` and
`σⱼ = sign(sⱼ ⋅ n̂)`, i.e. a hit on the back side of the detector counts like one on its front
side. Hits without a polarization ([`Ray`](@ref)s, [`GaussianBeamlet`](@ref)s) have no vector `H`
and give the scalar form of the same flux,

```math
I = \\frac{1}{2} \\operatorname{Re}\\big[E^* F\\big], \\qquad E = \\sum_j E_j, \\qquad
F = \\sum_j \\frac{n_j}{Z_0} \\cos\\theta_j\\, E_j, \\qquad \\cos\\theta_j = |\\mathbf{s}_j \\cdot \\hat{\\mathbf{n}}|.
```

A single hit under the angle `θ` to the normal gives `I = n|E|² cos θ / (2Z₀)`, the power of the
beam per detector area, so that [`optical_power`](@ref) is the power of the beam for any angle of
incidence. No projection factor is applied to the field. For two scalar hits, the cross term is
`Re(E₁* E₂) (n₁ cos θ₁ + n₂ cos θ₂) / (2Z₀)`, i.e. hits under different angles interfere with the
mean of their `n cos θ`.

`I` is in \\[W/m²\\] for beamlet hits. For ray hits, whose fields are unscaled, it has the same
form but no absolute unit.
"""
function intensity(pd::Detector, s::AbstractDetectorSampling = RegularGrid();
        progress::Bool = true)
    grid = sample_points(pd, s)
    return (grid.coordinates..., intensity(pd, grid.points; progress))
end

"""
    intensity(pd::Detector, points::AbstractArray{<:AbstractVector}; progress = true)

The irradiance on the [`Detector`](@ref) `pd` in \\[W/m²\\] at the local detector coordinates
`points` in \\[m\\], an array of vectors `(x, z)` (e.g. `Point2`s), see
[`intensity`](@ref)`(pd, s)`. Returns the values only, in the shape of `points`.
"""
function intensity(pd::Detector, points::AbstractArray{<:AbstractVector};
        progress::Bool = true)
    # downstream normal: the detector normal points against the incoming beams
    n̂ = Point3(-orientation(pd)[:, 2])
    return _sample((hs, p) -> intensity(hs, p, n̂), pd, points, progress)
end

"""
    intensity(hits::AbstractVector{<:AbstractDetectorHit}, p::AbstractVector, n̂::AbstractVector)

Irradiance \\[W/m²\\] of the coherent sum of the detector `hits` at the global point `p` through
a plane with the unit normal `n̂`, see [`intensity`](@ref)`(pd, s)`.
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
    optical_power(pd::Detector, s::AbstractDetectorSampling = RegularGrid(); progress = true)

Calculates the total optical power on `pd` in \\[W\\] as the integral of its [`intensity`](@ref)
over the window of the sampling `s` in the local detector coordinates,

```math
P = \\iint I(x, z)\\, \\mathrm{d}x\\, \\mathrm{d}z \\approx \\sum_i w_i\\, I(x_i, z_i),
```

with the points `(xᵢ, zᵢ)` and the weights `wᵢ` \\[m²\\] of the sampling, see
[`sample_points`](@ref): the trapezoidal rule for [`RegularGrid`](@ref) and
[`RectilinearGrid`](@ref). Since `I` is the flux through the detector, `P` is the power of the
beams for any angle of incidence.

The integral runs over the window of the sampling, not over the active area of the detector:
the window is not cut at the edges of the detector. The default window of beamlet hits holds
the power of each beamlet except for a part of `10⁻⁸`, see `crop_factor` of
[`RegularGrid`](@ref). For ray hits, `P` has no absolute unit, see [`intensity`](@ref).
"""
function optical_power(pd::Detector, s::AbstractDetectorSampling = RegularGrid();
        progress::Bool = true)
    grid = sample_points(pd, s)
    I = intensity(pd, grid.points; progress)
    # lazy: neither the weights nor their product with I are materialized
    return sum(Broadcast.instantiate(Broadcast.broadcasted(*, grid.weights, I)))
end

# Samples `f(hits, p)` at the global points `p` of the local detector coordinates `points`,
# threaded over chunks of the points. All readouts of a detector share this loop.
function _sample(f, pd::Detector, points::AbstractArray, progress::Bool)
    # PD local coordinate axes (left-handed coords.), flipped x due to the rotated pd mesh
    local_x = Point3(-orientation(pd)[:, 1])
    local_z = Point3(orientation(pd)[:, 3])
    return _sample(f, _detector_hits(pd), position(pd), local_x, local_z, points, progress)
end

# function barrier: `hits` has a concrete type here
function _sample(f, hits::AbstractVector, origin, local_x, local_z, points::AbstractArray,
        progress::Bool)
    point(q) = origin + q[1] * local_x + q[2] * local_z
    values = similar(points, typeof(f(hits, point(first(points)))))
    # equal chunks of consecutive points, independent of the shape; one tick per chunk
    nparts = min(length(points), max(100, 4 * Threads.nthreads()))
    parts = collect(Iterators.partition(eachindex(values, points), cld(length(points), nparts)))
    _with_progress(progress, length(parts), "Detector field: ") do prog
        Threads.@threads for part in parts
            @inbounds for k in part
                values[k] = f(hits, point(points[k]))
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
