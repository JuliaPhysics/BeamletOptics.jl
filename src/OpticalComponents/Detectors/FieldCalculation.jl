"""
    electric_field(pd::Detector, sampling::AbstractDetectorSampling = RegularGrid(); progress = true)

Samples the electric field of the hits of the [`Detector`](@ref) `pd` at the points of
`sampling` in its local (x, z) plane, by default a [`RegularGrid`](@ref) in a window around
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
`electric_field(pd, RegularGrid(num_points = 200, crop_factor = 5))`, see [`RegularGrid`](@ref).

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
    [`intensity`](@ref)`(pd, sampling)` with the same [`RectilinearGrid`](@ref) `sampling` after
    each trace and add the results.
"""
function electric_field(pd::Detector, sampling::AbstractDetectorSampling = RegularGrid();
        progress::Bool = true)
    grid = sample_points(pd, sampling)
    return (grid.coordinates..., electric_field(pd, grid.points; progress))
end

"""
    electric_field(pd::Detector, points::AbstractArray{<:AbstractVector}; progress = true)

The field of the hits of the [`Detector`](@ref) `pd` at the local detector coordinates
`points` in \\[m\\], an array of vectors `(x, z)` (e.g. `Point2`s), see
[`electric_field`](@ref)`(pd, sampling)`. Returns the values only, in the shape of `points`.
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
- [`AstigmaticGaussianBeamletHit`](@ref): the field vector of the chief ray of the segment, `E₀`,
  times the general astigmatic Gaussian beam of its parabasal rays,
  ``\\mathbf{E}_j(\\mathbf{p}) = \\mathbf{E}_0\\, g(l) \\kappa^{-1/4} \\exp(i k_0 [L(l) + n\\, W(\\mathbf{p})]) \\exp(-\\mathbf{r}^T \\mathbf{S}^{-1} \\mathbf{r})``,
  with the Gouy factor `g`, which also holds the change of the amplitude with the beam area, the
  transverse phase `W`, the spot matrix `S` at the offset `r` from the chief ray and
  `κ = det S / |a|²` of the complex beam area `a`, see [Astigmatic polarized beamlets](@ref).

Beamlet hits give the physical field in \\[V/m\\]. Ray hits are raw and unscaled: `|E|²` does not
carry the power of the ray, see the note on scaling in [`electric_field`](@ref)`(pd, sampling)`.

A new hit type implements this method and the other functions listed in
[`BeamletOptics.AbstractDetectorHit`](@ref); the readouts of a [`Detector`](@ref) on a grid
([`electric_field`](@ref), [`intensity`](@ref), [`optical_power`](@ref)) then work for it.
"""
electric_field(hit::AbstractDetectorHit, ::AbstractVector) =
    error(lazy"electric_field not defined for $(typeof(hit))")

electric_field(hit::RayHit, p::AbstractVector) = ray_phasor(hit, p)

electric_field(hit::PolarizedRayHit, p::AbstractVector) = polarization(hit) * ray_phasor(hit, p)

function electric_field(hit::GaussianBeamletHit, p::AbstractVector)
    v = p - hit.p0
    l1 = _pseudo_dot(v, direction(hit))

    # Transverse distance r
    r = norm(v - l1 * direction(hit))

    # Beam parameters of the segment
    w, curv, ψ, _ = gauss_parameters(hit.segment, l1)
    ψ += hit.ψ0

    # Optical path up to the point, curvature term with n k0
    k0 = wavenumber(wavelength(hit))
    n = refractive_index(hit)
    opl = hit.l0 + l1 + hit.Δl0 + (n - 1) * l1
    # Field in the medium, smaller by √n
    return hit.E_amp / (sqrt(n) * w) * exp(-r^2 / w^2) * cis(k0 * opl + ψ + n * k0 * r^2 * curv / 2)
end

electric_field(hit::AstigmaticGaussianBeamletHit, p::AbstractVector) =
    hit.E0 * astigmatic_phasor(hit, p)

# Phase factor of a ray hit at `p`, from its optical path and the optical distance along the ray
ray_phasor(hit::AbstractRayHit, p) = cis(wavenumber(hit) *
    (optical_path_length(hit) + refractive_index(hit) * dot(p - hit_point(hit), direction(hit))))

"""
    intensity(pd::Detector, sampling::AbstractDetectorSampling = RegularGrid(); progress = true)

Samples the irradiance on the [`Detector`](@ref) `pd`, the power per detector area in \\[W/m²\\],
at the points of `sampling`, by default a [`RegularGrid`](@ref) around the hits. Returns
the coordinates of the sampling followed by the irradiance at its points, for the grids
`(xs, zs, I)`; `progress` as for [`electric_field`](@ref)`(pd, sampling)`.

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
function intensity(pd::Detector, sampling::AbstractDetectorSampling = RegularGrid();
        progress::Bool = true)
    grid = sample_points(pd, sampling)
    return (grid.coordinates..., intensity(pd, grid.points; progress))
end

"""
    intensity(pd::Detector, points::AbstractArray{<:AbstractVector}; progress = true)

The irradiance on the [`Detector`](@ref) `pd` in \\[W/m²\\] at the local detector coordinates
`points` in \\[m\\], an array of vectors `(x, z)` (e.g. `Point2`s), see
[`intensity`](@ref)`(pd, sampling)`. Returns the values only, in the shape of `points`.
"""
function intensity(pd::Detector, points::AbstractArray{<:AbstractVector};
        progress::Bool = true)
    # Downstream normal
    n̂ = -direction(pd)
    return _sample((hs, p) -> intensity(hs, p, n̂), pd, points, progress)
end

"""
    intensity(hits::AbstractVector{<:AbstractDetectorHit}, p::AbstractVector, n̂::AbstractVector)

Irradiance \\[W/m²\\] of the coherent sum of the detector `hits` at the global point `p` through
a plane with the unit normal `n̂`, see [`intensity`](@ref)`(pd, sampling)`.
"""
intensity(hits::AbstractVector{<:AbstractDetectorHit}, ::AbstractVector, ::AbstractVector) =
    error(lazy"intensity not defined for $(eltype(hits))")

function intensity(hits::AbstractVector{<:Union{PolarizedRayHit, AstigmaticGaussianBeamletHit}},
        p::AbstractVector, n̂::AbstractVector)
    E = H = zero(Point3{Complex{eltype(p)}})
    for hit in hits
        Ej = electric_field(hit, p)
        # Back side counts like front side
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
    optical_power(pd::Detector, sampling::AbstractDetectorSampling = RegularGrid(); progress = true)

Calculates the total optical power on `pd` in \\[W\\] as the integral of its [`intensity`](@ref)
over the window of `sampling` in the local detector coordinates,

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
function optical_power(pd::Detector, sampling::AbstractDetectorSampling = RegularGrid();
        progress::Bool = true)
    grid = sample_points(pd, sampling)
    I = intensity(pd, grid.points; progress)
    # Lazy weighted sum
    return sum(Broadcast.instantiate(Broadcast.broadcasted(*, grid.weights, I)))
end

# Samples `f(hits, p)` at the global points `p` of the local detector coordinates `points`,
# threaded over chunks of the points. All readouts of a detector share this loop.
function _sample(f, pd::Detector, points::AbstractArray, progress::Bool)
    local_x, local_z = local_axes(pd)
    return _sample(f, _detector_hits(pd), position(pd), local_x, local_z, points, progress)
end

# Function barrier
function _sample(f, hits::AbstractVector, origin, local_x, local_z, points::AbstractArray,
        progress::Bool)
    point(q) = origin + q[1] * local_x + q[2] * local_z
    values = similar(points, typeof(f(hits, point(first(points)))))
    # Equal chunks of consecutive points, one tick per chunk
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

# Field of an astigmatic beamlet hit at `p` per unit amplitude: Gouy factor and phase
function astigmatic_phasor(hit::AstigmaticGaussianBeamletHit{G}, p::AbstractArray) where {G}
    dir = direction(hit)
    v = p - hit.p0
    l1 = _pseudo_dot(v, dir)
    r_vec = v - l1 * dir

    h1_z = hit.h1 + l1 * hit.u1
    h2_z = hit.h2 + l1 * hit.u2
    area_z = _pseudo_cross2d(h1_z, h2_z, dir)
    # Regularization of the area at caustics
    if abs(area_z) < 1e-25
        area_z = Complex{G}(1e-25, 1e-25)
    end
    # Gouy factor from the continuous area argument
    gouy = hit.gouy0 / (sqrt(1 - hit.ρ1 * l1) * sqrt(1 - hit.ρ2 * l1))

    n = refractive_index(hit)
    phase_corr = (n - 1) * l1
    z_total = hit.l0 + l1
    # Transverse term with n k0
    return gouy * cis(hit.k0 * (z_total + hit.Δl0 + phase_corr)) *
           transverse_profile(h1_z, hit.u1, h2_z, hit.u2, r_vec, dir, area_z, n * hit.k0)
end
