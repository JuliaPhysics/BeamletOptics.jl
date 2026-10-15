"""
    electric_field(pd::Detector; kwargs...)

Compute a two‐dimensional electric field based on incoming rays or beams as captured by a [`Detector`](@ref).
The returned E-field map is sampled on a regular `n×n` grid in the detector's local (x,z)-plane.
Note that the `pd` local coordinates are given in a (x, z) basis where the normal vector forms a left-handed system.

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
    The values are raw/unscaled and not equal to a Strehl ratio. This feature is not
    yet added. In future versions a pupil finder along with a Strehl estimator will be added.

## Beamlet specific keyword arguments

- `num_spots::Int=50`
  Number of hit spots used to determine bounding box

# Returns

A tuple `(xs, zs, E)` where
- `xs::LinRange{T}` and `zs::LinRange{T}` are the sampled coordinates
  in the detector's local x and z axes,
- `E::Matrix{Complex{T}}` is the corresponding raw/unscaled intensity map,
  except for [`PolarizedRayHit`](@ref)s, see below.

# [`PolarizedRayHit`](@ref) hits

For [`PolarizedRay`](@ref)s, the per-ray `E0` field vectors are added **coherently as 3D vectors**
in global coordinates (not projected onto the detector plane), so `E::Matrix{Point3{Complex{T}}}`.
No obliquity/projection factor is applied: unlike the scalar ray case, the relative projection between
rays is already encoded in their vector directions. As with the scalar case, the result is raw/unscaled
(`E0` carries Fresnel/Jones amplitude factors but not ray-tube area or pupil-sampling density).

!!! note "Unpolarized light"
    Since polarization states with orthogonal `E0` do not interfere, an unpolarized PSF is **not**
    obtained from a single coherent trace. Instead, trace twice with orthogonal input polarizations
    (e.g. `E0 = [1,0,0]` and `E0 = [0,0,1]`) and incoherently add the resulting intensities,
    i.e. `I_total = intensity(E_x) .+ intensity(E_z)`.
"""
electric_field(d::Detector; kwargs...) = electric_field(d, hits(d); kwargs...)

function electric_field(d::Detector, ::Nothing; kwargs...)
    throw(ErrorException("No hits available on detector."))
end

"""
    beamlet_hit_field(hit::AbstractBeamletHit, p::AbstractArray)

Complex scalar field of a single beamlet `hit` at the global point `p` \\[m\\]: the
physical field amplitude in \\[V/m\\], *without* the `√|cos θ|` projection factor that
`electric_field(::Detector, ::Vector{<:AbstractBeamletHit})` multiplies each term with
before summing (that factor makes `|E|²` integrate to the power per detector area; it is
not part of the physical field). It is kept as its own function (rather than inlined in that
loop) so that other code that needs the *individual* hit contributions at a point — e.g.
to build a 3D field vector via [`beamlet_hit_polarization`](@ref), as the
`BeamletOpticsOpticsBaseExt` package extension does to assemble the tangential magnetic
field — does not have to duplicate the beamlet field formulas.
"""
function beamlet_hit_field(hit::GaussianBeamletHit{G}, p::AbstractArray) where {G}
    v = p - hit.p0
    l1 = _pseudo_dot(v, hit.d0)

    # Transverse distance r
    r_vec = v - l1 * hit.d0
    r = norm(r_vec)

    # Distance along beam
    z = hit.l0 + l1

    return electric_field(hit.gauss, r, z; hint = (hit.p0 + l1 * hit.d0, hit.id))
end

function beamlet_hit_field(hit::AstigmaticGaussianBeamletHit{G}, p::AbstractArray) where {G}
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
    ψ = gouy * cis(hit.k0 * (z_total + hit.n_eff * w + hit.Δl + phase_corr))

    return hit.E_ref_amp * ψ
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
    E = polarization(rays(hit.agb.c)[hit.id])
    iszero(hit.E_ref_amp) && return zero(E)
    return E / hit.E_ref_amp
end

function electric_field(
        pd::Detector,
        hits::Vector{GaussianBeamletHit{G}};
        # kwargs
        n::Int = 100,
        crop_factor::Real = 3.0,
        num_spots::Int = 50,
        x_min = Inf,
        x_max = Inf,
        z_min = Inf,
        z_max = Inf,
        x0_shift::Real = 0,
        z0_shift::Real = 0,
        progress::Bool = true,
        kwargs...
) where {G}
    # Calculate autolims
    _x_min, _x_max, _z_min, _z_max = calc_local_lims(pd; crop_factor, num_spots)
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
    # Preallocate e-field matrix
    field = zeros(Complex{G}, n, n)
    # PD local coordinate axes (left-handed coords.)
    local_x = Point3(-orientation(pd)[:, 1])  # flipped sign due to rotated pd mesh
    local_z = Point3(orientation(pd)[:, 3])
    origin = position(pd)
    # Calculate field superposition
    _with_progress(progress, length(zs), "Detector field: ") do prog
        Threads.@threads for j in eachindex(zs)
            z_grid = zs[j]
            # Hoist z-grid math out of inner loop
            p_row = origin + z_grid * local_z

            @inbounds for i in eachindex(xs)
                x_grid = xs[i]
                p1 = p_row + x_grid * local_x

                acc = Complex{G}(0.0)
                for hit in hits
                    acc += beamlet_hit_field(hit, p1) * hit.sqrt_proj
                end
                field[i, j] = acc
            end
            _tick!(prog)
        end
    end
    return xs, zs, field
end

function electric_field(
        pd::Detector,
        hits::Vector{AstigmaticGaussianBeamletHit{G}};
        # kwargs
        n::Int = 100,
        crop_factor::Real = 3.0,
        num_spots::Int = 50,
        x_min = Inf,
        x_max = Inf,
        z_min = Inf,
        z_max = Inf,
        x0_shift::Real = 0,
        z0_shift::Real = 0,
        progress::Bool = true,
        kwargs...
) where {G}
    # Calculate autolims
    _x_min, _x_max, _z_min, _z_max = calc_local_lims(pd; crop_factor, num_spots)
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
    # Preallocate e-field matrix
    field = zeros(Complex{G}, n, n)
    # PD local coordinate axes (left-handed coords.)
    local_x = Point3(-orientation(pd)[:, 1])
    local_z = Point3(orientation(pd)[:, 3])
    origin = position(pd)
    # Calculate field superposition
    _with_progress(progress, length(zs), "Detector field: ") do prog
        Threads.@threads for j in eachindex(zs)
            z_grid = zs[j]
            # Hoist z-grid math out of inner loop
            p_row = origin + z_grid * local_z

            @inbounds for i in eachindex(xs)
                x_grid = xs[i]
                p1 = p_row + x_grid * local_x

                acc = Complex{G}(0.0)
                for hit in hits
                    acc += beamlet_hit_field(hit, p1) * hit.sqrt_proj
                end
                field[i, j] = acc
            end
            _tick!(prog)
        end
    end
    return xs, zs, field
end

function electric_field(
        pd::Detector,
        hits::Vector{RayHit{R}};
        # kwargs
        n::Int = 100,
        crop_factor::Real = 1,
        center::AbstractCenterAlgorithm = Centroid(),
        x_min = Inf,
        x_max = Inf,
        z_min = Inf,
        z_max = Inf,
        x0_shift::Real = 0,
        z0_shift::Real = 0,
        progress::Bool = true,
        kwargs...
) where {R}
    # automatically calculate limits
    _x_min, _x_max, _z_min, _z_max = calc_local_lims(pd; crop_factor, center)
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
    # Buffer field
    field = zeros(Complex{R}, n, n)

    # PD local coordinate axis
    orient = orientation(pd)
    @views e1, e2 = Point3(-orient[:, 1]), Point3(orient[:, 3])
    origin_pd = position(pd)

    hit_data = map(hits) do hit
        dir = direction(hit)
        p_hit = position(hit) + length(hit) * dir
        proj = projection_factor(hit)
        k = wavenumber(hit)
        opl = optical_path_length(hit)
        (p_hit, dir, proj, k, opl)
    end

    _with_progress(progress, length(zs), "Detector field: ") do prog
        Threads.@threads for j in eachindex(zs)
            z = zs[j]
            @inbounds for i in eachindex(xs)
                x = xs[i]
                # Global detector surface point coordinate
                p = origin_pd + x * e1 + z * e2
                # Add all field contributions
                acc = zero(complex(R))
                for (p_hit, dir, proj, k, opl) in hit_data
                    l = dot(p - p_hit, dir)
                    acc += proj * cis(k * (opl + l))
                end
                field[i, j] = acc
            end
            _tick!(prog)
        end
    end

    return xs, zs, field
end

function electric_field(
        pd::Detector,
        hits::Vector{PolarizedRayHit{R}};
        # kwargs
        n::Int = 100,
        crop_factor::Real = 1,
        center::AbstractCenterAlgorithm = Centroid(),
        x_min = Inf,
        x_max = Inf,
        z_min = Inf,
        z_max = Inf,
        x0_shift::Real = 0,
        z0_shift::Real = 0,
        progress::Bool = true,
        kwargs...
) where {R}
    # automatically calculate limits
    _x_min, _x_max, _z_min, _z_max = calc_local_lims(pd; crop_factor, center)
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
    # Buffer field (global E-field vector per grid point)
    field = fill(zero(Point3{Complex{R}}), n, n)

    # PD local coordinate axis
    orient = orientation(pd)
    @views e1, e2 = Point3(-orient[:, 1]), Point3(orient[:, 3])
    origin_pd = position(pd)

    hit_data = map(hits) do hit
        dir = direction(hit)
        p_hit = position(hit) + length(hit) * dir
        E0 = polarization(hit)
        k = wavenumber(hit)
        opl = optical_path_length(hit)
        (p_hit, dir, E0, k, opl)
    end

    _with_progress(progress, length(zs), "Detector field: ") do prog
        Threads.@threads for j in eachindex(zs)
            z = zs[j]
            @inbounds for i in eachindex(xs)
                x = xs[i]
                # Global detector surface point coordinate
                p = origin_pd + x * e1 + z * e2
                # Coherently add all vector field contributions
                acc = zero(Point3{Complex{R}})
                for (p_hit, dir, E0, k, opl) in hit_data
                    l = dot(p - p_hit, dir)
                    acc += E0 * cis(k * (opl + l))
                end
                field[i, j] = acc
            end
            _tick!(prog)
        end
    end

    return xs, zs, field
end

"""
    intensity(pd::Detector, Z::Number = Z_vacuum; kwargs...)

Calculates the intensity distribution on the [`Detector`](@ref) via

```math
I = \\frac{|E|^2}{2 \\cdot Z}
```

where E is the electric field value and Z is the wave impedance. In general, vacuum wave impedance is assumed.
This function returns a tuple `(x, y, I)`. For more information on the available
keyword arguments, refer to the [`electric_field`](@ref) documentation.
"""
function intensity(pd::Detector, Z::Number = Z_vacuum; kwargs...)
    x, y, E = electric_field(pd; kwargs...)
    return x, y, intensity.(E, Z)
end

"""
    optical_power(pd::Detector; kwargs...)

Calculates the total optical power on `pd` in [W] by integration over the local intensity.
For more information on keyword argument options, refer to the [`electric_field`](@ref) docs.
"""
function optical_power(pd::Detector; kwargs...)
    x, y, I = intensity(pd; kwargs...)
    return trapz((x, y), I)
end
