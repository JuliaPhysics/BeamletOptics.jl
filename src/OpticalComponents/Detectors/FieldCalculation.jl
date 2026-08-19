"""
    electric_field(pd::Detector; kwargs...)

Compute a two‐dimensional electric field based on incoming rays or beams as captured by a [`Detector`](@ref).
The returned E-field map is sampled on a regular `n×n` grid in the detector's local (x,z)-plane.
Note that the `pd` local coordinates are given in a (x, z) basis where the normal vector forms a left-handed system.

!!! note "Resetting detectors"
    Be sure to call `empty!(pd)` before each new measurement if reusing the same detector.

# Keyword Arguments

The following generic kwargs can be used for all hit types:

- `n::Int=100`
  Number of sample points per axis.
- `crop_factor::Real=1`
  Scales the width of the sampling window returned by
  [`calc_local_lims`](@ref); values >1 expand, <1 shrink.
- `x_min, x_max, z_min, z_max`
  Manually override the sampling bounds in the local x or z directions.
  If left as `Inf`, the bounds from `calc_local_lims` are used.
- `x0_shift::Real=0, z0_shift::Real=0`
  Applies a constant offset to the entire x or z coordinate arrays,
  useful for recentring or testing alignment.

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
- `E::Matrix{Complex{T}}` is the corresponding raw/unscaled intensity map
"""
electric_field(d::Detector; kwargs...) = electric_field(d, hits(d); kwargs...)

function electric_field(d::Detector, ::Nothing; kwargs...)
    throw(ErrorException("No hits available on detector."))
end

function electric_field(
        pd::Detector,
        hits::Vector{GaussianBeamletHit{G}};
        # kwargs
        n::Int = 100,
        crop_factor::Real = 1.5,
        num_spots::Int = 50,
        x_min = Inf,
        x_max = Inf,
        z_min = Inf,
        z_max = Inf,
        x0_shift::Real = 0,
        z0_shift::Real = 0,
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
    for hit in hits
        # Pre-compute hit constants
        p0 = hit.p0
        d0 = hit.d0
        gauss = hit.gauss
        l0 = hit.l0
        id = hit.id
        sqrt_proj = hit.sqrt_proj

        Threads.@threads for j in eachindex(zs)
            z_grid = zs[j]
            p_row = origin + z_grid * local_z

            @inbounds for i in eachindex(xs)
                x_grid = xs[i]
                p1 = p_row + x_grid * local_x

                v = p1 - p0
                l1 = _pseudo_dot(v, d0)
                r_vec = v - l1 * d0
                r = norm(r_vec)
                z = l0 + l1

                field[i, j] += electric_field(gauss, r, z; hint = (p0 + l1 * d0, id)) * sqrt_proj
            end
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

    # Calculate field superposition with hoisted hit invariants
    for hit in hits
        d0 = hit.d0
        h1 = hit.h1
        u1 = hit.u1
        h2 = hit.h2
        u2 = hit.u2
        p0 = hit.p0
        l0 = hit.l0
        k0 = hit.k0
        Δl = hit.Δl
        n_eff = hit.n_eff
        area_ref = hit.area_ref
        E_ref_amp = hit.E_ref_amp
        sqrt_proj = hit.sqrt_proj

        Threads.@threads for j in eachindex(zs)
            z_grid = zs[j]
            p_row = origin + z_grid * local_z

            @inbounds for i in eachindex(xs)
                x_grid = xs[i]
                p1 = p_row + x_grid * local_x

                v = p1 - p0
                l1 = _pseudo_dot(v, d0)
                r_vec = v - l1 * d0

                h1_z = h1 + l1 * u1
                h2_z = h2 + l1 * u2
                area_z = _pseudo_cross2d(h1_z, h2_z, d0)
                if abs(area_z) < 1e-25
                    area_z = Complex{G}(1e-25, 1e-25)
                end

                ξ1 = _pseudo_cross2d(h1_z, r_vec, d0)
                ξ2 = _pseudo_cross2d(h2_z, r_vec, d0)
                w = (ξ1 * _pseudo_dot(u2, r_vec) - ξ2 * _pseudo_dot(u1, r_vec)) / (2 * area_z)

                phase_corr = (n_eff - 1) * l1
                z_total = l0 + l1
                ψ = sqrt(area_ref / area_z) * cis(k0 * (z_total + w + Δl + phase_corr))

                field[i, j] += (E_ref_amp * ψ) * sqrt_proj
            end
        end
    end
    return xs, zs, field
end

function electric_field(
        pd::Detector,
        hits::Vector{<:AbstractRayHit{R}};

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
