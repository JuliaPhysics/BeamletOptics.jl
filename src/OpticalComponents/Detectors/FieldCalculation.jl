"""
    electric_field(pd::Detector; kwargs...)

Compute a two‐dimensional electric field based on incoming rays or beams as captured by a [`Detector`](@ref).
The returned E-field map is sampled on a regular `n×n` grid in the detector’s local (x,z)-plane.
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

- `center::Symbol=:centroid`
  How the sampling window is centred.
  `:centroid` uses the projection‑weighted centroid,
  `:bbox` uses the geometric mid‑point of the bounding box.

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
  in the detector’s local x and z axes,
- `E::Matrix{Complex{T}}` is the corresponding raw/unscaled intensity map
"""
electric_field(d::Detector; kwargs...) = electric_field(d, hits(d); kwargs...)

electric_field(d::Detector, ::Nothing; kwargs...) = throw(ErrorException("No hits available on detector."))

function electric_field(
        pd::Detector,
        hits::Vector{GaussianBeamletHit{G}};
        # kwargs
        n::Int=100,
        crop_factor::Real=1.5,
        num_spots::Int=50,
        x_min = Inf,
        x_max = Inf,
        z_min = Inf,
        z_max = Inf,
        x0_shift::Real=0,
        z0_shift::Real=0
    ) where G
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
    local_x = Point3(-orientation(pd)[:,1])  # flipped sign due to rotated pd mesh
    local_z = Point3(orientation(pd)[:,3])
    origin = position(pd)
    # Calculate field superposition
    Threads.@threads for j in eachindex(zs) # FIXME row column major order?
        z = zs[j]
        @inbounds for i in eachindex(xs)
            x = xs[i] 
            # Transform point p on PD into world coordinates
            p1 = origin + x * local_x + z * local_z
            # Add E-field contribution for each hit
            for hit in hits
                p0 = position(hit)
                d0 = direction(hit)
                # important: calculate l0 with geometric length for local r, z
                l0 = hit.l0
                proj = projection_factor(hit)
                # Find projection of p1 onto Gaussian optical axis, i.e. local r and z
                l1 = dot(p1 - p0, d0)
                p2 = p0 + l1 * d0
                r_gb = norm(p1 - p2)
                z_gb = l0 + l1
                # Add field contribution, projection factor accounts for beam spot stretching
                field[i, j] += electric_field(hit.gauss, r_gb, z_gb) * sqrt(proj)
            end
        end
    end
    return xs, zs, field
end

function electric_field(
        pd::Detector,
        hits::Vector{RayHit{R}};
        # kwargs
        n::Int=100,
        crop_factor::Real=1,
        center::Symbol=:centroid,
        x_min = Inf,
        x_max = Inf,
        z_min = Inf,
        z_max = Inf,
        x0_shift::Real=0,
        z0_shift::Real=0
    ) where R
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

    Threads.@threads for j in eachindex(zs)
        z = zs[j]
        @inbounds for i in eachindex(xs)
            x = xs[i]
            # Global detector surface point coordinate
            p = origin_pd + x * e1 + z * e2
            # Add all field contributions
            acc = zero(complex(R))
            @inbounds @simd for hit in hits
                l = dot(p - hit_point(hit), direction(hit))
                acc += projection_factor(hit) * cis(wavenumber(hit) * (optical_path_length(hit) + l))
            end
            field[i,j] = acc
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