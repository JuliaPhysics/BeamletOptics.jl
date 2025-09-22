electric_field(d::Detector; kwargs...) = electric_field(d, hits(d); kwargs...)

function electric_field(
        pd::Detector,
        hits::Vector{GaussianBeamletHit{G}};
        n::Int=100,
        kwargs...
    ) where G
    # Preallocate transforms
    T = transpose(orientation(shape(pd)))
    p = position(shape(pd))
    # Preallocate x-/y-buffers and field
    xs = ys = LinRange(-pd.edgln/2, pd.edgln/2, n)
    field = zeros(Complex{G}, n, n)
    # Calculate field superposition
    Threads.@threads for j in eachindex(ys) # FIXME row column major order?
        y = ys[j]
        @inbounds for i in eachindex(xs)
            x = xs[i] 
            # Transform point p on PD into world coordinates
            p1 = Point3(
                T[1, 1] * x + T[1, 3] * y + p[1],
                T[2, 1] * x + T[2, 3] * y + p[2],
                T[3, 1] * x + T[3, 3] * y + p[3]
            )
            # Add E-field contribution for each hit
            for hit in hits
                p0 = position(hit)
                d0 = direction(hit)
                l0 = length(hit)
                proj = projection_factor(hit)
                # Find projection of p1 onto Gaussian optical axis, i.e. local r and z
                l1 = dot(p1 - p0, d0)
                p2 = p0 + l1 * d0
                r = norm(p1 - p2)
                z = l0 + l1
                # Add field contribution, projection factor accounts for beam spot stretching
                field[i, j] += electric_field(hit.gauss, r, z) * sqrt(proj)
            end
        end
    end
    return xs, ys, field
end

function electric_field(
        pd::Detector,
        hits::Vector{RayHit{R}};
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
    @views e1, e2 = Point3(orient[:, 1]), Point3(orient[:, 3])
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
                l = dot(p - position(hit), direction(hit))
                acc += projection_factor(hit) * cis(wavenumber(hit) * (optical_path_length(hit) + l))
            end
            field[i,j] = acc
        end
    end

    return xs, zs, field
end

function intensity(pd::Detector, Z::Number = Z_vacuum; kwargs...)
    x, y, E = electric_field(pd; kwargs...)
    return x, y, intensity.(E, Z)
end

"""
    optical_power(pd::Detector)

Calculates the total optical power on `pd` in [W] by integration over the local intensity.
"""
function optical_power(pd::Detector)
    x, y, I = intensity(pd)
    return trapz((x, y), I)
end