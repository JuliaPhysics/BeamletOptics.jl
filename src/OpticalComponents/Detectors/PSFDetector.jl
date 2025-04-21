struct PSFData{T}
    hit::Point3{T}
    dir::Point3{T}
    opl::T
    proj::T
    k::T
end

position(data::PSFData) = data.hit
direction(data::PSFData) = data.dir
optical_path_length(data::PSFData) = data.opl
wavenumber(data::PSFData) = data.k
projection_factor(data::PSFData) = data.proj

function PSFData(hit::AbstractArray{H}, dir::AbstractArray{D}, opl::O, proj::P, wavenumber::K) where {H, D, O, P, K}
    T = promote_type(H, D, O, P, K)
    return PSFData{T}(T.(hit), T.(dir), T(opl), T(proj), T(wavenumber))
end

struct PSFDetector{T} <: AbstractDetector{T, Mesh{T}}
    shape::Mesh{T}
    data::Vector{PSFData{T}}
    solved::Ref{Bool}
end

empty!(psf::PSFDetector) = (empty!(psf.data); psf.solved[] = false)

Base.push!(psf::PSFDetector, new::PSFData) = push!(psf.data, new)

function PSFDetector(width::W) where W <: AbstractFloat
    # Spawn mesh, align with neg. y-axis, empty data field
    shape = QuadraticFlatMesh(width)
    zrotate3d!(shape, π)
    data = Vector{PSFData{W}}()
    return PSFDetector(shape, data, Ref(false))
end

function interact3d(::AbstractSystem, psf::PSFDetector, beam::Beam{T, Ray{T}}, ray::Ray{T}) where T
    # Global hit pos
    hit_3D = position(ray) + length(ray) * direction(ray)
    # Global hit dir
    dir_3D = direction(ray)
    # Optical path length and phase
    l = optical_path_length(beam)
    λ = wavelength(ray)
    # Projection factor
    proj = abs(dot(dir_3D, normal3d(intersection(ray))))
    push!(psf, PSFData(hit_3D, dir_3D, l, proj, 2π/λ))
    return nothing
end

function calc_local_pos(psf::PSFDetector{T}) where T
    hits_2D = Vector{Point2{T}}(undef, length(psf.data))
    for (i, h) in enumerate(psf.data)
        hit = position(h)
        loc_pos = hit - position(psf)
        @views x = dot(loc_pos, Point3(orientation(psf)[:,1]))
        @views z = dot(loc_pos, Point3(orientation(psf)[:,3]))
        hits_2D[i] = Point2(T(x), T(z))
    end
    return hits_2D
end

"""
    calc_local_lims(psf; crop_factor=1, center=:centroid)

Compute a symmetric [x_min,x_max]×[z_min,z_max] box around the PSF’s weighted centroid.

• If `center==:centroid` (the default), uses
    x0 = ∑ wᵢ·xᵢ / ∑ wᵢ,  z0 = ∑ wᵢ·yᵢ / ∑ wᵢ
  with wᵢ = projection_factor.
• If `center==:bbox`, falls back to the midpoint of [min,max].

Returns `(x_min, x_max, z_min, z_max)`.
"""
function calc_local_lims(psf::PSFDetector{T};
                         crop_factor::Real=one(T),
                         center::Symbol = :centroid) where T

    # get all hit‐points in local (x,y)
    hits = calc_local_pos(psf)
    xs = getindex.(hits,1)
    zs = getindex.(hits,2)

    # choose center
    if center == :centroid
        w_sum = sum(projection_factor, psf.data)
        x0 = sum(x -> (projection_factor(x[1]) * x[2]), zip(psf.data, xs)) / w_sum
        z0 = sum(x -> (projection_factor(x[1]) * x[2]), zip(psf.data, zs)) / w_sum
    else
        x0 = (minimum(xs) + maximum(xs)) / 2
        z0 = (minimum(zs) + maximum(zs)) / 2
    end

    # half‐widths
    dx = maximum(x->abs(x - x0), xs)
    dy = maximum(y->abs(y - z0), zs)
    hwx, hwy = dx*crop_factor, dy*crop_factor

    return x0 - hwx, x0 + hwx, z0 - hwy, z0 + hwy
end

function intensity(psf::PSFDetector{T};
        n::Int=100,
        crop_factor::Real=1,
        normalization=:strehl,
        x_min = Inf,
        x_max = Inf,
        z_min = Inf,
        z_max = Inf,
        x0_shift::Real=0,
        z0_shift::Real=0) where T
    # automatically calculate limits
    _x_min, _x_max, _z_min, _z_max = calc_local_lims(psf; crop_factor)
    if x_min != Inf && x_max != Inf
        _x_min = x_min
        _x_max = x_max
    end
    if z_min != Inf && z_max != Inf
        _z_min = x_min
        _z_max = x_max
    end
    xs = LinRange(_x_min, _x_max, n) .+ x0_shift
    zs = LinRange(_z_min, _z_max, n) .+ z0_shift
    # Buffer field
    field = zeros(Complex{T}, n, n)

    # PD local coordinate axis
    orient = orientation(psf)
    @views e1, e2 = Point3(orient[:, 1]), Point3(orient[:, 3])
    origin_pd = position(psf)

    Threads.@threads for j in eachindex(zs)
        z = zs[j]
        @inbounds for i in eachindex(xs)
            x = xs[i]
            # Global detector surface point coordinate
            p = origin_pd + x * e1 + z * e2
            # Add all field contributions
            acc = zero(complex(T))
            @inbounds @simd for h in psf.data
                l = dot(p - position(h), direction(h))
                acc += projection_factor(h) * cis(wavenumber(h) * (optical_path_length(h) + l))
            end
            field[i,j] = acc
        end
    end

    _intensity = abs2.(field)

    if normalization === :strehl
        S = strehl(psf)
        _intensity .*= S / maximum(_intensity)
    end

    return xs, zs, _intensity
end

function strehl(psf::PSFDetector{T}) where {T}
    k = psf.data[1].k
    ops = optical_path_length.(psf.data)

    μ = sum(ops) / length(ops)
    num = zero(complex(T))
    den = zero(T)
    for h in psf.data
        w = projection_factor(h)
        num += w * cis(k * (optical_path_length(h) - μ))
        den += w
    end
    return abs2(num / den)
end
