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
        x = dot(loc_pos, orientation(psf)[:,1])
        z = dot(loc_pos, orientation(psf)[:,3])
        hits_2D[i] = Point2(T(x), T(z))
    end
    return hits_2D
end

function calc_local_lims(psf::PSFDetector{T}; crop_factor::Real=1) where T
    hits_2D = calc_local_pos(psf)
    x_min = 0
    x_max = 0
    y_min = 0
    y_max = 0    
    for (i, p) in enumerate(hits_2D)
        # ignore the center beam due to small angle spot diagram "bug"
        # refer to https://github.com/StackEnjoyer/BeamletOptics.jl/issues/11
        if i == 1
            continue
        end
        # replace 0 vals
        if i == 2
            x_min = p[1]
            x_max = p[1]
            y_min = p[2]
            y_max = p[2]
            continue
        end
        p[1] < x_min ? x_min = p[1] : nothing
        p[1] > x_max ? x_max = p[1] : nothing
        p[2] < y_min ? y_min = p[2] : nothing
        p[2] > y_max ? y_max = p[2] : nothing
    end
    # Apply crop factor around box center
    x_halfwidth = 0.5(x_max - x_min)
    y_halfwidth = 0.5(y_max - y_min)
    x_center = x_min + x_halfwidth
    y_center = y_min + y_halfwidth
    x_min = x_center - x_halfwidth * crop_factor
    x_max = x_center + x_halfwidth * crop_factor
    y_min = y_center - y_halfwidth * crop_factor
    y_max = y_center + y_halfwidth * crop_factor
    return x_min, x_max, y_min, y_max
end


function intensity(psf::PSFDetector{T}; n::Int=100, crop_factor::Real=1) where T
    # automatically calculate limits
    x_min, x_max, y_min, y_max = calc_local_lims(psf; crop_factor)
    xs = LinRange(x_min, x_max, n)
    ys = LinRange(y_min, y_max, n)
    # Buffer field
    field = zeros(Complex{T}, n, n)
    
    # PD local coordinate axis
    orient = orientation(psf)
    @views e1, e2 = orient[:, 1], orient[:, 3]
    origin_pd = position(psf)
    
    Threads.@threads for j in eachindex(ys)
        y = ys[j]
        @inbounds for i in eachindex(xs)
            x = xs[i]
            # Global detector surface point coordinate
            p = origin_pd + x * e1 + y * e2
            # Add all field contributions
            acc = zero(Complex{T})
            @inbounds @simd for h in psf.data
                l = dot(p - position(h), direction(h))
                acc += projection_factor(h) * cis(wavenumber(h) * (optical_path_length(h) + l))
            end
            field[i,j] = acc
        end
    end
    return xs, ys, abs2.(field)
end