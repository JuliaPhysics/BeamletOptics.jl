struct PSFData{T}
    hit::Point3{T}
    dir::Point3{T}
    opl::T
    proj::T
    λ::T
end

position(data::PSFData) = data.hit
direction(data::PSFData) = data.dir
optical_path_length(data::PSFData) = data.opl
wavelength(data::PSFData) = data.λ
projection_factor(data::PSFData) = data.proj

function PSFData(hit::AbstractArray{H}, dir::AbstractArray{D}, opl::O, proj::P, lambda::L) where {H, D, O, P, L}
    T = promote_type(H, D, O, P, L)
    return PSFData{T}(T.(hit), T.(dir), T(opl), T(proj), T(lambda))
end

mutable struct PSFDetector{T} <: AbstractDetector{T, Mesh{T}}
    const shape::Mesh{T}
    data::Vector{PSFData{T}}
    solved::Bool
end

function reset_detector!(psf::PSFDetector{T}) where T
    psf.data = Vector{PSFData{T}}()
    psf.solved = false
    return nothing
end

Base.push!(psf::PSFDetector{T}, new::PSFData{T}) where T = push!(psf.data, new)

function PSFDetector(width::W) where W<:AbstractFloat
    # Spawn mesh, align with neg. y-axis, empty data field
    shape = QuadraticFlatMesh(width)
    zrotate3d!(shape, π)
    data = Vector{PSFData{W}}()
    return PSFDetector(shape, data, false)
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
    push!(psf, PSFData(hit_3D, dir_3D, l, proj, λ))
    return nothing
end

function intensity(psf::PSFDetector, sz::Real=1e-5, n::Int=100)
    xs = LinRange(-sz/2, sz/2, n)
    ys = LinRange(-sz/2, sz/2, n)
    field = zeros(Complex, n, n)
    
    orient = orientation(psf)
    e1 = @view orient[:, 1]   # local x-axis
    e2 = @view orient[:, 3]   # local y-axis
    origin_pd = position(psf)
    
    for j in eachindex(ys)
        y = ys[j]
        for i in eachindex(xs)
            x = xs[i]
            # Global detector surface point coordinate
            p1 = origin_pd + x * e1 + y * e2
            # Add all field contributions
            for hit in psf.data
                # Calculate local offset and phase
                l1 = dot(p1 - position(hit), direction(hit))
                k = 2π / wavelength(hit)
                amplitude = projection_factor(hit) * exp(im * k * (optical_path_length(hit) + l1))
                field[i, j] += amplitude
            end
        end
    end
    return xs, ys, abs2.(field)
end