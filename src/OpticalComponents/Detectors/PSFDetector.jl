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

reset!(psf::PSFDetector) = (empty!(psf.data); psf.solved[] = false)    

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

function intensity(psf::PSFDetector{T}, sz::Real=T(1e-5), n::Int=100) where T
    xs = LinRange(-sz/2, sz/2, n)
    ys = LinRange(-sz/2, sz/2, n)
    field = zeros(Complex{T}, n, n)
    
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
            @inbounds @simd for h in det.data
                l = dot(p - position(h), direction(h))
                acc += projection_factor(h) * cis(wavenumber(h) * (optical_path_length(h) + l))
            end
            field[i,j] = acc
        end
    end
    return xs, ys, abs2.(field)
end