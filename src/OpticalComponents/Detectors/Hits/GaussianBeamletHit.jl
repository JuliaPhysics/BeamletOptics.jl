"""
    GaussianBeamletHit{T} <: AbstractBeamletHit{T}

Stores [`GaussianBeamlet`] hit data.
"""
struct GaussianBeamletHit{T} <: AbstractBeamletHit{T}
    # ray pos and dir
    pos_c::Point3{T}
    pos_w::Point3{T}
    pos_d::Point3{T}
    dir_c::Point3{T}
    dir_w::Point3{T}
    dir_d::Point3{T}
    # surface parameters
    nml::Point3{T}
    hit::Point3{T}
    projection::T
    # beam parameters
    λ::T
    n::T
    l0::T
    E0::Complex{T}
    Δϕ::T
    w0::T
    z0::T
end

function GaussianBeamletHit(gauss::GaussianBeamlet, index::Int)
    # select rays at detector
    chief = gauss.chief.rays[index]
    div = gauss.divergence.rays[index]
    waist = gauss.waist.rays[index]
    # fetch beam ray parameters
    pos_c = position(chief)
    pos_w = position(waist)
    pos_d = position(div)
    dir_c = direction(chief)
    dir_w = direction(div)
    dir_d = direction(waist)
    normal = normal3d(intersection(chief))
    # hit point and projection factor
    hit_pos = position(chief) + length(chief) * direction(chief)
    proj = abs(dot(direction(chief), normal3d(intersection(chief))))
    # E-field parameters (base length and optical index phase)
    l0 = length(gauss) - length(chief)
    Δl = optical_path_length(gauss) - length(gauss)
    z0 = length(gauss)
    Δϕ = Δl / wavelength(gauss) * 2π
    return GaussianBeamletHit(
        pos_c, pos_w, pos_d,
        dir_c, dir_w, dir_d,
        normal, hit_pos, proj, 
        wavelength(gauss), refractive_index(chief), l0, 
        electric_field(gauss), Δϕ, beam_waist(gauss), z0
    )
end

position(hit::GaussianBeamletHit) = hit.pos_c
direction(hit::GaussianBeamletHit) = hit.dir_c
normal3d(hit::GaussianBeamletHit) = hit.nml

hit_point(hit::GaussianBeamletHit) = hit.hit

projection_factor(hit::GaussianBeamletHit) = hit.projection

gauss_parameters(hit::GaussianBeamletHit, p0::AbstractArray) = gauss_parameters(
    p0, hit.pos_w, hit.pos_d, hit.dir_c, hit.dir_w, hit.dir_d, hit.λ, hit.n
)

function electric_field(hit::GaussianBeamletHit, r::Real, z::Real)
    point = hit_point(hit) + direction(hit) * (hit.z0 - z)


    
    w, R, ψ, w0 = gauss_parameters(hit, point)
    k = wavenumber(hit.λ)
    # Calculate new local field strength based on E0*w0 = const.
    E0 = hit.E0 * (hit.w0 / w0)
    return electric_field(r, z, E0, w0, w, k, ψ, R) * exp(im*hit.Δϕ)
end

