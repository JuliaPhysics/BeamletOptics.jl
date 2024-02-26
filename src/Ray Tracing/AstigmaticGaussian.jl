mutable struct AstigmaticGaussianBeamlet{T} <: AbstractBeam{T, PolarizedRay{T}}
    id::UUID
    c::Beam{T, PolarizedRay{T}}       # chief
    wx::Beam{T, Ray{T}}               # waist x
    wy::Beam{T, Ray{T}}               # waist y
    dx::Beam{T, Ray{T}}               # div. x
    dy::Beam{T, Ray{T}}               # div. y
    parent::Nullable{AstigmaticGaussianBeamlet{T}}
    children::Vector{AstigmaticGaussianBeamlet{T}}
end

function AstigmaticGaussianBeamlet(
    c::Beam{T, PolarizedRay{T}},
    wx::Beam{T, Ray{T}},
    wy::Beam{T, Ray{T}},
    dx::Beam{T, Ray{T}},
    dy::Beam{T, Ray{T}}
    ) where {T <: Real}
return AstigmaticGaussianBeamlet{T}(uuid4(), c, wx, wy, dx, dy,
    nothing,
    Vector{AstigmaticGaussianBeamlet{T}}())
end

function AstigmaticGaussianBeamlet(
    position,
    direction,
    λ = 1000e-9,
    w0 = 1e-3;
    M2 = 1,
    E0 = [0, 0, 1])
    # Create orthogonal vectors for construction purposes (right-handed)
    s1 = [0,1,0] # normal3d(direction)
    s2 = cross(direction, s1)
    # Divergence angle in rad
    θ = divergence_angle(λ, w0, M2)
    # Waist rays
    wx = Ray(position + s1 * w0, direction, λ)
    wy = Ray(position + s2 * w0, direction, λ)
    # Divergence ray
    div_dir_x = normalize(direction + s1 * tan(θ))
    div_dir_y = normalize(direction + s2 * tan(θ))
    dx = Ray(position, div_dir_x, λ)
    dy = Ray(position, div_dir_y, λ)
    # Chief ray
    c = PolarizedRay(position, direction, λ, E0)
    return AstigmaticGaussianBeamlet(Beam(c), Beam(wx), Beam(wy), Beam(dx), Beam(dy))
end

_last_beam_intersection(agb::AstigmaticGaussianBeamlet) = intersection(last(rays(agb.c)))

point_on_beam(agb::AstigmaticGaussianBeamlet, t::Real) = point_on_beam(agb.c, t)

function gauss_parameters(agb::AstigmaticGaussianBeamlet, z::Real)
    p0, i = point_on_beam(agb, z)

    n = refractive_index(agb.c.rays[i])
    λ = wavelength(agb.c.rays[i])

    wx, ~, ~, ~, Hnx = gauss_parameters(agb.c.rays[i], agb.wx.rays[i], agb.dx.rays[i], p0)
    wy, ~, ~, ~, Hny = gauss_parameters(agb.c.rays[i], agb.wy.rays[i], agb.dy.rays[i], p0)

    return wx, wy, Hnx*n, Hny*n
end

function gauss_parameters(agb::AstigmaticGaussianBeamlet{G}, zs::AbstractArray) where {G}
    wx = Vector{G}(undef, length(zs))
    wy = Vector{G}(undef, length(zs))
    Hnx = Vector{G}(undef, length(zs))
    Hny = Vector{G}(undef, length(zs))
    for (i, z) in enumerate(zs)
        wx[i], wy[i], Hnx[i], Hny[i] = gauss_parameters(agb, z)
    end
    return wx, wy, Hnx, Hny
end
