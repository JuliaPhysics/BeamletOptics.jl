mutable struct AstigmaticGaussianBeamlet{T} <: AbstractBeam{T, PolarizedRay{T}}
    c::Beam{T, PolarizedRay{T}}       # chief
    wxp::Beam{T, Ray{T}}              # waist x +
    wxm::Beam{T, Ray{T}}              # waist x -
    wyp::Beam{T, Ray{T}}              # waist y +
    wym::Beam{T, Ray{T}}              # waist y -
    dxp::Beam{T, Ray{T}}              # div. x +
    dxm::Beam{T, Ray{T}}              # div. x -
    dyp::Beam{T, Ray{T}}              # div. y +
    dym::Beam{T, Ray{T}}              # div. y -
    parent::Nullable{AstigmaticGaussianBeamlet{T}}
    children::Vector{AstigmaticGaussianBeamlet{T}}
end

function AstigmaticGaussianBeamlet(
    c::Beam{T, PolarizedRay{T}},
    wxp::Beam{T, Ray{T}},
    wxm::Beam{T, Ray{T}},
    wyp::Beam{T, Ray{T}},
    wym::Beam{T, Ray{T}},
    dxp::Beam{T, Ray{T}},
    dxm::Beam{T, Ray{T}},
    dyp::Beam{T, Ray{T}},
    dym::Beam{T, Ray{T}}
    ) where {T <: Real}
return AstigmaticGaussianBeamlet{T}(c, wxp, wxm, wyp, wym, dxp, dxm, dyp, dym,
    nothing,
    Vector{AstigmaticGaussianBeamlet{T}}())
end

function AstigmaticGaussianBeamlet(
    position,
    direction,
    λ = 1000e-9,
    w0 = 1e-3;
    M2 = 1,
    E0 = [0, 0, 1],
    support = nothing,
    z0 = 0)
    # Create orthogonal vectors for construction purposes (right-handed)
    direction = normalize(direction)
    if isnothing(support)
        s1 = normal3d(direction)
    else
        if !isorthogonal(direction, support)
            error("Ray direction and support vector must be orthogonal!")
        end
        s1 = normalize(support)
    end
    s2 = cross(direction, s1)
    # Divergence angle in rad
    θ = divergence_angle(λ, w0, M2)
    # Waist rays
    wxp = Ray(position + s1 * w0 + z0 * direction, direction, λ)
    wxm = Ray(position - s1 * w0 + z0 * direction, direction, λ)
    wyp = Ray(position + s2 * w0 + z0 * direction, direction, λ)
    wym = Ray(position - s2 * w0 + z0 * direction, direction, λ)
    # Divergence ray
    div_dir_xp = normalize(direction + s1 * tan(θ))
    div_dir_xm = normalize(direction - s1 * tan(θ))
    div_dir_yp = normalize(direction + s2 * tan(θ))
    div_dir_ym = normalize(direction - s2 * tan(θ))
    dxp = Ray(position + div_dir_xp * z0, div_dir_xp, λ)
    dxm = Ray(position + div_dir_xm * z0, div_dir_xm, λ)
    dyp = Ray(position + div_dir_yp * z0, div_dir_yp, λ)
    dym = Ray(position + div_dir_ym * z0, div_dir_ym, λ)
    # Chief ray
    c = PolarizedRay(position + z0 * direction, direction, λ, E0)
    return AstigmaticGaussianBeamlet(
        Beam(c),
        Beam(wxp),
        Beam(wxm),
        Beam(wyp),
        Beam(wym),
        Beam(dxp),
        Beam(dxm),
        Beam(dyp),
        Beam(dym)
    )
end

Base.length(agb::AstigmaticGaussianBeamlet) = length(agb.c)
optical_path_length(agb::AstigmaticGaussianBeamlet) = optical_path_length(agb.c)

isentering(agb::AstigmaticGaussianBeamlet, id::Int) = isentering(rays(agb.c)[id])

_last_beam_intersection(agb::AstigmaticGaussianBeamlet) = intersection(last(rays(agb.c)))

point_on_beam(agb::AstigmaticGaussianBeamlet, t::Real) = point_on_beam(agb.c, t)