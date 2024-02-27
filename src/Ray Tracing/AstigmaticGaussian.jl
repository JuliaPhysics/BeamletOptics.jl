mutable struct AstigmaticGaussianBeamlet{T} <: AbstractBeam{T, PolarizedRay{T}}
    id::UUID
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
return AstigmaticGaussianBeamlet{T}(uuid4(), c, wxp, wxm, wyp, wym, dxp, dxm, dyp, dym,
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
    support = nothing)
    # Create orthogonal vectors for construction purposes (right-handed)
    if isnothing(support)
        s1 = normal3d(direction)
    else
        s1 = normalize(support)
    end
    s2 = cross(direction, s1)
    # Divergence angle in rad
    θ = divergence_angle(λ, w0, M2)
    # Waist rays
    wxp = Ray(position + s1 * w0, direction, λ)
    wxm = Ray(position - s1 * w0, direction, λ)
    wyp = Ray(position + s2 * w0, direction, λ)
    wym = Ray(position - s2 * w0, direction, λ)
    # Divergence ray
    div_dir_xp = normalize(direction + s1 * tan(θ))
    div_dir_xm = normalize(direction - s1 * tan(θ))
    div_dir_yp = normalize(direction + s2 * tan(θ))
    div_dir_ym = normalize(direction - s2 * tan(θ))
    dxp = Ray(position, div_dir_xp, λ)
    dxm = Ray(position, div_dir_xm, λ)
    dyp = Ray(position, div_dir_yp, λ)
    dym = Ray(position, div_dir_ym, λ)
    # Chief ray
    c = PolarizedRay(position, direction, λ, E0)
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

_last_beam_intersection(agb::AstigmaticGaussianBeamlet) = intersection(last(rays(agb.c)))

point_on_beam(agb::AstigmaticGaussianBeamlet, t::Real) = point_on_beam(agb.c, t)

function gauss_parameters(agb::AstigmaticGaussianBeamlet, z::Real)
    p0, i = point_on_beam(agb, z)

    n = refractive_index(agb.c.rays[i])
    λ = wavelength(agb.c.rays[i])

    wxp, ~, ~, ~, Hnxp = gauss_parameters(agb.c.rays[i], agb.wxp.rays[i], agb.dxp.rays[i], p0)
    wxm, ~, ~, ~, Hnxm = gauss_parameters(agb.c.rays[i], agb.wxm.rays[i], agb.dxm.rays[i], p0)
    wyp, ~, ~, ~, Hnyp = gauss_parameters(agb.c.rays[i], agb.wyp.rays[i], agb.dyp.rays[i], p0)
    wym, ~, ~, ~, Hnym = gauss_parameters(agb.c.rays[i], agb.wym.rays[i], agb.dym.rays[i], p0)

    # generate basis
    r = agb.dxp.rays[i]
    l = line_plane_distance3d(p0, direction(agb.c.rays[i]), position(r), direction(r))
    x = position(r) + l * direction(r) - p0
    x = normalize(x)

    r = agb.dyp.rays[i]
    l = line_plane_distance3d(p0, direction(agb.c.rays[i]), position(r), direction(r))
    y = position(r) + l * direction(r) - p0
    y = normalize(y)
    
    z = direction(agb.c.rays[i])

    return wxp, wxm, wyp, wym, [p0 x*wxp y*wyp]
end

function gauss_parameters(agb::AstigmaticGaussianBeamlet{G}, zs::AbstractArray) where {G}
    wxp = Vector{G}(undef, length(zs))
    wxm = Vector{G}(undef, length(zs))
    wyp = Vector{G}(undef, length(zs))
    wym = Vector{G}(undef, length(zs))
    for (i, z) in enumerate(zs)
        wxp[i], wxm[i], wyp[i], wym[i], ~ = gauss_parameters(agb, z)
    end
    return wxp, wxm, wyp, wym
end
