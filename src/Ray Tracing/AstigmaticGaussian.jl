mutable struct AstigmaticGaussianBeamlet{T} <: SCDI.AbstractBeam{T, SCDI.PolarizedRay{T}}
    id::UUID
    c::SCDI.Beam{T, SCDI.PolarizedRay{T}}       # chief
    wx::SCDI.Beam{T, SCDI.Ray{T}}               # waist x
    wy::SCDI.Beam{T, SCDI.Ray{T}}               # waist y
    dx::SCDI.Beam{T, SCDI.Ray{T}}               # div. x
    dy::SCDI.Beam{T, SCDI.Ray{T}}               # div. y
    parent::SCDI.Nullable{AstigmaticGaussianBeamlet{T}}
    children::Vector{AstigmaticGaussianBeamlet{T}}
end

function AstigmaticGaussianBeamlet(
    c::SCDI.Beam{T, SCDI.PolarizedRay{T}},
    wx::SCDI.Beam{T, SCDI.Ray{T}},
    wy::SCDI.Beam{T, SCDI.Ray{T}},
    dx::SCDI.Beam{T, SCDI.Ray{T}},
    dy::SCDI.Beam{T, SCDI.Ray{T}}
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
    s1 = [0,1,0] # SCDI.normal3d(direction)
    s2 = cross(direction, s1)
    # Divergence angle in rad
    θ = SCDI.divergence_angle(λ, w0, M2)
    # Waist rays
    wx = SCDI.Ray(position + s1 * w0, direction, λ)
    wy = SCDI.Ray(position + s2 * w0, direction, λ)
    # Divergence ray
    div_dir_x = normalize(direction + s1 * tan(θ))
    div_dir_y = normalize(direction + s2 * tan(θ))
    dx = SCDI.Ray(position, div_dir_x, λ)
    dy = SCDI.Ray(position, div_dir_y, λ)
    # Chief ray
    c = SCDI.PolarizedRay(position, direction, λ, E0)
    return AstigmaticGaussianBeamlet(SCDI.Beam(c), SCDI.Beam(wx), SCDI.Beam(wy), SCDI.Beam(dx), SCDI.Beam(dy))
end

SCDI._last_beam_intersection(agb::AstigmaticGaussianBeamlet) = SCDI.intersection(last(SCDI.rays(agb.c)))

SCDI.point_on_beam(agb::AstigmaticGaussianBeamlet, t::Real) = SCDI.point_on_beam(agb.c, t)

function SCDI.trace_system!(system::SCDI.AbstractSystem, agb::AstigmaticGaussianBeamlet; r_max::Int = 20)
    # placeholder solver
    SCDI.trace_system!(system, agb.c, r_max=r_max)
    SCDI.trace_system!(system, agb.wx, r_max=r_max)
    SCDI.trace_system!(system, agb.wy, r_max=r_max)
    SCDI.trace_system!(system, agb.dx, r_max=r_max)
    SCDI.trace_system!(system, agb.dy, r_max=r_max)
end

function SCDI.render_beam!(axis, agb::AstigmaticGaussianBeamlet; flen = 0.1)
    # placeholder render
    SCDI.render_beam!(axis, agb.c, flen = flen, color = :black)
    SCDI.render_beam!(axis, agb.wx, flen = flen, color = :green)
    SCDI.render_beam!(axis, agb.dx, flen = flen, color = :green)
    SCDI.render_beam!(axis, agb.wy, flen = flen, color = :red)
    SCDI.render_beam!(axis, agb.dy, flen = flen, color = :red)
end

function gauss_parameters(agb::AstigmaticGaussianBeamlet, z::Real)
    p0, i = SCDI.point_on_beam(agb, z)

    n = SCDI.refractive_index(agb.c.rays[i])
    λ = SCDI.wavelength(agb.c.rays[i])

    wx, ~, ~, ~, Hnx = SCDI.gauss_parameters(agb.c.rays[i], agb.wx.rays[i], agb.dx.rays[i], p0)
    wy, ~, ~, ~, Hny = SCDI.gauss_parameters(agb.c.rays[i], agb.wy.rays[i], agb.dy.rays[i], p0)

    return wx, wy, Hnx, Hny
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
