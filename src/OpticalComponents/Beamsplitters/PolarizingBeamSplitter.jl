
"""
    PolarizingBeamSplitter{T,S} <: AbstractBeamsplitter{T,S}

Ideal polarizing plate that separates an incoming `PolarizedRay` into
transmitted and reflected beams. The device is represented by a zero‑thickness
surface whose orientation sets the splitting axes:

- local `x`‑axis → transmitted (Ex) polarization component
- local `z`‑axis → reflected (Ez) polarization component

Rotate the underlying shape to align these axes with the desired global
polarization directions. Incoming rays are assumed to strike the plate from the
positive local `y` direction.
"""
struct PolarizingBeamSplitter{T, S <: AbstractShape{T}} <: AbstractBeamsplitter{T}
    shape::S
end

shape(pbs::PolarizingBeamSplitter) = pbs.shape

function PolarizingBeamSplitter(width::Real, height::Real)
    PolarizingBeamSplitter(RectangularFlatMesh(width, height))
end

@inline function _pbs_transmitted_beam(pbs::PolarizingBeamSplitter,
        ::Beam{T, R}, ray::R) where {T <: Real, R <: PolarizedRay{T}}
    pos = position(ray) + length(ray) * direction(ray)
    dir = direction(ray)
    J = XZBasis(one(T), zero(T), zero(T), zero(T))
    E0 = _calculate_global_E0(pbs, ray, dir, J)
    return Beam(PolarizedRay(pos, dir, wavelength(ray), E0))
end

@inline function _pbs_reflected_beam(pbs::PolarizingBeamSplitter,
        ::Beam{T, R}, ray::R) where {T <: Real, R <: PolarizedRay{T}}
    pos = position(ray) + length(ray) * direction(ray)
    normal = normal3d(intersection(ray))
    in_dir = direction(ray)
    out_dir = reflection3d(in_dir, normal)
    J = XZBasis(zero(T), zero(T), zero(T), -one(T))
    E0 = _calculate_global_E0(pbs, ray, out_dir, J)
    return Beam(PolarizedRay(pos, out_dir, wavelength(ray), E0))
end

function interact3d(::AbstractSystem, pbs::PolarizingBeamSplitter,
        beam::Beam{T, R}, ray::R) where {T <: Real, R <: PolarizedRay{T}}
    children!(
        beam, [_pbs_transmitted_beam(pbs, beam, ray),
            _pbs_reflected_beam(pbs, beam, ray)])
    return nothing
end

@inline function _pbs_transmitted_beam(
        pbs::PolarizingBeamSplitter, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    c = _pbs_transmitted_beam(pbs, agb.c, rays(agb.c)[ray_id])
    wxp = _beamsplitter_transmitted_beam(pbs, agb.wxp, rays(agb.wxp)[ray_id])
    wxm = _beamsplitter_transmitted_beam(pbs, agb.wxm, rays(agb.wxm)[ray_id])
    wyp = _beamsplitter_transmitted_beam(pbs, agb.wyp, rays(agb.wyp)[ray_id])
    wym = _beamsplitter_transmitted_beam(pbs, agb.wym, rays(agb.wym)[ray_id])
    dxp = _beamsplitter_transmitted_beam(pbs, agb.dxp, rays(agb.dxp)[ray_id])
    dxm = _beamsplitter_transmitted_beam(pbs, agb.dxm, rays(agb.dxm)[ray_id])
    dyp = _beamsplitter_transmitted_beam(pbs, agb.dyp, rays(agb.dyp)[ray_id])
    dym = _beamsplitter_transmitted_beam(pbs, agb.dym, rays(agb.dym)[ray_id])
    return AstigmaticGaussianBeamlet(c, wxp, wxm, wyp, wym, dxp, dxm, dyp, dym)
end

@inline function _pbs_reflected_beam(
        pbs::PolarizingBeamSplitter, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    c = _pbs_reflected_beam(pbs, agb.c, rays(agb.c)[ray_id])
    wxp = _beamsplitter_reflected_beam(pbs, agb.wxp, rays(agb.wxp)[ray_id])
    wxm = _beamsplitter_reflected_beam(pbs, agb.wxm, rays(agb.wxm)[ray_id])
    wyp = _beamsplitter_reflected_beam(pbs, agb.wyp, rays(agb.wyp)[ray_id])
    wym = _beamsplitter_reflected_beam(pbs, agb.wym, rays(agb.wym)[ray_id])
    dxp = _beamsplitter_reflected_beam(pbs, agb.dxp, rays(agb.dxp)[ray_id])
    dxm = _beamsplitter_reflected_beam(pbs, agb.dxm, rays(agb.dxm)[ray_id])
    dyp = _beamsplitter_reflected_beam(pbs, agb.dyp, rays(agb.dyp)[ray_id])
    dym = _beamsplitter_reflected_beam(pbs, agb.dym, rays(agb.dym)[ray_id])
    return AstigmaticGaussianBeamlet(c, wxp, wxm, wyp, wym, dxp, dxm, dyp, dym)
end

function interact3d(
        ::AbstractSystem, pbs::PolarizingBeamSplitter, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    t = _pbs_transmitted_beam(pbs, agb, ray_id)
    r = _pbs_reflected_beam(pbs, agb, ray_id)
    children!(agb, [t, r])
    return nothing
end
