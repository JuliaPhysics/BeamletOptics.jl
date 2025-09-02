
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
struct PolarizingBeamSplitter{T,S<:AbstractShape{T}} <: AbstractBeamsplitter{T,S}
    shape::S
end

shape(pbs::PolarizingBeamSplitter) = pbs.shape

PolarizingBeamSplitter(width::Real, height::Real) =
    PolarizingBeamSplitter(RectangularFlatMesh(width, height))

@inline function _pbs_transmitted_beam(pbs::PolarizingBeamSplitter,
        ::Beam{T,R}, ray::R) where {T<:Real,R<:PolarizedRay{T}}
    pos = position(ray) + length(ray) * direction(ray)
    dir = direction(ray)
    J = XZBasis(one(T), zero(T), zero(T), zero(T))
    E0 = _calculate_global_E0(pbs, ray, dir, J)
    return Beam(PolarizedRay(pos, dir, wavelength(ray), E0))
end

@inline function _pbs_reflected_beam(pbs::PolarizingBeamSplitter,
        ::Beam{T,R}, ray::R) where {T<:Real,R<:PolarizedRay{T}}
    pos = position(ray) + length(ray) * direction(ray)
    normal = normal3d(intersection(ray))
    in_dir = direction(ray)
    out_dir = reflection3d(in_dir, normal)
    J = XZBasis(zero(T), zero(T), zero(T), one(T))    
    E0 = _calculate_global_E0(pbs, ray, out_dir, J)
    return Beam(PolarizedRay(pos, out_dir, wavelength(ray), E0))
end

function interact3d(::AbstractSystem, pbs::PolarizingBeamSplitter,
        beam::Beam{T,R}, ray::R) where {T<:Real,R<:PolarizedRay{T}}
    children!(beam, [_pbs_transmitted_beam(pbs, beam, ray),
                    _pbs_reflected_beam(pbs, beam, ray)])
    return nothing
end
