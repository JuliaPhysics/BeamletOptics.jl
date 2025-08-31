
using StaticArrays

"""
    PolarizingBeamSplitter{T,S} <: AbstractBeamsplitter{T,S}

Ideal polarizing plate that splits the incoming `PolarizedRay` into two beams.
The component is modelled as a zero thickness surface described by a 2D shape.
The local `x`-axis defines the transmitted polarization while the local `z`-axis
defines the reflected polarization.
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
    fast = SVector{3,T}(orientation(shape(pbs))[:,1])
    fast -= dir * dot(fast, dir)
    fast = normalize(fast)
    pol = SVector{3,Complex{T}}(polarization(ray))
    Ex = dot(pol, fast)
    E0 = fast * Ex
    return Beam(PolarizedRay(pos, dir, wavelength(ray), E0))
end

@inline function _pbs_reflected_beam(pbs::PolarizingBeamSplitter,
        ::Beam{T,R}, ray::R) where {T<:Real,R<:PolarizedRay{T}}
    pos = position(ray) + length(ray) * direction(ray)
    normal = normal3d(intersection(ray))
    in_dir = direction(ray)
    out_dir = reflection3d(in_dir, normal)
    slow = SVector{3,T}(orientation(shape(pbs))[:,3])
    slow -= out_dir * dot(slow, out_dir)
    slow = normalize(slow)
    pol = SVector{3,Complex{T}}(polarization(ray))
    Ez = dot(pol, slow)
    E0 = slow * Ez
    return Beam(PolarizedRay(pos, out_dir, wavelength(ray), E0))
end

function interact3d(::AbstractSystem, pbs::PolarizingBeamSplitter,
        beam::Beam{T,R}, ray::R) where {T<:Real,R<:PolarizedRay{T}}
    children!(beam, [_pbs_transmitted_beam(pbs, beam, ray),
                    _pbs_reflected_beam(pbs, beam, ray)])
    return nothing
end

