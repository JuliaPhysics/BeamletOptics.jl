"""
    Waveplate{T,S<:AbstractShape{T}} <: AbstractObject{T,S}

Represents an ideal thin waveplate introducing a relative phase delay `δ` between
its fast and slow axes. The fast axis is aligned with the local x-axis of the
`shape`, the plate normal with the local y-axis.
"""
struct Waveplate{T<:Real,S<:AbstractShape{T}} <: AbstractObject{T,S}
    shape::S
    δ::T
end

Waveplate(shape::S, δ::Real) where {T<:Real,S<:AbstractShape{T}} =
    Waveplate{T,S}(shape, T(δ))

shape(wp::Waveplate) = wp.shape
retardance(wp::Waveplate) = wp.δ

Waveplate(width::Real, height::Real, δ::Real) =
    Waveplate(RectangularFlatMesh(width, height), δ)
Waveplate(width::Real, δ::Real) = Waveplate(width, width, δ)

QuarterWaveplate(width::Real) = Waveplate(width, π/2)
HalfWaveplate(width::Real) = Waveplate(width, π)

function _waveplate_E0(E0, dir, fast, normal, δ)
    fast = normalize(fast - dot(fast, normal)*normal)
    slow = normalize(cross(normal, fast))
    Ef = dot(E0, fast)
    Es = dot(E0, slow)
    Et = E0 - Ef*fast - Es*slow
    return Ef*fast + Es*exp(im*δ)*slow + Et
end

function interact3d(::AbstractSystem, wp::Waveplate, ::Beam{T,R}, ray::R) where {T<:Real,R<:Ray{T}}
    pos = position(ray) + length(ray)*direction(ray)
    dir = direction(ray)
    return BeamInteraction{T,R}(nothing, Ray(pos, dir, nothing, wavelength(ray), refractive_index(ray)))
end

function interact3d(::AbstractSystem, wp::Waveplate, ::Beam{T,R}, ray::R) where {T<:Real,R<:PolarizedRay{T}}
    pos = position(ray) + length(ray)*direction(ray)
    dir = direction(ray)
    n = orientation(shape(wp))[:,2]
    fast = orientation(shape(wp))[:,1]
    E0 = _waveplate_E0(polarization(ray), dir, fast, n, retardance(wp))
    return BeamInteraction{T,R}(nothing,
        PolarizedRay(pos, dir, wavelength(ray), E0))
end
