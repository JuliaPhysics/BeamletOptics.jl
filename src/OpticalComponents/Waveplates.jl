"""
    Waveplate{T,S<:AbstractShape{T}} <: AbstractObject{T,S}

Ideal zero-thickness birefringent plate introducing a relative phase delay
`δ` between its fast and slow axes.  The local x-axis of the `shape`
corresponds to the fast axis, while the normal of the plate coincides with the
local y-axis.  Waveplates may be created with rectangular or circular shapes via
the provided constructors.
"""
struct Waveplate{T<:Real,S<:AbstractShape{T}} <: AbstractObject{T,S}
    shape::S
    δ::T
end

Waveplate(shape::S, δ::Real) where {T<:Real,S<:AbstractShape{T}} =
    Waveplate{T,S}(shape, T(δ))

shape(wp::Waveplate) = wp.shape
retardance(wp::Waveplate) = wp.δ

"""
    Waveplate(width, height, δ)

Construct a rectangular waveplate with edge lengths `width` and `height`
introducing a phase delay `δ`.
"""
Waveplate(width::Real, height::Real, δ::Real) =
    Waveplate(RectangularFlatMesh(width, height), δ)

"""
    Waveplate(width, δ)

Convenience constructor for a quadratic waveplate.
"""
Waveplate(width::Real, δ::Real) = Waveplate(width, width, δ)

"""
    RoundWaveplate(diameter, δ)

Construct a circular waveplate with `diameter` introducing a phase delay `δ`.
"""
RoundWaveplate(diameter::Real, δ::Real) =
    Waveplate(CircularFlatMesh(diameter/2), δ)

QuarterWaveplate(width::Real) = Waveplate(width, π/2)
HalfWaveplate(width::Real) = Waveplate(width, π)
RoundQuarterWaveplate(d::Real) = RoundWaveplate(d, π/2)
RoundHalfWaveplate(d::Real) = RoundWaveplate(d, π)

"""
    _waveplate_E0(E0, dir, fast, normal, δ)

Internal helper calculating the transformed electric field vector after passing
through a waveplate. `E0` denotes the incoming field vector, `dir` the
propagation direction, `fast` the fast axis and `normal` the plate normal.
`δ` is the retardance in radians.
"""
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
