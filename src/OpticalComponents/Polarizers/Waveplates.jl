using StaticArrays

"""
    Waveplate{T,S} <: AbstractJonesPolarizer{T,S}

The wave plate is modeled as a zero‑thickness element described by a 2D
`AbstractShape`. The surface normal points along the local `y`‑axis, and the
local `x`‑axis marks the fast axis of the retarder. Rotating the object therefore
rotates the fast axis in global coordinates. The plate does not deflect rays; it
merely imposes a phase delay between field components parallel to the fast
(`x`) and slow (`z`) axes.

# Fields

- `shape`: 2D shape describing the physical aperture of the plate
- `retardance`: phase delay between fast and slow axis in radians
"""
struct Waveplate{T,S<:AbstractShape{T}} <: AbstractJonesPolarizer{T,S}
    shape::S
    retardance::T
end

Waveplate(shape::S, retardance::Real) where {T,S<:AbstractShape{T}} =
    Waveplate{T,S}(shape, T(retardance))

shape(wp::Waveplate) = wp.shape

"""
    Waveplate(width, height, retardance)
    Waveplate(diameter, retardance)

Construct a wave plate centred at the origin and aligned to the `y`-axis. The
fast axis coincides with the local `x`-axis of the shape. Providing `width` and
`height` creates a rectangular plate while supplying a single `diameter`
argument creates a circular plate.
"""
Waveplate(width::Real, height::Real, retardance::Real) =
    Waveplate(RectangularFlatMesh(width, height), retardance)

Waveplate(diameter::Real, retardance::Real) =
    Waveplate(CircularFlatMesh(diameter/2), retardance)

"""
    HalfWaveplate(width, height)
    HalfWaveplate(diameter)

Convenience constructor for a `Waveplate` with retardance `π`.
"""
HalfWaveplate(width::Real, height::Real) = Waveplate(width, height, π)
HalfWaveplate(diameter::Real) = Waveplate(diameter, π)

"""
    QuarterWaveplate(width, height)
    QuarterWaveplate(diameter)

Convenience constructor for a `Waveplate` with retardance `π/2`.
"""
QuarterWaveplate(width::Real, height::Real) = Waveplate(width, height, π/2)
QuarterWaveplate(diameter::Real) = Waveplate(diameter, π/2)

"""
    interact3d(AbstractSystem, Waveplate, Beam, PolarizedRay)

Applies the Jones matrix of the wave plate to the polarization state of the
incoming ray. The plate does not alter the direction of propagation.
"""
function interact3d(::AbstractSystem, wp::Waveplate, ::Beam{T,R},
        ray::R) where {T<:Real, R<:PolarizedRay{T}}
    pos = position(ray) + length(ray) * direction(ray)
    dir = direction(ray)
    J = XZBasis(one(Complex{T}), zero(Complex{T}), zero(Complex{T}), exp(im * wp.retardance))
    E0 = _calculate_global_E0(wp, ray, dir, J)
    new_ray = PolarizedRay(pos, dir, wavelength(ray), E0)
    refractive_index!(new_ray, refractive_index(ray))
    return BeamInteraction{T,R}(nothing, new_ray)
end

