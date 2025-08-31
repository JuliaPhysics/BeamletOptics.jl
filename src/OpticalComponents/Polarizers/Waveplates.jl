using StaticArrays

"""
    Waveplate{T,S} <: AbstractObject{T,S}

Thin retardation plate acting only on the polarization state of a `PolarizedRay`.
The wave plate is modelled as a zero-thickness element described by a 2D
`AbstractShape`. The local `x`-axis of the shape represents the fast axis of the
plate. Rotating the object therefore rotates the fast axis accordingly.

# Fields

- `shape`: 2D shape describing the physical aperture of the plate
- `retardance`: phase delay between fast and slow axis in radians
"""
struct Waveplate{T,S<:AbstractShape{T}} <: AbstractObject{T,S}
    shape::S
    retardance::T
end

Waveplate(shape::S, retardance::Real) where {T,S<:AbstractShape{T}} =
    Waveplate{T,S}(shape, T(retardance))

shape(wp::Waveplate) = wp.shape

"""
    Waveplate(width, height, retardance)

Construct a rectangular wave plate centred at the origin and aligned to the
`y`-axis. The fast axis coincides with the local `x`-axis of the shape.
"""
Waveplate(width::Real, height::Real, retardance::Real) =
    Waveplate(RectangularFlatMesh(width, height), retardance)

"""
    HalfWaveplate(width, height)

Convenience constructor for a `Waveplate` with retardance `π`.
"""
HalfWaveplate(width::Real, height::Real) = Waveplate(width, height, π)

"""
    QuarterWaveplate(width, height)

Convenience constructor for a `Waveplate` with retardance `π/2`.
"""
QuarterWaveplate(width::Real, height::Real) = Waveplate(width, height, π/2)

"""
    interact3d(AbstractSystem, Waveplate, Beam, Ray)

Non‑polarized rays pass through a `Waveplate` without modification.
"""
function interact3d(::AbstractSystem, wp::Waveplate, ::Beam{T,R},
        ray::R) where {T<:Real, R<:Ray{T}}
    pos = position(ray) + length(ray) * direction(ray)
    return BeamInteraction{T,R}(nothing,
        Ray{T}(pos, direction(ray), nothing, wavelength(ray), refractive_index(ray)))
end

"""
    interact3d(AbstractSystem, Waveplate, Beam, PolarizedRay)

Applies the Jones matrix of the wave plate to the polarization state of the
incoming ray. The plate does not alter the direction of propagation.
"""
function interact3d(::AbstractSystem, wp::Waveplate, ::Beam{T,R},
        ray::R) where {T<:Real, R<:PolarizedRay{T}}
    pos = position(ray) + length(ray) * direction(ray)
    dir_vec = SVector{3,T}(direction(ray))
    # Fast axis is local x-axis of the shape, slow axis is local z-axis
    fast = normalize(SVector{3,T}(orientation(shape(wp))[:,1]))
    slow = normalize(SVector{3,T}(orientation(shape(wp))[:,3]))
    O_in = SMatrix{3,3,T}(vcat(fast', slow', dir_vec'))
    O_out = SMatrix{3,3,T}(hcat(fast, slow, dir_vec))
    J = @SMatrix [one(Complex{T}) 0 0;
                   0 exp(im*wp.retardance) 0;
                   0 0 one(Complex{T})]
    E_local = O_in * polarization(ray)
    E_local = J * E_local
    E0 = O_out * E_local
    new_ray = PolarizedRay(pos, direction(ray), wavelength(ray), E0)
    refractive_index!(new_ray, refractive_index(ray))
    return BeamInteraction{T,R}(nothing, new_ray)
end

