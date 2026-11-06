"""
    PolarizationFilter <: AbstractJonesPolarizer

Represents a zero-thickness, ideal polarization filter. 
"""
struct PolarizationFilter{T, S <: AbstractShape{T}} <: AbstractJonesPolarizer{T}
    shape::S
    JMat::GlobalJonesBasis{T}
    cutoff::T
end

function PolarizationFilter(shape::S, J::GlobalJonesBasis{TJ}, cs) where {TS, S<:AbstractShape{TS}, TJ}
    return PolarizationFilter{TS,S}(shape, GlobalJonesBasis{TS}(J), TS(cs))
end

"""
    PolarizationFilter(edge_length; cutoff_strength)

Spawns a thin, rectangular [`PolarizationFilter`](@ref). The `edge_length` has to be specified in [m].
The filter is aligned with the global y-axis and transmits along the x-axis, while blocking polarization components
along the global z-axis.

A [`PolarizedRay`](@ref) whose polarization vector has a norm of at most `cutoff_strength` behind the filter
is blocked, i.e. its beam ends here. Rays without polarization pass unchanged.
"""
function PolarizationFilter(edge_length::Real; cutoff_strength=eps())
    shape = QuadraticFlatMesh(edge_length)
    # Rotate normals against pos. y-axis
    zrotate3d!(shape, π)
    set_new_origin3d!(shape)
    return PolarizationFilter(shape, XZBasis(1, 0, 0, 0), cutoff_strength)
end

"""
    RoundPolarizationFilter(diameter; cutoff_strength)

Spawns a thin, round [`PolarizationFilter`](@ref) with the given `diameter` in [m].
The filter is centered at the origin, aligned with the global y-axis and transmits along the x-axis,
while blocking polarization components along the global z-axis.
"""
function RoundPolarizationFilter(diameter::Real; cutoff_strength=eps())
    shape = CircularFlatMesh(diameter / 2)
    return PolarizationFilter(shape, XZBasis(1, 0, 0, 0), cutoff_strength)
end

"""
    transmission_axis(pf::PolarizationFilter)

Returns the unit vector (in global coordinates) along which [`PolarizationFilter`](@ref) `pf` transmits
polarization, derived from its Jones matrix so that it stays correct for any filter orientation or custom
`GlobalJonesBasis`. The sign of the returned vector is arbitrary, since it represents an axis rather than
a direction.
"""
function transmission_axis(pf::PolarizationFilter)
    R = orientation(pf)
    n = R[:, 2]
    Q = I - n * n'
    J = static_data(pf.JMat)
    M = Q * R * J * R' * Q
    F = svd(M)
    return normalize(real(F.V[:, 1]))
end

"""
    JonesFilm{T, M} <: AbstractSurfaceModel

An ideal, infinitesimally thin polarizing film with the Jones matrix `J` (3×3, in global coordinates).
It filters the field of a [`PolarizedRay`](@ref) with `J` projected into the plane transverse to the
ray, see [`AbstractJonesPolarizer`](@ref), and stops the ray if the field falls to or below `cutoff`.
Rays without a field, e.g. the auxiliary rays of a beamlet, pass unfiltered.

The film acts between the media of the boundary, e.g. as a film cemented between two glass plates, see
[`LinearPolarizer`](@ref): the filtered field is refracted into the medium behind the film with the
[`fresnel_coefficients`](@ref). Beyond the critical angle of these media, the ray is totally reflected
before it reaches the film.

A [`PolarizationFilter`](@ref) carries this film in its current orientation, see
[`surface_model`](@ref).
"""
struct JonesFilm{T, M <: AbstractMatrix} <: AbstractSurfaceModel
    J::M
    cutoff::T
end

function surface_model(pf::PolarizationFilter)
    R = orientation(pf)
    return JonesFilm(R * static_data(pf.JMat) * transpose(R), pf.cutoff)
end

# The filtered field is transmitted, a blocked ray leaves no transmitted channel; beyond the critical
# angle the ray is totally reflected before it reaches the film
function channels(film::JonesFilm, ::Transition, ray::PolarizedRay, cs)
    transmitted, reflected = cs
    iszero(amplitude(transmitted)) && return cs
    # Jones matrix projected into the plane transverse to the ray
    d = direction(ray)
    Q = I - d * transpose(d)
    J = Q * film.J * Q
    norm(J * polarization(ray)) ≤ film.cutoff && return (transfer(transmitted, zero(transmitted.P), 0), reflected)
    return (transfer(transmitted, transmitted.P * J, amplitude(transmitted)), reflected)
end

# Rays without a field pass unfiltered
channels(::JonesFilm, ::Transition, ::AbstractRay, cs) = cs
