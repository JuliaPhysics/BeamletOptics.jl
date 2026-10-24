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
    interact3d(system::AbstractSystem, polfilter::PolarizationFilter, beam::Beam, ray::PolarizedRay)

Applies the Jones matrix of the filter to the field of the `ray`. The filter may lie on the boundary
between two media, e.g. a film cemented between two glass plates, see [`LinearPolarizer`](@ref): the
filtered field is then refracted into the medium behind the filter with the
[`fresnel_coefficients`](@ref) of the [`transition`](@ref). Beyond the critical angle of these media,
the `ray` is totally reflected before it reaches the filter. A ray whose field falls to or below the
`cutoff` of the filter is stopped.
"""
function interact3d(system::AbstractSystem,
        polfilter::PolarizationFilter,
        beam::Beam{T, R},
        ray::R) where {T <: Real, R <: PolarizedRay{T}}
    t = transition(system, ray)
    is_totally_reflected(ray, t) && return cross_boundary(t, beam, ray)
    E = _calculate_global_E0(polfilter, ray, direction(ray), polfilter.JMat)
    # Terminate blocked rays
    if norm(E) ≤ polfilter.cutoff
        return nothing
    end
    n1, n2 = medium_in(t), medium_out(t)
    new_dir, _ = refraction3d(direction(ray), normal3d(t), n1, n2)
    _, _, ts, tp = fresnel_coefficients(angle3d(direction(ray), -normal3d(t)), n2 / n1)
    P = _calculate_global_E0(direction(ray), new_dir, normal3d(intersection(ray)), SPBasis(ts, 0, 0, tp))
    E0 = _transverse(P * E, new_dir)
    return BeamInteraction{T, R}(hint_out(t),
        PolarizedRay{T}(hit_point(ray), new_dir, nothing, wavelength(ray), n2, E0))
end

# unpolarized rays (e.g. beamlet auxiliary rays) pass through the ideal filter, refracted at the
# boundary between the media on both sides of the filter
function interact3d(
        system::AbstractSystem,
        ::PolarizationFilter,
        beam::Beam{T, R},
        ray::R
    ) where {T <: Real, R <: Ray{T}}
    return cross_boundary(system, beam, ray)
end
