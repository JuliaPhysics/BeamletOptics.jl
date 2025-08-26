"""
    PolarizationFilter <: AbstractJonesPolarizer

Represents a zero-thickness, ideal polarization filter. 
"""
struct PolarizationFilter{T, S <: AbstractShape{T}} <: AbstractJonesPolarizer{T,S}
    shape::S
    JMat::GlobalJonesBasis{T}
    cutoff::T
end

function PolarizationFilter(shape::S, J::GlobalJonesBasis{TJ}, cs) where {TS, S<:AbstractShape{TS}, TJ}
    return PolarizationFilter{TS,S}(shape, GlobalJonesBasis{TS}(J), TS(cs))
end

"""
    PolarizationFilter(size; cutoff_strength)

Spawns a thin, rectangular [`PolarizationFilter`](@ref) with the edge length as specified via the `size` in [m].
The filter is aligned with the global y-axis and transmits along the x-axis.
"""
function PolarizationFilter(size::Real; cutoff_strength=eps())
    shape = QuadraticFlatMesh(size)
    # Rotate normals against pos. y-axis
    zrotate3d!(shape, π)
    set_new_origin3d!(shape)
    return PolarizationFilter(shape, XZBasis(1, 0, 0, 0), cutoff_strength)
end

function interact3d(::AbstractSystem,
        polfilter::PolarizationFilter,
        ::Beam{T, R},
        ray::R) where {T <: Real, R <: PolarizedRay{T}}
    npos = position(ray) + length(ray) * direction(ray)
    ndir = direction(ray)

    E0 = _calculate_global_E0(polfilter, ray, ndir, polfilter.JMat)

    # Terminate blocked rays
    # FIXME this needs to be reworked in the future
    if norm(E0) ≈ polfilter.cutoff
        return nothing
    end

    return BeamInteraction{T, R}(nothing,
        PolarizedRay{T}(npos, ndir, nothing, wavelength(ray), refractive_index(ray), E0))
end
