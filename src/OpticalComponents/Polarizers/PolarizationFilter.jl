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

function interact3d(::AbstractSystem,
        polfilter::PolarizationFilter,
        beam::Beam{T},
        ray::PolarizedRay{T}) where {T <: Real}
    int = isempty(beam.segments) ? intersect3d(shape(polfilter), ray) : intersection(last(beam.segments))
    isnothing(int) && return nothing
    npos = position(int)
    ndir = direction(ray)

    E0 = _calculate_global_E0(polfilter, ray, ndir, polfilter.JMat)

    # Terminate blocked rays
    if norm(E0) ≈ polfilter.cutoff || norm(E0) < polfilter.cutoff
        return nothing
    end

    return BeamInteraction(nothing,
        PolarizedRay(
            npos,
            ndir,
            wavelength(ray),
            refractive_index(ray),
            E0
        ))
end

function interact3d(system::AbstractSystem,
        polfilter::PolarizationFilter,
        agb::AstigmaticGaussianBeamlet{T},
        id::Int) where {T <: Real}
    # Chief ray interaction (handles polarization change)
    i_c = interact3d(system, polfilter, agb.c, rays(agb.c)[id])
    isnothing(i_c) && return nothing
    
    # Auxiliary rays only undergo geometric interaction with the filter shape.
    # They hit at their own transverse locations, preserving beam width/divergence.
    aux_ints = map(b -> begin
        interact3d(system, polfilter.shape, b, rays(b)[id])
    end, _aux_beams(agb))
    
    if any(isnothing, aux_ints)
        return nothing
    end
    
    return AstigmaticGaussianBeamletInteraction{T}(i_c, aux_ints...)
end
