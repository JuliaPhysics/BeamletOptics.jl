"""
    spot_diagram(d::Detector)

Returns an array of 2D points in local detector coordinates that represent the points of intersection for incoming beams.

# Beams

For [`Beam`](@ref)s, the point of intersection on the screen surface is stored.

# Beamlets

For [`GaussianBeamlet`](@ref)s, the projected 1/e² waist is returned. The number of points and radius can be adjusted via the 
`num_spots` and `crop_factor` keyword arguments.

A detector that no beam hit returns no points, an empty vector.
"""
spot_diagram(d::Detector) = spot_diagram(d, hits(d))

spot_diagram(::Detector{T}, ::Nothing; kwargs...) where T = Point2{T}[]

function spot_diagram(pd::Detector, ::Vector{<:AbstractRayHit{T}}; kwargs...) where T
    return calc_local_pos(pd)
end

function spot_diagram(pd::Detector, ::Vector{B}; kwargs...) where B<:AbstractBeamletHit
    return calc_local_pos(pd; kwargs...)
end