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

"""
    ray_hits(pd::Detector) -> (points, directions)

Where each ray crossed the [`Detector`](@ref) `pd` and in which direction it went, in the
coordinates of the detector. `points[i]` is the local `(x, z)` of hit `i` in \\[m\\], as in
[`spot_diagram`](@ref), and `directions[i]` its unit direction in the right-handed frame
`(x̂, n̂, ẑ)`, with the local axes `x̂`, `ẑ` and the normal `n̂ = -direction(pd)`, along which
beams arrive on the front side. For a detector in its initial orientation this frame is the
global one, so a beam along the optical axis has the direction `(0, 1, 0)`. A hit on the back
side, which the detector also records, has `d[2] < 0`. The angle of a ray to the detector normal,
in `[0, π/2]` for either side, is `atan(hypot(d[1], d[3]), abs(d[2]))` (accurate also near `0`,
unlike `acos`), its distance from the detector center `norm(p)`, e.g. to judge how far a beam is
from collimated.

For [`GaussianBeamlet`](@ref)s and [`AstigmaticGaussianBeamlet`](@ref)s, each hit gives the
point and direction of the chief ray of the beamlet, whereas [`spot_diagram`](@ref) gives its
waist contour. The pairs follow the order of the stored hits, which is not deterministic for
beam groups. A detector that no beam hit returns two empty vectors.
"""
ray_hits(pd::Detector) = ray_hits(pd, hits(pd))

ray_hits(::Detector{T}, ::Nothing) where T = (Point2{T}[], Point3{T}[])

function ray_hits(pd::Detector{T}, hits::Vector{<:AbstractDetectorHit}) where T
    x̂, ẑ = local_axes(pd)
    # Downstream normal
    n̂ = -direction(pd)
    points = calc_local_pos(position(pd), x̂, ẑ, Point3{T}.(hit_point.(hits)))
    directions = [Point3{T}(dot(direction(hit), x̂), dot(direction(hit), n̂),
                      dot(direction(hit), ẑ)) for hit in hits]
    return points, directions
end