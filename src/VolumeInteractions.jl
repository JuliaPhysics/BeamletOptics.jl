"""
    current_medium(system::AbstractSystem, hint::Nullable{Hint})

Determines the optical medium through which the ray is currently propagating.
Returns the ambient medium of `system` if `hint === nothing`, or the medium of `hint.object`.
"""
@inline current_medium(system::AbstractSystem, ::Nothing) = ambient_medium(system)
@inline current_medium(::AbstractSystem, hint::Hint) = medium_from(object(hint))

"""
    propagate_volume(system::AbstractSystem, medium::AbstractMedium, ray::AbstractRay, hint, opl_accum)

Propagates `ray` through `medium` in `system` until the next surface intersection or infinity.
Constructs and attaches the corresponding `AbstractSegment` (e.g. `LineSegment`, `OpenSegment`, `SplineSegment`)
to the returned resolved ray.

Returns a 3-tuple `(resolved_ray, hit_intersection, hit_object)`.
"""
function propagate_volume end

@inline function propagate_volume(
        system::AbstractSystem,
        ::AbstractMedium,
        ray::AbstractRay{T},
        hint::Nullable{Hint},
        opl_accum::T
) where {T <: Real}
    res = trace_one(system, ray, hint)
    if res === nothing
        seg = OpenSegment(position(ray), direction(ray), opl_accum)
        resolved_ray = ray isa PolarizedRay ? PolarizedRay(ray, seg) : Ray(ray, seg)
        return (resolved_ray, nothing, nothing)
    else
        hit, obj = res
        seg_opl = opl_accum + length(hit) * real(refractive_index(ray))
        hit_geom = hit isa MultiIntersection ? hit.hit : hit
        seg = LineSegment(position(ray), direction(ray), hit_geom, seg_opl)
        resolved_ray = ray isa PolarizedRay ? PolarizedRay(ray, seg) : Ray(ray, seg)
        return (resolved_ray, hit, obj)
    end
end

"""
    interact3d(model::AbstractVolumeModel, medium::AbstractMedium, ray::AbstractRay)
    interact3d(model::AbstractVolumeModel, object::AbstractObject, ray::AbstractRay)

Universal volumetric interaction function evaluating continuous or bulk volume physics (e.g. GRIN ray tracing, participating media).
"""
function interact3d(
    model::AbstractVolumeModel,
    medium::AbstractMedium,
    ray::AbstractRay{T}
) where {T <: Real}
    error("Volume interaction not implemented for $(typeof(model)) in $(typeof(medium)).")
end

function interact3d(
    model::AbstractVolumeModel,
    object::AbstractObject,
    ray::AbstractRay{T}
) where {T <: Real}
    error("Volume interaction not implemented for $(typeof(model)) in $(typeof(object)).")
end
