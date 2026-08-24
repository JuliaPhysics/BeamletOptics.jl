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
