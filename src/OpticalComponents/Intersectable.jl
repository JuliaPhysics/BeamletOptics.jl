"""
    IntersectableObject

A passive [`AbstractObject`](@ref) which can be hit by a beam. In this case, the object acts like a hard target and blocks the beam path.

# Fields

- `shape`: an [`AbstractShape`](@ref)
"""
struct IntersectableObject{T, S <: AbstractShape{T}} <: AbstractObject{T, S}
    shape::S
end

set_new_origin3d!(d::IntersectableObject) = set_new_origin3d!(d.shape)
interact3d(::AbstractSystem, ::IntersectableObject, ::AbstractBeam, ::AbstractRay) = nothing