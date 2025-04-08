"""
    NonInteractableObject

A passive [`AbstractObject`](@ref) which does not interact with the ray tracing simulation but can be moved via the kinematic API.

# Fields

- `shape`: an [`AbstractShape`](@ref)

!!! info "Usage"
    This type is intended mainly for visualization purposes, e.g. kinematic mount [`Mesh`](@ref)s, or similar applications.
    In essence, this object behaves fully transparent. The `intersect3d` and `interact3d` methods default to `nothing`.
"""
struct NonInteractableObject{T, S <: AbstractShape{T}} <: AbstractObject{T, S}
    shape::S
end

set_new_origin3d!(d::NonInteractableObject) = set_new_origin3d!(d.shape)
intersect3d(::NonInteractableObject, ::AbstractRay) = nothing
interact3d(::AbstractSystem, ::NonInteractableObject, ::AbstractBeam, ::AbstractRay) = nothing

MeshDummy(loadpath::String) = NonInteractableObject(Mesh(load(loadpath)))

render_object!(axis, dummy::NonInteractableObject; kwargs...) = render_dummy_mesh!(axis, dummy; kwargs...)

render_dummy_mesh!(::Any, ::NonInteractableObject; kwargs...) = nothing
