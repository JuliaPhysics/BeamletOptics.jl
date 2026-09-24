"""
    NonInteractableObject

A passive [`AbstractObject`](@ref) which does not interact with the ray tracing simulation but can be moved via the kinematic API.

# Fields

- `shape`: an [`AbstractShape`](@ref)

!!! info "Usage"
    This type is intended mainly for visualization purposes, e.g. kinematic mount [`Mesh`](@ref)s, or similar applications.
    In essence, this object behaves fully transparent. The `intersect3d` and `interact3d` methods default to `nothing`.
"""
struct NonInteractableObject{T, S <: AbstractShape{T}} <: AbstractObject{T}
    shape::S
end

set_new_origin3d!(d::NonInteractableObject) = set_new_origin3d!(d.shape)
intersect3d(::NonInteractableObject, ::AbstractRay) = nothing
interact3d(::AbstractSystem, ::NonInteractableObject, ::AbstractBeam, ::AbstractRay) = nothing

"""
    MeshDummy(loadpath::String)

Creates a [`NonInteractableObject`](@ref) with a [`Mesh`](@ref) loaded from the specified file path.
Useful for rendering background objects or geometry that does not interact with rays.
"""
MeshDummy(loadpath::String) = NonInteractableObject(Mesh(load(loadpath)))

"""
    KM100CPMount()

Returns a [`MeshDummy`](@ref) of a Thorlabs [KM100CP/M](https://www.thorlabs.com/thorproduct.cfm?partnumber=KM100CP/M)
kinematic mount for Ø1" optics on a post. The origin lies at the center of a mounted Ø1" mirror, i.e. a
[`RoundPlanoMirror`](@ref) spawned at the origin sits in the mount, and the post base is 81.8 mm below it.

The mesh is loaded from the documentation assets that ship with the package. Combine it with a mirror via
[`ObjectGroup`](@ref) to move both together.
"""
function KM100CPMount()
    mount = MeshDummy(joinpath(pkgdir(@__MODULE__), "docs", "src", "assets", "mirror_renders", "Mirror_Post.stl"))
    translate_to3d!(mount, [0, 0, -5.68e-2])
    set_new_origin3d!(mount)
    return mount
end