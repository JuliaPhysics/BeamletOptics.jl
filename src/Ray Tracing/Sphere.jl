abstract type AbstractSphere <: AbstractEntity end

struct Sphere{T} <: AbstractSphere # is this inheritance correct? see also Mesh{T}
    pos::Vector{T}
    radius::T
end

"""
    sphere(object::AbstractSphere)

Enforces that objects with meshes have to have the field mesh or implement `mesh`.
"""
sphere(object::AbstractSphere) = object.sphere

translate3d!(object::AbstractSphere, offset::Vector) = nothing
scale3d!(object::AbstractSphere, scale) = nothing
# Mesh intersection
intersect3d(object::AbstractSphere, ray::Ray) = nothing
# Utils
orthogonal3d(object::AbstractSphere, fID::Int) = nothing
reset_translation3d!(object::AbstractSphere) = nothing
reset_rotation3d!(object::AbstractSphere) = nothing
set_new_origin3d!(object::AbstractSphere) = nothing
# Render
render_object!(axis, object::AbstractSphere) = nothing
render_object_normals!(axis, object::AbstractSphere) = nothing
