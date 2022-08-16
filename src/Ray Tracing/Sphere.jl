abstract type AbstractSphere <: AbstractEntity end

mutable struct Sphere{T} <: AbstractSphere # is this inheritance correct? see also Mesh{T}
    pos::Vector{T}
    radius::T
end

"""
    translate3d!(sphere::Sphere, offset::Vector)

Mutating function that translates the center of origin of a `sphere` by the `offset` vector.
"""
function translate3d!(sphere::Sphere, offset::Vector)
    sphere.pos .+= offset
    return nothing
end

"""
    scale3d!(sphere::Sphere, scale)

Resizes the radius of a `sphere` by the `scale` factor.
"""
function scale3d!(sphere::Sphere, scale)
    sphere.radius *= scale
    return nothing
end

"""
    reset_translation3d!(sphere::Sphere)

Resets all previous translations and returns the `sphere` back to the global origin.
"""
function reset_translation3d!(sphere::Sphere)
    sphere.pos = zeros(eltype(mesh.pos), 3)
    return nothing
end

"""
    sphere(object::AbstractSphere)

Enforces that sphere-like objects have to have the field sphere or implement `sphere`.
"""
sphere(object::AbstractSphere) = object.sphere

translate3d!(object::AbstractSphere, offset::Vector) = translate3d!(sphere(object), offset)
scale3d!(object::AbstractSphere, scale) = scale3d!(sphere(object), scale)
# Mesh intersection
intersect3d(object::AbstractSphere, ray::Ray) = intersect3d(sphere(object), ray)
# Mesh interaction
interact(object::AbstractSphere, beam::Beam, ~) = false
# Utils
orthogonal3d(object::AbstractSphere, fID::Int) = nothing
reset_translation3d!(object::AbstractSphere) = reset_translation3d!(sphere(object))
reset_rotation3d!(object::AbstractSphere) = nothing
set_new_origin3d!(object::AbstractSphere) = nothing
# Render
render_object!(axis, object::AbstractSphere) = nothing
render_object_normals!(axis, object::AbstractSphere) = nothing
