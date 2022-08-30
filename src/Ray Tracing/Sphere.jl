"""
AbstractSphere <: AbstractObject

A generic type for an object whose volume can be described by a sphere. Must have a field `sphere` of type `Sphere`. See also `Sphere{T}`.
"""
abstract type AbstractSphere{T} <: AbstractObject{T} end

"Enforces that `sphere`-like objects have to have the field `sphere` or implement `sphere()`."
sphere(object::AbstractSphere) = object.sphere
sphere!(object::AbstractSphere, sphere) = nothing

mutable struct Sphere{T} <: AbstractSphere{T}
    pos::Vector{T}
    radius::T
end

orientation(sphere::Sphere{T}) where T = Matrix{T}(I, 3, 3)
orientation!(sphere::Sphere, dir) = nothing 