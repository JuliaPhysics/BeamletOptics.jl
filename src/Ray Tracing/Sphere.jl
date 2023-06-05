"""
AbstractSphere <: AbstractObject

A generic type for an object whose volume can be described by a sphere. Must have a field `sphere` of type `Sphere`. See also `Sphere{T}`.
"""
abstract type AbstractSphere{T} <: AbstractShape{T} end

"Enforces that `sphere`-like objects have to have the field `sphere` or implement `sphere()`."
sphere(object::AbstractSphere) = object.sphere
sphere!(::AbstractSphere, ::Any) = nothing

mutable struct Sphere{T} <: AbstractSphere{T}
    id::UUID
    pos::Vector{T}
    radius::T
end

orientation(::Sphere{T}) where T = Matrix{T}(I, 3, 3)
orientation!(::Sphere, ::Any) = nothing 