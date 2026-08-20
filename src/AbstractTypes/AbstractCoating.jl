"""
    AbstractCoating <: AbstractObject

Abstract supertype for optical boundary containers that pair geometry with an [`AbstractSurfaceModel`](@ref).
"""
abstract type AbstractCoating{T} <: AbstractObject{T} end

"""
    Coating{T, S<:AbstractShape{T}, M<:AbstractSurfaceModel} <: AbstractCoating{T}

A container representing an optical boundary interface.
- `shape`: The boundary geometry (e.g. partial-volume SDF, mesh patch, or NURBS face).
- `model`: The boundary interaction physics (e.g. `ARCoating`, `FresnelInterface`, `IdealMirror`).
"""
struct Coating{T <: Real, S <: AbstractShape{T}, M <: AbstractSurfaceModel} <: AbstractCoating{T}
    shape::S
    model::M
end

Coating(shape::AbstractShape{T}, model::AbstractSurfaceModel) where {T} = Coating{T, typeof(shape), typeof(model)}(shape, model)

shape(c::Coating) = c.shape
surface_model(c::Coating) = c.model

Base.position(c::Coating) = position(c.shape)
orientation(c::Coating) = orientation(c.shape)
