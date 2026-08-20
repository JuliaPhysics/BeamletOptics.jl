"""
    AbstractCoating <: AbstractObject

Abstract supertype for optical boundary containers that pair geometry with an [`AbstractSurfaceModel`](@ref).
"""
abstract type AbstractCoating{T} <: AbstractObject{T} end

