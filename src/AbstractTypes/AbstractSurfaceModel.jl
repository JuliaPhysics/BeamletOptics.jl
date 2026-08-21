"""
    AbstractSurfaceModel

Abstract supertype governing the physical behavior at an optical boundary interface.
"""
abstract type AbstractSurfaceModel end

"""
    FresnelInterface <: AbstractSurfaceModel

A bare dielectric interface governed by Snell's law of refraction and Fresnel reflection equations.
"""
struct FresnelInterface <: AbstractSurfaceModel end

"""
    IdealMirror <: AbstractSurfaceModel

A perfect specular reflective surface with 100% reflectance across all wavelengths.
"""
struct IdealMirror <: AbstractSurfaceModel end

"""
    AbsorbingSurface <: AbstractSurfaceModel

A perfectly absorbing surface that terminates incoming rays (e.g. apertures, baffles, beam dumps).
"""
struct AbsorbingSurface <: AbstractSurfaceModel end

"""
    DetectorSurface{D} <: AbstractSurfaceModel

A surface model that records incoming light (power, phase, wavefront, or electric field distribution).
"""
struct DetectorSurface{D} <: AbstractSurfaceModel
    detector_data::D
end

DetectorSurface() = DetectorSurface(nothing)

"""
    CoatedSurface{C} <: AbstractSurfaceModel

A surface model representing a thin-film optical coating (AR, HR, partial beamsplitter, dichroic).
"""
struct CoatedSurface{C} <: AbstractSurfaceModel
    coating_model::C
end

"""
    GratingSurface{T<:Real, G} <: AbstractSurfaceModel

A diffraction grating surface characterized by a groove density, diffraction order, and ruling direction.
"""
struct GratingSurface{T <: Real, G} <: AbstractSurfaceModel
    groove_density::T
    order::Int
    groove_vector::G
end

GratingSurface(groove_density::Real, order::Int = 1) = 
    GratingSurface(Float64(groove_density), order, nothing)

surface_model(::AbstractObject) = FresnelInterface()
surface_model(::Nothing) = FresnelInterface()
surface_model(m::AbstractSurfaceModel) = m
