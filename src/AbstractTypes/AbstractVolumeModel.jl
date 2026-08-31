"""
    AbstractVolumeModel

Abstract supertype governing continuous or volumetric optical interactions (e.g. GRIN media, participating media, volume scattering).
"""
abstract type AbstractVolumeModel end

"""
    VolumeInteraction <: AbstractVolumeModel

Stub type representing a volumetric interaction model.
"""
struct VolumeInteraction <: AbstractVolumeModel end
