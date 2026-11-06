"""
    AbstractJonesPolarizer <: AbstractObject
    
Represents infinitesimally thin components that change the polarization state of incoming [`PolarizedRay`](@ref)s via global Jones matrix calculus.
Rather than using the generic Yun ray tracing scheme as referred to in the `PolarizedRay` docs, this element interacts with
the global E-field vector `E0` by using a [`GlobalJonesBasis`](@ref) and projecting the entries into the transverse plane defined
by the incoming ray direction and orthogonal E-field vector. This approach is partially inspired by the publication:

**Jan Korger et al., "The polarization properties of a tilted polarizer," Opt. Express 21, 27032-27042 (2013)**

!!! warning
    It is assumed that the ray direction of propagation is not changed during the interaction.

# Implementation reqs.

Subtypes of `AbstractJonesPolarizer` should implement all supertype requirements.

## Interaction logic

The [`GlobalJonesBasis`](@ref) tracks the rotation in 3D-space via the [`orientation`](@ref) of the attached [`AbstractShape`](@ref).
The polarization matrix `P` is calculated by projecting the previous matrix into the incoming orthogonal plane of polarization.
Refer to the [`JonesFilm`](@ref) for more information.

!!! info
    The validity of this approach is still under consideration for non-normal incidence.
"""
abstract type AbstractJonesPolarizer{T} <: AbstractObject{T} end

include("PolarizationFilter.jl")
include("LinearPolarizer.jl")