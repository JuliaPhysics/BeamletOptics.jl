"""
    AbstractBeamsplitter <: AbstractObject

A generic type to represent an [`AbstractObject`](@ref) that splits incoming beams by reflection and transmission.

# Implementation reqs.

Subtypes of `AbstractBeamsplitter` should implement all supertype requirements.

## Interaction logic

After intersection with the [`AbstractShape`](@ref) at which the beam splitting occurs, the `interact3d` function 
should appended the transmitted and reflected sub-beams to the parent beam via the [`children!`](@ref) function.
The `interact3d` function should then return `nothing` in order to stop the tracing of the parent beam.

!!! info "Appending convention"
    As a convention, when splitting an incoming beam, the order of `children` appended to the parent beam should be
      1. transmitted beam
      2. reflected beam

## Functions

- `interact3d`: see above
- `_beamsplitter_transmitted_beam`: optional helper function
- `_beamsplitter_reflected_beam`: optical helper function
"""
abstract type AbstractBeamsplitter{T, S <: AbstractShape{T}} <: AbstractObject{T, S} end

# order of inclusion matters
include("ThinBeamsplitter.jl")
include("PlateBeamsplitter.jl")
include("CubeBeamsplitter.jl")
include("Compensators.jl")