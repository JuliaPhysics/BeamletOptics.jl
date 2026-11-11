"""
    AbstractBeamsplitter <: AbstractObject

A generic type to represent an [`AbstractObject`](@ref) that splits incoming beams by reflection and transmission.

# Implementation reqs.

Subtypes of `AbstractBeamsplitter` should implement all supertype requirements.

## Interaction logic

A splitting surface is a thin layer whose [`surface_model`](@ref) splits, e.g. a
[`BeamsplitterCoating`](@ref) (see [`splitting_of`](@ref)): the solver appends the transmitted and the
reflected beam to the parent beam via [`children!`](@ref) and stops the parent.

!!! info "Appending convention"
    The `children` of a split beam are ordered
      1. transmitted beam
      2. reflected beam

## Functions

- `surface_model`: the splitting model of a thin layer, e.g. of a [`ThinBeamsplitter`](@ref); a
  composite holds such a layer, e.g. a [`CubeBeamsplitter`](@ref)
"""
abstract type AbstractBeamsplitter{T} <: AbstractObject{T} end

# order of inclusion matters
include("ThinBeamsplitter.jl")
include("PlateBeamsplitter.jl")
include("CubeBeamsplitter.jl")
include("Compensators.jl")