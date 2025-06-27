"""
    AbstractDetector <: AbstractObject

A generic representation of a detector that evaluates [`AbstractBeam`](@ref) data during interaction.
Refer to [`Photodetector`](@ref) for more information.

# Implementation reqs.

Subtypes of `AbstractDetector` should implement all supertype requirements as well as:

## Functions

- `interact3d`: see e.g. [`Photodetector`](@ref) for reference
- `empty!`: resets data stored in the detector, see below

# Additional information

The information provided below applies to the standard functional implementation of this type and may be overwritten
by specialized subtypes.

## Data mutability

In order to model field superposition effects, the concrete implementation of an `AbstractDetector` should be a **mutable struct**
with a `const`ant `shape` field. This is necessary, since (sub-)beams will interact sequentially with the detector during [`solve_system!`](@ref).
Only if the data can be accumulated sequentially, multiple beam interactions can be captured for a complex system, e.g. an interferometer.

## Data reset

Since e.g. E-field data is supposed to be accumulated by mutability of the detector data, the burden of resetting the data for a new solver call
is placed on the user. This function should be called `empty!`.
"""
abstract type AbstractDetector{T, S <: AbstractShape{T}} <: AbstractObject{T, S} end

"""
    empty!(detector)

Resets the field data of the `detector`. Must be implemented for each concrete subtype of [`AbstractDetector`](@ref).
"""
function empty!(::D) where D <: AbstractDetector
    @warn "Detector reset logic for $D not implemented"
    return nothing
end

# include concrete detectors
include("Photodetector.jl")
include("Spotdetector.jl")
include("PSFDetector.jl")
