"""
    AbstractDetector <: AbstractObject

A generic representation of a detector that evaluates [`AbstractBeam`](@ref) data during and/or after interaction.
Refer to e.g. [`Detector`](@ref) for more information.

# Implementation reqs.

Subtypes of `AbstractDetector` should implement all supertype requirements as well as:

## Functions

- `interact3d`: see e.g. [`Detector`](@ref) for reference
- `empty!`: resets data stored in the detector, see below. [`initialize!`](@ref) calls it for all detectors of a system

# Additional information

The information provided below applies to the standard functional implementation of this type and may be overwritten
by specialized subtypes.

## Data mutability

In order to model field superposition effects, the concrete implementation of an `AbstractDetector` should be a **mutable struct**
with a `const`ant `shape` field. This is necessary, since (sub-)beams will interact sequentially with the detector during [`solve_system!`](@ref).
Only if the data can be accumulated sequentially, multiple beam interactions can be captured for a complex system, e.g. an interferometer.

## Data reset

Since e.g. E-field data is supposed to be accumulated by mutability of the detector data, the burden of resetting the data for a new solver call
is placed on the user. This function should be called `empty!`. A single detector is reset with `empty!(detector)`, all
detectors of a [`System`](@ref) or [`ObjectGroup`](@ref) with [`initialize!`](@ref), which calls `empty!` on every `AbstractDetector`.
The solver does not reset detectors unless asked to (keyword `initialize` of [`solve_system!`](@ref)), since the beams of a beam group
and several sources solved one after another are meant to superpose.

Every kinematic function empties a detector, since the stored data belong to the pose of the detector during the solve, see
[`translate3d!`](@ref)`(::Movable, ::AbstractDetector, offset)`. A subtype gets this from its `empty!` method without further code.
"""
abstract type AbstractDetector{T} <: AbstractObject{T} end

"""
    empty!(detector)

Resets the field data of the `detector`. Must be implemented for each concrete subtype of [`AbstractDetector`](@ref).
"""
function empty!(::D) where {D <: AbstractDetector}
    @warn "Detector reset logic for $D not implemented"
    return nothing
end

initialize!(detector::AbstractDetector) = (empty!(detector); nothing)

"""
    translate3d!(::Movable, detector::AbstractDetector, offset)

Empties the `detector` (see `empty!`) and moves it by `offset`. Its hits lie on the surface
where it was during the solve, hence they are discarded when it moves. All kinematic functions
are derived from this one and [`rotate3d!`](@ref)`(::Movable, ::AbstractDetector, R)`, so each of
them empties the detector, also when it moves as a member of an [`ObjectGroup`](@ref).
"""
function translate3d!(::Movable, detector::AbstractDetector, offset)
    empty!(detector)
    translate3d!(shape_trait_of(detector), detector, offset)
    return nothing
end

"""
    rotate3d!(::Movable, detector::AbstractDetector, R::AbstractMatrix)

Empties the `detector` (see `empty!`) and rotates it by `R` about its `position`, see
[`translate3d!`](@ref)`(::Movable, ::AbstractDetector, offset)`.
"""
function rotate3d!(::Movable, detector::AbstractDetector, R::AbstractMatrix)
    empty!(detector)
    rotate3d!(shape_trait_of(detector), detector, R)
    return nothing
end
