"""
    AbstractSystem

A generic representation of a system of optical elements.

# Implementation reqs.

Subtypes of `AbstractSystem` must implement the following:

## Fields:

- `objects`: a vector or tuple of [`AbstractObject`](@ref)s that make up the system
- `n`: (optional) [`RefractiveIndex`](@ref) of the surrounding medium, default value is 1.0

## Functions:

- `refractive_index`: returns the [`RefractiveIndex`](@ref) `n` of the system medium, see above
"""
abstract type AbstractSystem end

refractive_index(::AbstractSystem, λ::Real) = 1.0

"""
    interact3d(::AbstractSystem, object::AbstractObject, ::AbstractBeam)

Defines the optical interaction between an incoming/outgoing beam/ray of light and an optical element, must return an [`AbstractInteraction`](@ref) or `nothing`.
The default behavior is that no interaction occurs, i.e. return of `nothing`, which should stop the system tracing procedure.
Refer to the [`AbstractInteraction`](@ref) typedocs for more information on the return type value.
"""
function interact3d(::AbstractSystem, ::ObjectType, ::AbstractBeam, ::RayType) where {ObjectType<:AbstractObject, RayType<:AbstractRay}
    @warn lazy"No interact3d method defined for:" ObjectType RayType
    return nothing
end

"""
    Hint

A `Hint` can be passed as part of an [`AbstractInteraction`](@ref) and will inform the tracing algorithm about which [`AbstractObject`](@ref)
in the [`AbstractSystem`](@ref) will be hit next.

!!! info 
    The hint does not need to result in a guaranteed [`Intersection`](@ref).

# Fields

- `object`: the object that might or will be intersected as next
- `shape`: the underlying shape that will be intersected next, i.e. `shape(object)`, relevant for multi-shape objects
"""
struct Hint
    object::AbstractObject
    shape::AbstractShape
end

Hint(o::AbstractObject) = Hint(o, shape(o))

object(h::Hint) = h.object
shape(h::Hint) = h.shape

"""
    AbstractInteraction

Describes how an [`AbstractBeam`](@ref) and an [`AbstractObject`](@ref) interact with each other.
This data type stores information from the [`interact3d`](@ref) function and provides it to the solver.
The solver can use this data to extend the [`AbstractBeam`](@ref).

# Implementation reqs.

Subtypes of `AbstractInteraction` must implement the following:

## Fields

- `hint`: a nullable [`Hint`](@ref) for the solver (optional but recommended)

## Beam data

It is required that concrete implementations of this type provide some form of data on how to extend the beam.
For instance, refer to [`BeamInteraction`](@ref) and [`GaussianBeamletInteraction`](@ref).
"""
abstract type AbstractInteraction end

hint(i::AbstractInteraction)::Nullable{Hint} = i.hint
hint(::Nothing) = nothing
