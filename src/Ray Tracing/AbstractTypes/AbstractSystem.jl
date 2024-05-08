"""
    AbstractSystem

A generic representation of a system of optical elements.

# Implementation reqs.

Subtypes of `AbstractBeam` must implement the following:

## Fields:

- `objects`: a vector or tuple of [`AbstractObject`](@ref)s that make up the system
- `n`: (optional) refractive index of the surrounding medium, default value is 1.0

## Functions:

- `refractive_index`: returns the refractive index `n`, see above
"""
abstract type AbstractSystem end

refractive_index(::AbstractSystem) = 1.0

"""
    interact3d(::AbstractSystem, object::AbstractObject, ::AbstractBeam)

Defines the optical interaction between an incoming/outgoing beam/ray of light and an optical element, must return an [`Interaction`](@ref) or `nothing`.
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
in the [`AbstractSystem`](@ref) will be hit next. The hint does not need to result in a guaranteed [`Intersection`](@ref).

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

Describes how an `AbstractBeam` and an `AbstractSystem` interact with each other.
Contains the parameters for the continuation of the `AbstractBeam`.
Can store an optional `id` which hints at the next object intersection.
"""
abstract type AbstractInteraction end

hint(i::AbstractInteraction)::Nullable{Hint} = i.hint
hint(::Nothing) = nothing