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

# A system is not moved as a whole, its objects are
kinematic_trait_of(::AbstractSystem) = Static()

"""
    interact3d(::AbstractSystem, object::AbstractObject, ::AbstractBeam, ::AbstractRay)

Defines the optical interaction between an incoming/outgoing beam/ray of light and an optical element, must return an [`AbstractInteraction`](@ref) or `nothing`.
The default behavior is that no interaction occurs, i.e. return of `nothing`, which should stop the system tracing procedure.
Refer to the [`AbstractInteraction`](@ref) typedocs for more information on the return type value.

A [`MultiShape`](@ref) object receives this call for whichever of its parts was hit. Its method
finds that part, e.g. via `shape(intersection(ray))`, and hands the beam to it, see
[`CubeBeamsplitter`](@ref). Parts that do not interact as one belong in an
[`ObjectGroup`](@ref) or another [`AbstractObjectGroup`](@ref) instead, where each part
interacts on its own.
"""
function interact3d(::AbstractSystem, object::AbstractObject, ::AbstractBeam, ray::AbstractRay)
    warn_no_interact3d(shape_trait_of(object), object, ray)
    return nothing
end

# Warning of the interact3d fallback
warn_no_interact3d(::AbstractShapeTrait, object, ray) =
    @warn lazy"No interact3d method defined for:" ObjectType = typeof(object) RayType = typeof(ray)
warn_no_interact3d(::MultiShape, object, ray) = @warn lazy"""
    No interact3d method defined for the MultiShape object $(nameof(typeof(object))).
    If its parts interact as one (e.g. cemented, with a Hint across a shared face), define
    interact3d and hand the beam to the part that was hit, like CubeBeamsplitter.
    If its parts interact on their own, make it an ObjectGroup or a subtype of
    AbstractObjectGroup, which needs no interact3d.""" RayType = typeof(ray)

"""
    Hint

A `Hint` can be passed as part of an [`AbstractInteraction`](@ref) and will inform the tracing algorithm about which [`AbstractObject`](@ref)
in the [`AbstractSystem`](@ref) will be hit next.

!!! info 
    The `Hint` does not need to result in a guaranteed [`Intersection`](@ref). However, if the hinted shape is intersected, it will
    be immediately assumed as the correct global intersection.

# Fields

- `object`: the object that might or will be intersected next
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
For instance, refer to the[`BeamInteraction`](@ref) and [`GaussianBeamletInteraction`](@ref).
"""
abstract type AbstractInteraction end

hint(i::AbstractInteraction)::Nullable{Hint} = i.hint
hint(::Nothing) = nothing

hint!(i::AbstractInteraction, new_hint::Nullable{Hint}) = (i.hint = new_hint)
