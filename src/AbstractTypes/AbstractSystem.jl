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
    interact3d(::AbstractSystem, object::AbstractObject, beam::AbstractBeam, at)

Defines the optical interaction between an incoming/outgoing beam/ray of light and an optical element, must return an [`AbstractInteraction`](@ref) or `nothing`.
The default behavior is that no interaction occurs, i.e. return of `nothing`, which should stop the system tracing procedure.
Refer to the [`AbstractInteraction`](@ref) typedocs for more information on the return type value.

The media on both sides of the boundary at the hit are given by [`transition`](@ref)`(system, ray)`,
also where several elements touch, e.g. at the cemented faces of a doublet.

The beam interacts `at` its ray that hit the object, i.e. a [`Beam`](@ref) at the `ray` itself and a
beamlet at the index `ray_id` of the segment of its component rays, see `intersection(beam, at)`.

A [`MultiShape`](@ref) object whose parts are objects needs no method: the call is handed to the part
that was hit, see [`part_of`](@ref), e.g. to the front or back lens of a [`DoubletLens`](@ref). A
`MultiShape` object of plain shapes defines its own method.
"""
interact3d(system::AbstractSystem, object::AbstractObject, beam::AbstractBeam, at::Union{AbstractRay, Int}) =
    interact3d(system, shape_trait_of(object), object, beam, at)

function interact3d(::AbstractSystem, ::SingleShape, object::AbstractObject, ::AbstractBeam,
        ray::AbstractRay)
    _warn_no_interact3d(object, ray)
    return nothing
end

function interact3d(system::AbstractSystem, ::MultiShape, object::AbstractObject,
        beam::AbstractBeam, at::Union{AbstractRay, Int})
    part = part_of(intersection(beam, at)::Intersection)
    part === object && return interact3d(system, SingleShape(), object, beam, at)
    return interact3d(system, part, beam, at)
end

# Warning of the fallback above
_warn_no_interact3d(object, ray) =
    @warn lazy"No interact3d method defined for:" ObjectType = typeof(object) RayType = typeof(ray)

"""
    Hint

A `Hint` can be passed as part of an [`AbstractInteraction`](@ref) and will inform the tracing algorithm about which [`AbstractObject`](@ref)
in the [`AbstractSystem`](@ref) will probably be hit next.

!!! info
    The `Hint` only speeds up the tracing: the hinted object is tested first, and objects beyond its hit are skipped.
    It does not change which object is hit, and touching elements (e.g. cemented lenses) need no `Hint`, since the
    solver finds all surfaces at the hit point, see [`coincident_hits`](@ref) and [`transition`](@ref).

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
