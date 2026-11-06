"""
    AbstractSystem

A generic representation of a system of optical elements.

# Implementation reqs.

Subtypes of `AbstractSystem` must implement the following:

## Fields:

- `objects`: a vector or tuple of [`AbstractObject`](@ref)s and object groups ([`AbstractObjectGroup`](@ref)) that make up the system.
  The solver does not read this field itself but calls the function `objects`, see below

## Functions:

- `objects(system)`: returns the objects of the system with object groups flattened into their objects. There is no fallback
- `refractive_index(system, λ::Real)` (optional): returns the refractive index of the surrounding medium at the wavelength `λ`.
  The default is `1.0`; a system does not store a refractive index

The solver tests a ray against the top-level objects and groups of the system. For a custom system these are the result of `objects(system)`.
The built-in [`System`](@ref) tests its groups as a whole, a [`StaticSystem`](@ref) flattens them when it is built.
"""
abstract type AbstractSystem end

refractive_index(::AbstractSystem, λ::Real) = 1.0

# A system is not moved as a whole, its objects are
kinematic_trait_of(::AbstractSystem) = Static()

"""
    interact3d(::AbstractSystem, object::AbstractObject, beam::AbstractBeam, at)

Defines the optical interaction between an incoming/outgoing beam/ray of light and an optical element, must return an [`AbstractInteraction`](@ref) or `nothing`.
Without a method of its own, a body (see [`medium`](@ref)) or a thin layer with a model (see
[`surface_model`](@ref)) acts at the boundary at its hit, see `interact3d(system, t::Transition, beam, at)`:
a new body or thin layer declares its medium or model and needs no `interact3d` method. Any other
object warns and returns `nothing`, which stops the beam.
Refer to the [`AbstractInteraction`](@ref) typedocs for more information on the return type value.

The media on both sides of the boundary at the hit are given by [`transition`](@ref)`(system, ray)`,
also where several elements touch, e.g. at the cemented faces of a doublet.

The beam interacts `at` its ray that hit the object, i.e. a [`Beam`](@ref) at the `ray` itself and a
beamlet at the index `ray_id` of the segment of its component rays, see `intersection(beam, at)`.

A [`MultiShape`](@ref) object whose parts are objects needs no method: the call is handed to the part
that was hit, see [`part_of`](@ref), e.g. to the front or back lens of a [`DoubletLens`](@ref). A
`MultiShape` object of plain shapes defines its own method.
"""
interact3d(system::AbstractSystem, object::AbstractObject, beam::AbstractBeam,
        at::Union{AbstractRay, Int}) =
    interact3d(system, shape_trait_of(object), object, beam, at)

# A body or a thin layer acts at the boundary, see `surface_model`; any other object blocks the beam
function interact3d(system::AbstractSystem, ::SingleShape, object::AbstractObject, beam::AbstractBeam,
        ray::AbstractRay)
    if isnothing(medium(object)) && isnothing(surface_model(object))
        warn_no_interact3d(shape_trait_of(object), object, ray)
        return nothing
    end
    return interact3d(system, transition(system, ray), beam, ray)
end

function interact3d(system::AbstractSystem, ::MultiShape, object::AbstractObject,
        beam::AbstractBeam, at::Union{AbstractRay, Int})
    part = part_of(intersection(beam, at)::Intersection)
    part === object && return interact3d(system, SingleShape(), object, beam, at)
    return interact3d(system, part, beam, at)
end

# Warning of the interact3d fallback
warn_no_interact3d(::AbstractShapeTrait, object, ray) = @warn lazy"""
    No interact3d method defined for $(nameof(typeof(object))): define interact3d, or declare its
    medium (a body) or surface_model (a thin layer).""" RayType = typeof(ray)
warn_no_interact3d(::MultiShape, object, ray) = @warn lazy"""
    No interact3d method defined for the MultiShape object $(nameof(typeof(object))).
    Its parts are plain shapes: define interact3d for it, or build it from objects with
    their own interact3d (e.g. lenses), to which the solver hands the beam.""" RayType = typeof(ray)

"""
    Hint

A `Hint` can be passed as part of an [`AbstractInteraction`](@ref) and will inform the tracing algorithm about which [`AbstractObject`](@ref)
in the [`AbstractSystem`](@ref) will probably be hit next.

!!! info
    The `Hint` only speeds up the tracing: the hinted object is tested first, and objects beyond its hit are skipped.
    It does not change which object is hit, and touching elements (e.g. cemented lenses) need no `Hint`, since the
    solver finds all surfaces at the hit point, see [`coincident_hits`](@ref) and [`transition`](@ref).

# Fields

- `object`: the object that might or will be intersected next. The solver calls `interact3d` for this object, hence for a part of a
  [`MultiShape`](@ref) object it is the composite itself, not the part
- `shape`: the underlying shape that will be intersected next, i.e. `shape(object)`, relevant for multi-shape objects

`Hint(object)` uses `shape(object)`, which is a tuple for a [`MultiShape`](@ref) object and fails. Pass the hit part explicitly there,
`Hint(object, part)`, with the composite as `object` and the shape that will be hit as `part`.
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

- `hint`: a nullable [`Hint`](@ref) for the solver. The field is required, since the solver reads it for every interaction, but it may be `nothing`

## Beam data

It is required that concrete implementations of this type provide some form of data on how to extend the beam.
For instance, refer to the[`BeamInteraction`](@ref) and [`GaussianBeamletInteraction`](@ref).
"""
abstract type AbstractInteraction end

hint(i::AbstractInteraction)::Nullable{Hint} = i.hint
hint(::Nothing) = nothing

hint!(i::AbstractInteraction, new_hint::Nullable{Hint}) = (i.hint = new_hint)
