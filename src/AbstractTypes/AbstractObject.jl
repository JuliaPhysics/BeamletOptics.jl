"""
    AbstractObject

A generic type for 2D/3D objects that can be used to model optical elements. The geometry of the object is represented via an [`AbstractShape`](@ref).
The optical effect that occurs between the object and an incoming ray/beam of light is modeled via its [`interact3d`](@ref) method.

# Implementation reqs.

Subtypes of `AbstractObject` must implement the following:

## Shape trait

An `AbstractObject` can consist of a single [`AbstractShape`](@ref), e.g. a lens element, or a collection of functionally dependant shapes, e.g. a cube beamsplitter.
In order to model this, the API implementation of an `AbstractObject` requires the definition of the **shape trait**. This trait allows the dispatch onto specialized methods to handle
the kinematic interface and tracing methods for objects consisting of one or more shapes.

- `shape_trait_of`: defines the shape type of the `AbstractObject`, refer to [`AbstractShapeTrait`](@ref) for more information

!!! info "Default shape trait"
    Unless specified otherwise, the `shape_trait_of` an `AbstractObject` is defined as [`SingleShape`](@ref). This requires `object.shape` as a dedicated field.
    For [`MultiShape`](@ref)s the getter function `shape(object)` must return a tuple of all shapes that make up the object.

## Getters/setters

All kinematic functions defined for the [`AbstractShape`](@ref) can also be called for a `AbstractObject`. In this case, the shape trait will define how the specific movement function
is dispatched.

## Kinematic

`AbstractObject`s are [`BeamletOptics.Movable`](@ref) with an [`BeamletOptics.Oriented`](@ref) frame, see
[`BeamletOptics.AbstractKinematicTrait`](@ref). The primitives `translate3d!(::Movable, object, offset)` and
`rotate3d!(::Movable, object, R::AbstractMatrix)` forward to the shape trait. A subtype can opt out of the
kinematic API via `BeamletOptics.kinematic_trait_of(::Foo) = BeamletOptics.Static()`.

A container (e.g. an [`ObjectGroup`](@ref)) moves its objects with the carried methods
`translate3d!(object, offset, carried)` and `rotate3d!(object, R, pivot, carried)`, not with the primitives, see
`translate3d!(::AbstractShape, offset, carried)`. A subtype that overrides the primitives to do more than move its
shape (e.g. an [`AbstractDetector`](@ref) empties itself) must override these carried methods as well, otherwise
that extra step is skipped whenever the object moves as part of a container:

```julia
function BeamletOptics.translate3d!(x::Foo, offset, carried::IdSet)
    extra_step!(x)
    BeamletOptics.translate3d!(BeamletOptics.shape_trait_of(x), x, offset, carried)
    return nothing
end
```

## Functions:

- [`interact3d`](@ref): defines the optical interaction, the return type must be `Nothing` or an [`AbstractInteraction`](@ref)
- [`initialize!`](@ref) (optional): resets the data that the object stores during a solve, e.g. the hits of a detector. By default nothing is reset.
"""
abstract type AbstractObject{T <: Real} end

"""
    initialize!(x)

Brings `x` into the state that a new solve starts from, by discarding the data that its objects
have stored during previous calls of [`solve_system!`](@ref). `x` can be a [`System`](@ref) or
[`StaticSystem`](@ref), an [`ObjectGroup`](@ref) (nested groups are included) or a single object.
Currently this empties all [`Detector`](@ref)s. Objects that store no such data are left
unchanged. Returns `nothing`.

`solve_system!` does not call this function unless its keyword `initialize` is set: a
[`Detector`](@ref) accumulates the hits of several solves on purpose, such that several sources
solved one after another superpose. Initialize before solving again after the setup was changed,
e.g. in a parameter scan. The beams are not touched, `solve_system!` resets them itself.

A new object type that stores data during a solve implements `initialize!(object)`. For an
[`AbstractDetector`](@ref) it calls `empty!(detector)`.

# Examples

```julia
for shift in shifts
    translate_to3d!(mirror, shift)
    initialize!(system)             # instead of empty!(pd) for every detector
    solve_system!(system, beam)
    P = optical_power(pd)
end
```

or in one call via `solve_system!(system, beam; initialize = true)`.
"""
initialize!(::AbstractObject) = nothing

"Default trait"
shape_trait_of(::AbstractObject) = SingleShape()

kinematic_trait_of(::AbstractObject) = Movable(Oriented())

"""
    shape(::AbstractObject)

Returns all component shapes of the object for a [`MultiShape`](@ref) or a single shape for a [`SingleShape`](@ref).
E.g. for a custom multi-shape object the user of this API needs to define:

`BeamletOptics.shape(obj::MyObject) = (obj.front, obj.back)`

!!! note "Shared parts"
    A part may be listed several times, e.g. through nested [`MultiShape`](@ref)s, or be shared with
    another object, e.g. two fields that hold the same shape or an operand of a composite SDF that is
    also listed on its own. The kinematic API moves every part once per motion, compared by identity
    (`===`), see the carried motion in [`BeamletOptics.AbstractKinematicTrait`](@ref).
"""
shape(object::AbstractObject) = shape(shape_trait_of(object), object)

"""
    position(object) -> Point3

Returns the current `position` of the `object` in R³ as a `Point3` where (x, y, z) in a right-hand coordinate system.

In general, `position(object)` returns `position(shape(object))` unless specified otherwise.
"""
Base.position(object::AbstractObject) = position(shape_trait_of(object), object)
position!(object::AbstractObject, pos) = position!(shape_trait_of(object), object, pos)

"""
    orientation(object) -> Matrix

Returns the current `orientation` of the `object` in R³ as a matrix.
The matrix represents the local fixed-body coordinate system.

In general, `orientation(object)` returns `orientation(shape(object))` unless specified otherwise.
"""
orientation(object::AbstractObject) = orientation(shape_trait_of(object), object)
orientation!(object::AbstractObject, dir) = orientation!(shape_trait_of(object), object, dir)

translate3d!(::Movable, object::AbstractObject, offset) = translate3d!(shape_trait_of(object), object, offset)

rotate3d!(::Movable, object::AbstractObject, R::AbstractMatrix) = rotate3d!(shape_trait_of(object), object, R)
