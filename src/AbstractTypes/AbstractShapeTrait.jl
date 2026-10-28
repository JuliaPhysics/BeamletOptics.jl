"""
    AbstractShapeTrait

The shape trait defines how many shapes an [`AbstractObject`](@ref) consists of. 
Two different traits are defined:

1. [`SingleShape`](@ref): the `AbstractObject` consists of a single [`AbstractShape`](@ref)
2. [`MultiShape`](@ref): the `AbstractObject` consists of two or more [`AbstractShape`](@ref)s

Refer to the respective documentation for more information
"""
abstract type AbstractShapeTrait end

function shape_getter_not_implemented_error(T::Type, O::Type)
    throw(ErrorException("Shape getter of $O not implemented! Refer to $T documentation"))
end

"""
    SingleShape <: AbstractShapeTrait

Represents that the [`AbstractObject`](@ref) consists of a single underlying shape.

# AbstractObject implementation reqs.

If `shape_trait_of(::Foo) = SingleShape()` is defined, `Foo` must implement the following: 

## Fields

- `shape`: a single concrete [`AbstractShape`](@ref), e.g. a [`CylinderSDF`](@ref)
"""
struct SingleShape <: AbstractShapeTrait end

shape(::SingleShape, object::AbstractObject) = object.shape

Base.position(::SingleShape, object::AbstractObject) = position(shape(object))
position!(::SingleShape, object::AbstractObject, pos) = position!(shape(object), pos)

orientation(::SingleShape, object::AbstractObject) = orientation(shape(object))
orientation!(::SingleShape, object::AbstractObject, dir) = orientation!(shape(object), dir)

translate3d!(::SingleShape, object::AbstractObject, offset) = translate3d!(shape(object), offset)

rotate3d!(::SingleShape, object::AbstractObject, R::AbstractMatrix) = rotate3d!(shape(object), R)

"""
    MultiShape <: AbstractShapeTrait

Represents that the [`AbstractObject`](@ref) consists of a two or more [`AbstractShape`](@ref)s.

# AbstractObject implementation reqs.

If `shape_trait_of(::Foo) = MultiShape()` is defined, `Foo` must implement the following: 

## Functions

- `shape(::Foo)`: a getter function that returns a `Tuple` of all relevant shapes, e.g. `(foo.front, foo.middle, foo.back)`

All shapes returned by `shape(::Foo)` must be movable, see [`BeamletOptics.AbstractKinematicTrait`](@ref).

The parts may be shared: a part listed several times, or also held by another part (e.g. an
object on the same shape), is moved once per rigid motion, compared by identity (`===`), see
the carried motion in [`BeamletOptics.AbstractKinematicTrait`](@ref). The kinematic primitives of
`MultiShape` objects are their carried methods `translate3d!(object, offset, carried)` and
`rotate3d!(object, R, pivot, carried)`; a `MultiShape` type with its own kinematics overrides
these carried methods, since its parent moves it with them.

# Additional information

!!! warn "Kinematic center"
    Unless specified otherwise by dispatching `position` / `position!` and `orientation` / `orientation!`
    onto custom `pos` and `dir` data fields, the **position and orientation of the first element** returned
    by `shape(object)` will be used as the **kinematic center** for e.g. `translate3d!`.
""" 
struct MultiShape <: AbstractShapeTrait end

shape(::MultiShape, ::O) where O <: AbstractObject = shape_getter_not_implemented_error(MultiShape, O)

Base.position(::MultiShape, object::AbstractObject) = position(first(shape(object)))
position!(::MultiShape, object::AbstractObject, ::Any) = nothing

orientation(::MultiShape, object::AbstractObject) = orientation(first(shape(object)))
orientation!(::MultiShape, object::AbstractObject, ::Any) = nothing

"""
    translate3d!(::MultiShape, object, offset)

Moves all parts of the [`MultiShape`](@ref) `object` along the specified `offset` vector, each
part once. The parts are carried along, see [`MultiShape`](@ref).
"""
function translate3d!(::MultiShape, object::AbstractObject, offset)
    translate3d!(object, offset, IdSet{Any}())
    return nothing
end

"""
    rotate3d!(::MultiShape, object, R::AbstractMatrix)

All parts of the [`MultiShape`](@ref) `object` are rotated around the pivot center via the rotation matrix `R`,
each part once. The parts are carried along, see [`MultiShape`](@ref).
"""
function rotate3d!(::MultiShape, object::AbstractObject, R::AbstractMatrix)
    rotate3d!(object, R, position(object), IdSet{Any}())
    return nothing
end

# Carried motion of objects: dispatch on the shape trait

function translate3d!(object::AbstractObject, offset, carried::IdSet)
    translate3d!(shape_trait_of(object), object, offset, carried)
    return nothing
end

function rotate3d!(object::AbstractObject, R::AbstractMatrix, pivot::AbstractVector, carried::IdSet)
    rotate3d!(shape_trait_of(object), object, R, pivot, carried)
    return nothing
end

# The shape is carried with the same set, so that an operand of a composite SDF that is also
# listed elsewhere moves once
function translate3d!(::SingleShape, object::AbstractObject, offset, carried::IdSet)
    _mark_carried!(carried, object) || return nothing
    translate3d!(shape(object), offset, carried)
    return nothing
end

function rotate3d!(::SingleShape, object::AbstractObject, R::AbstractMatrix, pivot::AbstractVector,
        carried::IdSet)
    _mark_carried!(carried, object) || return nothing
    rotate3d!(shape(object), R, pivot, carried)
    return nothing
end

function translate3d!(::MultiShape, object::AbstractObject, offset, carried::IdSet)
    _mark_carried!(carried, object) || return nothing
    # Translate tracking vector
    position!(object, position(object) .+ offset)
    # Recursively translate all subparts
    for subpart in shape(object)
        translate3d!(subpart, offset, carried)
    end
    return nothing
end

function rotate3d!(::MultiShape, object::AbstractObject, R::AbstractMatrix, pivot::AbstractVector,
        carried::IdSet)
    _mark_carried!(carried, object) || return nothing
    # Read the tracking position first, by default it is the position of the first subpart
    p = position(object)
    orientation!(object, R * orientation(object))
    position!(object, pivot + R * (p - pivot))
    # Recursively rotate all subgroups and objects around the pivot point
    for subpart in shape(object)
        rotate3d!(subpart, R, pivot, carried)
    end
    return nothing
end
