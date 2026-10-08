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

!!! warning
    Whenever your object consists of nested structures (e.g. other [`MultiShape`](@ref)s) 
    it is the responsibility of the user to ensure that **each atomic shape, i.e. `SingleShape`, is only listed once**.
    Failure to ensure this can lead to spurious behaviour when using the kinematic API.
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

"""
    AbstractObjectGroup

Container type for groups of optical elements, based on a tree-like data structure. Intended for easier kinematic handling of connected elements.
See also [`ObjectGroup`](@ref) for a concrete implementation.

A group is no [`AbstractObject`](@ref): it has no shape, no [`AbstractShapeTrait`](@ref) and no
[`interact3d`](@ref) method. A ray or beam interacts with the object of the group that it hits, never
with the group. A component that consists of several parts which interact as one, e.g. a cemented
doublet, is therefore an `AbstractObject` with the [`MultiShape`](@ref) trait and not a group.

# Implementation reqs.

Subtypes of `AbstractObjectGroup` must implement the following:

## Fields:

- `objects`: a tuple of [`AbstractObject`](@ref)s or additional subgroups, allows for hierarchical structures
- `center`: a `Point3{T}` which is regarded as the reference origin (pivot) of the group
- `dir`: a `SMatrix{3,3,T,9}` that describes the local coordinate system of the group

Since the kinematic API modifies `center` and `dir`, the subtype must be a `mutable struct`.

## Functions:

If the fields above do not exist, the following getters/setters must be dispatched:

- [`objects`](@ref): getter for the `objects` field or equivalent return type
- `position` / `position!`: gets or sets the reference origin (pivot)
- [`orientation`](@ref) / `orientation!`: gets or sets the orientation matrix

## Kinematic

An `AbstractObjectGroup` is a container: it takes the kinematic class of its `objects`, i.e.
[`BeamletOptics.Movable`](@ref) with an [`BeamletOptics.Oriented`](@ref) frame for movable objects, see
[`BeamletOptics.AbstractKinematicTrait`](@ref). The constructors must check that the `objects` are either
all static or all movable. The primitives `translate3d!(::Movable, group, offset)` and
`rotate3d!(::Movable, group, R::AbstractMatrix)` move all `objects` and the pose of the group.
"""
abstract type AbstractObjectGroup{T} end

"""
    ObjectOrGroup

Union of [`AbstractObject`](@ref) and [`AbstractObjectGroup`](@ref): what a [`System`](@ref) or a
group stores and what [`render!`](@ref) draws. The two types share no supertype.
"""
const ObjectOrGroup = Union{AbstractObject, AbstractObjectGroup}

"""
    objects(group::AbstractObjectGroup)

Exposes all objects/subgroups stored within the group.
"""
objects(group::AbstractObjectGroup) = group.objects

AbstractTrees.children(group::AbstractObjectGroup) = objects(group)

kinematic_trait_of(group::AbstractObjectGroup) = _container_trait(objects(group))

Base.position(group::AbstractObjectGroup) = group.center
position!(group::AbstractObjectGroup{T}, pos) where {T} = (group.center = Point3{T}(pos))

orientation(group::AbstractObjectGroup) = group.dir
orientation!(group::AbstractObjectGroup{T}, dir) where {T} = (group.dir = SMatrix{3, 3, T, 9}(dir))

function initialize!(group::AbstractObjectGroup)
    foreach(initialize!, objects(group))
    return nothing
end

"""
    translate3d!(::Movable, group::AbstractObjectGroup, offset)

Moves all objects of the `group` and its `center` by `offset`.
"""
function translate3d!(::Movable, group::AbstractObjectGroup, offset)
    position!(group, position(group) .+ offset)
    foreach(member -> translate3d!(member, offset), objects(group))
    return nothing
end

"""
    rotate3d!(::Movable, group::AbstractObjectGroup, R::AbstractMatrix)

Rotates all objects of the `group` by `R` about the group `center` and updates the group
[`orientation`](@ref) to `R * orientation(group)`.
"""
function rotate3d!(::Movable, group::AbstractObjectGroup, R::AbstractMatrix)
    orientation!(group, R * orientation(group))
    pivot = position(group)
    foreach(member -> rotate3d!(member, R, pivot), objects(group))
    return nothing
end
