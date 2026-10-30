"""
    AbstractKinematicTrait

The kinematic trait defines whether an entity can be moved by the kinematic API
([`translate3d!`](@ref), [`rotate3d!`](@ref), ...). It is returned by
[`kinematic_trait_of`](@ref) and follows the same dispatch
pattern as [`AbstractShapeTrait`](@ref): every public function e.g. `translate3d!(x, args...)` forwards to
`translate3d!(kinematic_trait_of(x), x, args...)`.

Two traits are defined:

1. [`Static`](@ref): `x` cannot be moved, every kin. function throws an `ArgumentError`
2. [`Movable`](@ref): `x` can be moved; its frame
   ([`Oriented`](@ref) or [`Directed`](@ref))
   selects orientation- or direction-based behaviour

The trait is opt-in: shapes, objects, rays, beams, beam groups and object groups declare it on
their abstract type, any other type with its own `kinematic_trait_of` method. For a type without
a method, `kinematic_trait_of` and every kinematic function throw an `ErrorException`.

A subtype of a movable abstract type can opt out of the kinematic API by declaring its own
`kinematic_trait_of` method, e.g.

```julia
struct FixedMirror{T, S <: AbstractShape{T}} <: AbstractObject{T}
    shape::S
end

kinematic_trait_of(::FixedMirror) = Static()
```

Every kinematic function on a `FixedMirror` now throws an `ArgumentError`; its getters
(`position`, `orientation`, ...) stay available.

# Containers

All members of a container ([`ObjectGroup`](@ref), [`UnionSDF`](@ref)/[`DifferenceSDF`](@ref),
[`AbstractBeamGroup`](@ref)) must be either all `Static` or all `Movable`. Frames can still be mixed.
**This must be checked in the respective constructors**; nesting is then covered because a nested container already has a uniform
class from its own constructor. The container then takes that class: `Static` if all members
are `Static`, else `Movable(`[`Oriented`](@ref)`())`. All shapes returned
by `shape(object)` of a [`MultiShape`](@ref) object must be movable, since `MultiShape` objects
(e.g. a `CubeBeamsplitter`) have no shared constructor to enforce this.

A container moves each of its members once per motion, compared by identity (`===`). Two members
that reach the same shape, or a static value at any level, throw an `ArgumentError` before
anything moves, see `shape`.
"""
abstract type AbstractKinematicTrait end

"""
    Static <: AbstractKinematicTrait

Represents a value that cannot be moved. Every kinematic function throws an `ArgumentError`; the
getters (e.g. `position`, [`orientation`](@ref)) stay available.
"""
struct Static <: AbstractKinematicTrait end

"""
    AbstractKinematicFrame

Reference frame of a [`Movable`](@ref) value, either [`Oriented`](@ref) or [`Directed`](@ref).
"""
abstract type AbstractKinematicFrame end

"""
    Oriented <: AbstractKinematicFrame

Reference frame of a movable value with a full local coordinate system, e.g. shapes, objects and beam
groups.

# Implementation reqs.

- `position` / `position!`
- [`orientation`](@ref) / `orientation!`: a right-handed orthonormal 3x3 matrix

[`direction`](@ref) returns the local y-axis, i.e. `orientation(x)[:, 2]`.
[`reset_rotation3d!`](@ref) rotates `x` by `transpose(orientation(x))` about its `position`
and then sets the orientation to the identity.
"""
struct Oriented <: AbstractKinematicFrame end

"""
    Directed <: AbstractKinematicFrame

Reference frame of a movable value that only has a direction but no orientation, e.g. rays and beams.

# Implementation reqs.

- `position`
- [`direction`](@ref): the direction getter of the type

[`reset_rotation3d!`](@ref) throws an `ArgumentError`, use [`align3d!`](@ref) instead.
"""
struct Directed <: AbstractKinematicFrame end

"""
    Movable{F <: AbstractKinematicFrame} <: AbstractKinematicTrait

Represents a value that can be moved by the kinematic API.
The `frame` field ([`Oriented`](@ref) or [`Directed`](@ref))
selects orientation- or direction-based behaviour of the derived functions.

Rotations are applied about the own `position` of the value unless a `pivot` is passed
explicitly. A rotation `axis` can have any non-zero length (it is normalized internally).

# Kinematic functions

- [`translate3d!`](@ref): moves `x` by an `offset` vector (primitive, implemented per type)
- [`rotate3d!`](@ref): rotates `x` by a rotation matrix `R` about `position(x)`, or about a
  given `pivot` if one is passed (the `R`-only form is the primitive, implemented per type;
  the axis/angle and pivot forms are derived from it)
- [`translate_to3d!`](@ref): moves `x` so that `position(x)` coincides with a target point
- [`xrotate3d!`](@ref) / [`yrotate3d!`](@ref) / [`zrotate3d!`](@ref): rotate `x` by an angle
  `θ` [rad] about the global x-, y- or z-axis through `position(x)`
- [`align3d!`](@ref): rotates `x` about `position(x)` so that [`direction`](@ref)`(x)` is
  aligned with a target vector
- [`reset_translation3d!`](@ref): moves `x` so that `position(x)` is the global origin;
  available for every `Movable`, including rays, beams and beam groups
- [`reset_rotation3d!`](@ref): only available for [`Oriented`](@ref)
  values, where it rotates `x` back to identity `orientation`; for
  [`Directed`](@ref) values (rays, beams) it throws an
  `ArgumentError`, use [`align3d!`](@ref) instead
- [`direction`](@ref): for [`Oriented`](@ref) values this is the local
  y-axis, `orientation(x)[:, 2]`; [`Directed`](@ref) values provide
  their own getter

Sources (rays, beams and beam groups) are additionally reset to their untraced start state by
every one of the functions above: a moved but already-traced source would otherwise keep rays or
child beams belonging to the old, now geometrically wrong, light path. Moving a beam that is
not the root of its beam tree (e.g. one created by a beamsplitter interaction) throws an
`ArgumentError`; move its root beam instead.

# Implementation reqs.

If `kinematic_trait_of(::Foo) = Movable(...)` is declared, `Foo` must implement the following:

- `position`
- `translate3d!(::Movable, x::Foo, offset)`: moves `x` by `offset`
- `rotate3d!(::Movable, x::Foo, R::AbstractMatrix)`: rotates `x` by `R` about `position(x)`

All other kin. functions listed above are derived from these two primitives.
"""
struct Movable{F <: AbstractKinematicFrame} <: AbstractKinematicTrait
    frame::F
end

"""
    kinematic_trait_of(x) -> AbstractKinematicTrait

Returns the [`AbstractKinematicTrait`](@ref) of `x`. A type takes part in the kinematic API only
if it opts in with a typed method: movable types declare [`Movable`](@ref)`(`[`Oriented`](@ref)`())`
or `Movable(`[`Directed`](@ref)`())`, static types [`Static`](@ref)`()`. Throws an
`ErrorException` for any other type.
"""
kinematic_trait_of(x) = error(lazy"kinematic_trait_of not defined for $(typeof(x))")

"""
    check_kinematic_members(members)

Throws an `ArgumentError` if some `members` of a container are [`Static`](@ref)
and some are [`Movable`](@ref). Empty collections pass. The constructor of an own
[`AbstractObjectGroup`](@ref) calls it on its objects.
"""
function check_kinematic_members(members::Union{Tuple, AbstractArray})
    allequal(m -> _container_trait(kinematic_trait_of(m)), members) && return nothing
    throw(ArgumentError("container members must be either all static or all movable, see `kinematic_trait_of`"))
end

"""
    _container_trait(members)
    _container_trait(trait::AbstractKinematicTrait)

Returns the kinematic trait of a container with homogeneous `members` (see
`check_kinematic_members`), or of a container whose members have the trait `trait`:
`Static()` for static members, else `Movable(Oriented())`. An empty container is movable.
"""
function _container_trait(members::Union{Tuple, AbstractArray})
    isempty(members) && return Movable(Oriented())
    return _container_trait(kinematic_trait_of(first(members)))
end

_container_trait(::Static) = Static()
_container_trait(::Movable) = Movable(Oriented())

"""
    _move_rigidly!(frames, parts, offset)
    _move_rigidly!(frames, parts, R::AbstractMatrix, pivot::AbstractVector)

Moves a container by `offset`, or rotates it by `R` about `pivot`: first the pose of every frame,
then every part with the public [`translate3d!`](@ref) or [`rotate3d!`](@ref). The caller passes
each part once.
"""
function _move_rigidly!(frames::Union{Tuple, AbstractVector}, parts::Union{Tuple, AbstractVector}, offset)
    foreach(frame -> position!(frame, position(frame) .+ offset), frames)
    foreach(part -> translate3d!(part, offset), parts)
    return nothing
end

function _move_rigidly!(frames::Union{Tuple, AbstractVector}, parts::Union{Tuple, AbstractVector},
                        R::AbstractMatrix, pivot::AbstractVector)
    for frame in frames
        p = position(frame)
        orientation!(frame, R * orientation(frame))
        position!(frame, pivot + R * (p - pivot))
    end
    foreach(part -> rotate3d!(part, R, pivot), parts)
    return nothing
end

static_error(x) = throw(ArgumentError(lazy"$(typeof(x)) is static and cannot be moved, see `kinematic_trait_of`"))

# Entry points

"""
    translate3d!(x, offset)

Moves `x` by the `offset` vector. See [`BeamletOptics.Movable`](@ref).
"""
translate3d!(x, offset) = translate3d!(kinematic_trait_of(x), x, offset)

"""
    rotate3d!(x, R::AbstractMatrix)
    rotate3d!(x, axis::AbstractVector, θ::Real)
    rotate3d!(x, R::AbstractMatrix, pivot::AbstractVector)
    rotate3d!(x, axis::AbstractVector, θ::Real, pivot::AbstractVector)

Rotates `x` by the rotation matrix `R` (or by the angle `θ` [rad] about `axis`) about its own
`position`, or about `pivot` if given. See [`BeamletOptics.Movable`](@ref).
"""
rotate3d!(x, R::AbstractMatrix) = rotate3d!(kinematic_trait_of(x), x, R)

rotate3d!(x, axis::AbstractVector, θ::Real) = rotate3d!(kinematic_trait_of(x), x, axis, θ)

rotate3d!(x, R::AbstractMatrix, pivot::AbstractVector) = rotate3d!(kinematic_trait_of(x), x, R, pivot)

function rotate3d!(x, axis::AbstractVector, θ::Real, pivot::AbstractVector)
    return rotate3d!(kinematic_trait_of(x), x, axis, θ, pivot)
end

"""
    translate_to3d!(x, target)

Moves `x` such that its `position` coincides with `target`. See [`BeamletOptics.Movable`](@ref).
"""
translate_to3d!(x, target) = translate_to3d!(kinematic_trait_of(x), x, target)

"""
    xrotate3d!(x, θ)

Rotates `x` by the angle `θ` [rad] about the global x-axis through its `position`. See [`BeamletOptics.Movable`](@ref).
"""
xrotate3d!(x, θ) = xrotate3d!(kinematic_trait_of(x), x, θ)

"""
    yrotate3d!(x, θ)

Rotates `x` by the angle `θ` [rad] about the global y-axis through its `position`. See [`BeamletOptics.Movable`](@ref).
"""
yrotate3d!(x, θ) = yrotate3d!(kinematic_trait_of(x), x, θ)

"""
    zrotate3d!(x, θ)

Rotates `x` by the angle `θ` [rad] about the global z-axis through its `position`. See [`BeamletOptics.Movable`](@ref).
"""
zrotate3d!(x, θ) = zrotate3d!(kinematic_trait_of(x), x, θ)

"""
    align3d!(x, target)

Rotates `x` about its `position` such that its [`direction`](@ref) is aligned with `target`.
See [`BeamletOptics.Movable`](@ref).
"""
align3d!(x, target) = align3d!(kinematic_trait_of(x), x, target)

"""
    reset_translation3d!(x)

Moves `x` such that its `position` is the global origin. The result is exact, because `p + (-p)`
is zero in IEEE arithmetic. See [`BeamletOptics.Movable`](@ref).
"""
reset_translation3d!(x) = reset_translation3d!(kinematic_trait_of(x), x)

"""
    reset_rotation3d!(x)

Rotates `x` about its `position` such that its [`orientation`](@ref) is the identity. The
inverse of the orientation is its transpose, which is exact for all angles, including π. Only
available for [`BeamletOptics.Oriented`](@ref) values, see [`BeamletOptics.Movable`](@ref).
"""
reset_rotation3d!(x) = reset_rotation3d!(kinematic_trait_of(x), x)

"""
    direction(x)

Returns the direction of `x`. For [`BeamletOptics.Oriented`](@ref) values this is the local
y-axis, i.e. `orientation(x)[:, 2]`; rays and beams return their own direction.
"""
direction(x) = direction(kinematic_trait_of(x), x)

# Static branch: every kinematic functions throws

translate3d!(::Static, x, offset) = static_error(x)
rotate3d!(::Static, x, R::AbstractMatrix) = static_error(x)
rotate3d!(::Static, x, axis::AbstractVector, θ::Real) = static_error(x)
rotate3d!(::Static, x, R::AbstractMatrix, pivot::AbstractVector) = static_error(x)
rotate3d!(::Static, x, axis::AbstractVector, θ::Real, pivot::AbstractVector) = static_error(x)
translate_to3d!(::Static, x, target) = static_error(x)
xrotate3d!(::Static, x, θ) = static_error(x)
yrotate3d!(::Static, x, θ) = static_error(x)
zrotate3d!(::Static, x, θ) = static_error(x)
align3d!(::Static, x, target) = static_error(x)
reset_translation3d!(::Static, x) = static_error(x)
reset_rotation3d!(::Static, x) = static_error(x)
direction(::Static, x) = static_error(x)

# Movable branch: primitive fallbacks, must be implemented per type

translate3d!(::Movable, x, offset) = error(lazy"translate3d!(::Movable, x, offset) not defined for $(typeof(x))")

rotate3d!(::Movable, x, R::AbstractMatrix) = error(lazy"rotate3d!(::Movable, x, R::AbstractMatrix) not defined for $(typeof(x))")

# Movable branch: generic derived verbs, inner calls use the public entry points

function translate_to3d!(::Movable, x, target)
    translate3d!(x, target - position(x))
    return nothing
end

function rotate3d!(::Movable, x, axis::AbstractVector, θ::Real)
    rotate3d!(x, rotate3d(axis, θ))
    return nothing
end

function xrotate3d!(::Movable, x, θ)
    T = eltype(position(x))
    rotate3d!(x, Point3{T}(one(T), zero(T), zero(T)), θ)
    return nothing
end

function yrotate3d!(::Movable, x, θ)
    T = eltype(position(x))
    rotate3d!(x, Point3{T}(zero(T), one(T), zero(T)), θ)
    return nothing
end

function zrotate3d!(::Movable, x, θ)
    T = eltype(position(x))
    rotate3d!(x, Point3{T}(zero(T), zero(T), one(T)), θ)
    return nothing
end

function align3d!(::Movable, x, target)
    rotate3d!(x, align3d(direction(x), target))
    return nothing
end

function reset_translation3d!(::Movable, x)
    # Translate by the negative position
    translate3d!(x, -position(x))
    return nothing
end

function reset_rotation3d!(::Movable{Oriented}, x)
    # Rotate by the inverse orientation
    rotate3d!(x, transpose(orientation(x)))
    # Remove precision artifacts
    orientation!(x, one(orientation(x)))
    return nothing
end

function reset_rotation3d!(::Movable{Directed}, x)
    throw(ArgumentError(lazy"$(typeof(x)) has a direction but no orientation; use align3d! instead"))
end

function rotate3d!(::Movable, x, R::AbstractMatrix, pivot::AbstractVector)
    v = position(x) - pivot
    rotate3d!(x, R)
    translate3d!(x, R * v - v)
    return nothing
end

function rotate3d!(::Movable, x, axis::AbstractVector, θ::Real, pivot::AbstractVector)
    rotate3d!(x, rotate3d(axis, θ), pivot)
    return nothing
end

direction(::Movable{Oriented}, x) = Point3(orientation(x)[:, 2])
