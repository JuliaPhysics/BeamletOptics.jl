"""
    AbstractBeamGroup

Provides a generic container type interface for bundles of [`Beam`](@ref)s.
This interface assumes that there exists a central beam around which the bundle propagates,
e.g. akin to an optical axis.

# AbstractBeamGroup implementation reqs.

Subtypes of `AbstractBeamGroup` must implement the following:

## Fields:

- `beams`: a vector or tuple of root [`AbstractBeam`](@ref)s, e.g. [`Beam`](@ref)s or [`AstigmaticGaussianBeamlet`](@ref)s
- `center`: a `Point3{T}` which is regarded as the source position, i.e. the reference origin (pivot) of the group
- `orientation`: a `SMatrix{3,3,T,9}` that describes the local coordinate system of the group

The columns of `orientation` are the local x-, y- and z-axis. They must form a right-handed
orthonormal basis, like the `orientation` of an [`ObjectGroup`](@ref). The local y-axis
(second column) is the central source direction, i.e. the optical axis of the group. The
local x-axis (first column) is the reference vector of the azimuthal sampling (the `basis`
of the source constructors). Unlike a direction vector, the orientation therefore also
tracks a roll of the group about its own optical axis.

Since the kinematic API modifies `center` and `orientation`, the subtype must be a `mutable struct`.

## Functions:

If the fields above do not exist, the following getters/setters must be dispatched:

- `beams`: getter for the `beams` field or equivalent return type
- `position` / `position!`: gets or sets the source position (pivot)
- [`orientation`](@ref) / `orientation!`: gets or sets the orientation matrix
- `wavelength`: getter for the common wavelength of the beam bundle

The central source direction is derived from the orientation via [`direction`](@ref).

## Kinematic

An `AbstractBeamGroup` is a container: it takes the kinematic class of its `beams`, i.e.
[`BeamletOptics.Movable`](@ref) with an [`BeamletOptics.Oriented`](@ref) frame for movable beams, see
[`BeamletOptics.AbstractKinematicTrait`](@ref). The constructors check that the `beams` are either all
static or all movable. The following logic is applied to

- [`translate3d!`](@ref): all beams and the group `center` are translated by the offset vector
- [`rotate3d!`](@ref): all beams are rotated around the `center` point with respect to their relative position, `orientation` is rotated
- [`reset_translation3d!`](@ref) / [`reset_rotation3d!`](@ref): moves the group back to the origin or its standard orientation
  (identity: the central source direction/optical axis is the global +y-axis and the azimuthal sampling
  reference vector `basis` is the global +x-axis); the beams are always reset, even if the group is already at the identity
- [`set_pivot3d!`](@ref): moves the group `center` (the pivot used above) without moving or resetting the `beams`

Except for `set_pivot3d!`, every command resets each beam to its untraced start state.
"""
abstract type AbstractBeamGroup{T <: Real, R <: AbstractRay{T}} end

kinematic_trait_of(bg::AbstractBeamGroup) = _container_trait(beams(bg))

beams(bg::AbstractBeamGroup) = bg.beams

Base.length(bg::AbstractBeamGroup) = length(beams(bg))
Base.iterate(bg::AbstractBeamGroup, state...) = iterate(beams(bg), state...)
Base.getindex(bg::AbstractBeamGroup, i::Int) = getindex(beams(bg), i)

Base.position(bg::AbstractBeamGroup) = bg.center
position!(bg::AbstractBeamGroup{T}, pos) where {T} = (bg.center = Point3{T}(pos))

"""
    orientation(bg::AbstractBeamGroup) -> SMatrix{3,3}

Returns the orientation matrix of the beam group `bg`. Its columns are the local x-axis (the
azimuthal sampling reference vector), the local y-axis (the central source direction, see
[`direction`](@ref)) and the local z-axis, forming a right-handed orthonormal basis.

The orientation tracks the full rotational state of the group, including a roll about its
own optical axis. It is updated by [`rotate3d!`](@ref) and returned to the identity by
[`reset_rotation3d!`](@ref).
"""
orientation(bg::AbstractBeamGroup) = bg.orientation

"""
    orientation!(bg::AbstractBeamGroup, M)

Overwrites the orientation matrix of the beam group `bg` with `M`, see [`orientation`](@ref).

!!! warning
    This only sets the bookkeeping of the group, the `beams` are **not** moved. Use
    [`rotate3d!`](@ref) to rotate a group including its beams. `M` is not validated and must
    be a right-handed orthonormal matrix.
"""
orientation!(bg::AbstractBeamGroup{T}, M) where {T} = (bg.orientation = SMatrix{3, 3, T, 9}(M))

wavelength(bg::AbstractBeamGroup) = wavelength(first(rays(first(beams(bg)))))

"""
    translate3d!(::Movable, bg::AbstractBeamGroup, offset)

Moves all beams of the group and its `center` by `offset`. Every beam is reset to its untraced start state.
"""
function translate3d!(::Movable, bg::AbstractBeamGroup, offset)
    foreach(b -> translate3d!(b, offset), beams(bg))
    position!(bg, position(bg) + offset)
    return nothing
end

"""
    rotate3d!(::Movable, bg::AbstractBeamGroup, R::AbstractMatrix)

Rotates all beams of the group by `R` about the group `center` and updates the group
[`orientation`](@ref) to `R * orientation(bg)`. Every beam is reset to its untraced start state.
"""
function rotate3d!(::Movable, bg::AbstractBeamGroup, R::AbstractMatrix)
    p = position(bg)
    foreach(b -> rotate3d!(b, R, p), beams(bg))
    orientation!(bg, R * orientation(bg))
    return nothing
end

"""
    set_pivot3d!(bg::AbstractBeamGroup, pivot)

Moves the kinematic pivot (`center`) of the beam group `bg` to `pivot`, without moving or
resetting any of its `beams`. The pivot is the reference point used by [`rotate3d!`](@ref)
and [`reset_translation3d!`](@ref); the group [`orientation`](@ref) is unchanged.

Unlike [`translate3d!`](@ref) or [`rotate3d!`](@ref), this only changes bookkeeping and does
not reset already traced beams.
"""
function set_pivot3d!(bg::AbstractBeamGroup, pivot)
    position!(bg, pivot)
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", bg::AbstractBeamGroup)
    println(io, "Subtype of AbstractBeamGroup")
    println(io, "   # of beams: $(length(beams(bg)))")
    return nothing
end
