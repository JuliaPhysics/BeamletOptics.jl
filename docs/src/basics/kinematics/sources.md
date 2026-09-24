# Moving sources

[`Ray`](@ref)s, [`Beam`](@ref)s, beamlets like the [`GaussianBeamlet`](@ref) and groups of beams can be repositioned with the commands listed in the [Kinematics](@ref) section.

!!! warning "Every move resets the source"
    A traced beam/beamlet whose start is moved would otherwise keep stale rays and child beams from the old, now geometrically wrong, light path. To avoid this, **every** call to one of the translation, rotation or reset commands first resets the moved beam (or, for a group, every one of its beams) to its untraced start state: all rays beyond the start ray (of every component beam of a beamlet) are dropped, and any child beams (e.g. created by a beamsplitter interaction) are removed. Call [`solve_system!`](@ref) again after moving a source to retrace it. [`set_pivot3d!`](@ref) is the exception, see below.

```julia
using BeamletOptics

beam = Beam([0, 0, 0], [0, 1, 0])
translate3d!(beam, [1, 0, 0])     # shift the start by [1, 0, 0], resets the beam
zrotate3d!(beam, deg2rad(90))     # rotate 90° about z through the (new) start position
```

## Rays and beams

The pivot used for rotations is the *start* of the beam, i.e. the position of its first/chief ray (the chief ray of a [`GaussianBeamlet`](@ref) or [`AstigmaticGaussianBeamlet`](@ref)). [`reset_translation3d!`](@ref) moves this start back to the global origin. Rays and single beams have no orientation, only a direction, so [`reset_rotation3d!`](@ref) throws an `ArgumentError` for them — use [`align3d!`](@ref) instead.

!!! info "Polarized rays"
    Rotating a [`PolarizedRay`](@ref) (including the chief ray of an [`AstigmaticGaussianBeamlet`](@ref)) also rotates its field vector `E0` together with the direction, so it stays orthogonal to the new propagation direction.

Only a root beam can be moved. Calling a kinematic command on a child beam, e.g. one created by a beamsplitter interaction during tracing, throws an `ArgumentError` — move the root beam instead.

## Beam groups

For a [`BeamletOptics.AbstractBeamGroup`](@ref), the pivot is the group's `center`, so that all member beams are rotated rigidly about the source origin. Rotating a group also updates its [`orientation`](@ref), so a roll about the group's own optical axis is tracked. [`reset_translation3d!`](@ref) moves the `center` back to the global origin, and [`reset_rotation3d!`](@ref) moves the group back to `orientation` = identity (optical axis along +y). Refer to [`BeamletOptics.AbstractBeamGroup`](@ref) for the orientation convention.

[`set_pivot3d!`](@ref) moves a group's `center` to a new pivot without moving or resetting its beams, e.g. to [`rotate3d!`](@ref) a group about one of its own beams instead of its default center.

```@docs; canonical=false
set_pivot3d!(::BeamletOptics.AbstractBeamGroup, ::Any)
```
