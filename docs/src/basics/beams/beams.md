```@setup beams
beam_showcase_dir = joinpath(@__DIR__, "..", "..", "assets", "beam_renders")

Main.DocUtils.conditional_include(joinpath(beam_showcase_dir, "beam_showcase.jl"))
```

# Beams

As mentioned in the [Rays](@ref) section, a beam within the context of this package serves as a data structure for storing collections of rays, forming the backbone of the simulation framework. Beams are intended to be designed as [AbstractTrees](https://github.com/JuliaCollections/AbstractTrees.jl) to allow for ray bifurcations, e.g. in the case of optical elements such as beamsplitters. The [`solve_system!`](@ref) function relies on this data structure to perform ray tracing computations within optical systems. 

To ensure compatibility and extensibility, beam types must adhere to the [`BeamletOptics.AbstractBeam`](@ref) interface. Refer to its documentation for more information.

## Basic beam

A minimal implementation of the [`BeamletOptics.AbstractBeam`](@ref) type is provided by the [`Beam`](@ref). It can be used to store a light path through an optical system. If the beam is split, its children will be recursively traced until all paths are solved.

```@docs; canonical=false
Beam
```

A ray tracing example through an arbitrary system using a [`Beam`](@ref) is shown below. Individual [`Ray`](@ref) segments are marked by their starting position and direction. The [Beam expander](@ref) and [Miniature microscope](@ref) tutorial covers the use of the [`Beam`](@ref) in more detail. 

![Beam structure](beam_showcase.png)

## Moving sources

[`BeamletOptics.AbstractRay`](@ref)s, [`BeamletOptics.AbstractBeam`](@ref)s (i.e. [`Beam`](@ref), [`GaussianBeamlet`](@ref) and [`AstigmaticGaussianBeamlet`](@ref)) and [`BeamletOptics.AbstractBeamGroup`](@ref)s can be repositioned with the same kinematic verbs used for optical elements (see [Moving optical elements](@ref)):

- Translation
    - [`translate3d!`](@ref)
    - [`translate_to3d!`](@ref)
- Rotation
    - [`rotate3d!`](@ref)
    - [`xrotate3d!`](@ref)
    - [`yrotate3d!`](@ref)
    - [`zrotate3d!`](@ref)
    - [`align3d!`](@ref)
- Reset commands
    - [`reset_translation3d!`](@ref)
    - [`reset_rotation3d!`](@ref) (beam groups only)
- Pivot commands (beam groups only)
    - [`set_pivot3d!`](@ref)

!!! warning "Every move resets the source"
    A traced beam/beamlet whose start is moved would otherwise keep stale rays and child beams from the old, now geometrically wrong, light path. To avoid this, **every** call to one of the translation, rotation or reset commands above first resets the moved beam (or, for a group, every one of its beams) to its untraced start state: all rays beyond the start ray (of every component beam of a beamlet) are dropped, and any child beams (e.g. created by a beamsplitter interaction) are removed. Call [`solve_system!`](@ref) again after moving a source to retrace it. [`set_pivot3d!`](@ref) is the exception, see below.

The pivot used for rotations is the *start* of the beam, i.e. the position of its first/chief ray (the chief ray of a [`GaussianBeamlet`](@ref) or [`AstigmaticGaussianBeamlet`](@ref)). For a [`BeamletOptics.AbstractBeamGroup`](@ref), the pivot is instead the group's `center`, so that all member beams are rotated rigidly about the source origin. Rotating a group also updates its [`orientation`](@ref), so a roll about the group's own optical axis is tracked. [`reset_translation3d!`](@ref) moves the `position` of any of these sources — a ray, a beam or a beam group's `center` — back to the global origin. [`reset_rotation3d!`](@ref) only works for a [`BeamletOptics.AbstractBeamGroup`](@ref), moving it back to `orientation` = identity (optical axis along +y); rays and single beams have no orientation, only a direction, so [`reset_rotation3d!`](@ref) throws an `ArgumentError` for them — use [`align3d!`](@ref) instead. [`set_pivot3d!`](@ref) moves a group's `center` to a new pivot without moving or resetting its beams, e.g. to `rotate3d!` a group about one of its own beams instead of its default center.

!!! info "Polarized rays"
    Rotating a [`PolarizedRay`](@ref) (including the chief ray of an [`AstigmaticGaussianBeamlet`](@ref)) also rotates its field vector `E0` together with the direction, so it stays orthogonal to the new propagation direction.

Only a root beam can be moved. Calling a kinematic verb on a child beam, e.g. one created by a beamsplitter interaction during tracing, throws an `ArgumentError` — move the root beam instead.

```julia
using BeamletOptics

beam = Beam([0, 0, 0], [0, 1, 0])
translate3d!(beam, [1, 0, 0])     # shift the start by [1, 0, 0], resets the beam
zrotate3d!(beam, deg2rad(90))     # rotate 90° about z through the (new) start position
```