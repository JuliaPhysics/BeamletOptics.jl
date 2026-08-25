# Optical components

Optical elements serve as the building blocks for optical systems in the context of this package, representing components such as mirrors, lenses, filters and so on. A collection of basic optical elements is provided with this package as is. They are tested for the correctness of their optical interactions and are verified to work with reasonable fidelity. For detailed documentation, refer to the following table of contents. 

## Component overview

```@contents
Pages = ["mirrors.md", "lenses.md", "beamsplitters.md", "detectors.md", "polarizers.md"]
Depth = 2
```

## Moving optical elements

Optical elements can move around freely in three-dimensional space, which enables the modeling of kinematics within optical setups. When objects are manipulated, they are translated and rotated around their self-defined center of gravity, which is represented as a ``\mathbb{R}^3``-vector and will be referred to as its [`position`](@ref). Additionally, the [`orientation`](@ref) of an object, defined as its local fixed coordinate system, is represented by an orthonormal matrix in ``\mathbb{R}^3``. If the object is rotated, this matrix can be used to calculate the inverse transform into global coordinates. 

!!! important "Optical system kinematics"
    Elements can be moved freely between each call of [`solve_system!`](@ref). However, during tracing it is assumed that all elements remain static.

For elements that implement the [`BeamletOptics.AbstractObject`](@ref) interface, the following movement commands are provided:

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
    - [`reset_rotation3d!`](@ref)

!!! important "Relative motion"
    Unless specified otherwise, the translation and rotation commands result in relative motions to the current position and orientation. This must be taken into account when trying to model a specific set of movements.

## Groups of optical elements

For the easier representation of a group of [`BeamletOptics.AbstractObject`](@ref)s that moves as one, the [`ObjectGroup`](@ref) can be used. A group accepts the same movement commands listed above, but is purely a kinematic convenience — it has no optical identity of its own, see [`BeamletOptics.AbstractObjectGroup`](@ref) and the [Core design](@ref) chapter. Refer to the [Lens groups](@ref) example for more information.

```@docs; canonical=false
ObjectGroup
```

## Listing available components

When using this package in the REPL, a tree view of all implemented [`BeamletOptics.AbstractObject`](@ref)s can be generated via the [`BeamletOptics.list_subtypes`](@ref) helper function. Note that this function is not able to determine all available constructors.

```@repl
using BeamletOptics # hide
BeamletOptics.list_subtypes(BeamletOptics.AbstractObject);
```