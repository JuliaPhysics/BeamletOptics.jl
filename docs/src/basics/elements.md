# Optical elements

Optical elements serve as the building blocks for optical systems in the context of this package, representing components such as mirrors, lenses, filters, etc. 

## Types of elements

Some optical elements are provided with this package, these include e.g.:

- Reflective optical elements
    - [`RoundPlanoMirror`](@ref)
- Refractive optical elements
    - [`SphericalLens`](@ref)
- Misc.
    - [`Photodetector`](@ref)
    - [`ThinBeamsplitter`](@ref)

For a detailed overview, refer to the [Optical components](@ref) section.

!!! info "Custom optical elements"
    In order to implement custom geometries and optical elements, visit the [API design](@ref) section.

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

For the easier representation of a group of [`BeamletOptics.AbstractObject`](@ref)s that moves as one, the [`ObjectGroup`](@ref) can be used. Refer to the [Lens groups](@ref) example for more information.

```@docs; canonical=false
ObjectGroup
```

