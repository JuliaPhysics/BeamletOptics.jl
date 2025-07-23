# Optical elements

Optical elements serve as the building blocks for optical systems in the context of this package, representing components such as mirrors, lenses, filters, etc. Unlike the surface/interface based representation of optical elements in other tools, they are treated as volumetric bodies in this simulation framework. Optical interactions between rays/beams and elements are defined based on the type of the element and the type of the incident beam/ray. **Note that optical elements will be referred to simply as objects moving forward**. A more detailed look into the design philosophy behind this approach is given in the [Objects and shapes](@ref) section below.

To ensure compatibility with the [API design](@ref), custom optical elements must adhere to the [`BeamletOptics.AbstractObject`](@ref) interface.

!!! info "Normal vector direction definition"
    Equations to calculate optical effects often rely on the normal vector at the ray intersection location to work correctly and point in a specific direction.
    It is important to ensure that this condition is fulfilled when spurious effects occur.

## Objects and shapes

In BMO, the distinction between an *object* and its geometric representation (*shape*) is a central design principle. This separation is intended to ensure flexibility and modularity in modeling optical components.

### Separation of geometry and optical interactions

Objects combine physical geometry with specific optical interactions. The geometry, represented by an [`BeamletOptics.AbstractShape`](@ref), defines the physical boundaries of the element. Shapes can be represented in various forms, such as [Meshes](@ref) or [Signed Distance Functions (SDFs)](@ref). For more information, refer to the respective documentation.

On the other hand, the optical behavior — how light interacts with the element — is defined by the [`BeamletOptics.AbstractObject`](@ref) type. This decoupling allows for independent development and extension of geometry representations and optical interaction models.

### Multi-shape objects

An [`BeamletOptics.AbstractObject`](@ref) can consist of multiple [`BeamletOptics.AbstractShape`](@ref)s or even multiple subsidiary [`BeamletOptics.AbstractObject`](@ref)s, facilitating the creation of composite optical elements. For example, a lens with an anti-reflective coating could be represented as the substrate and a seperate model for the coating, each with its own geometric and optical properties. In general, an `object` can be a [`BeamletOptics.SingleShape`](@ref) or a [`BeamletOptics.MultiShape`](@ref). Refer to the [Geometry representation](@ref) section for more information. 

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

