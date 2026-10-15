# Moving optical elements

Optical elements that implement the [`BeamletOptics.AbstractObject`](@ref) interface, as well as the [`BeamletOptics.AbstractShape`](@ref)s they are built from, can be moved with the commands listed in the [Kinematics](@ref) section. When objects are manipulated, they are translated and rotated around their self-defined center of gravity, i.e. their [`position`](@ref). Their [`orientation`](@ref) defines their local fixed coordinate system.

```julia
using BeamletOptics

mirror = RoundPlanoMirror(25.4e-3, 5e-3)
translate3d!(mirror, [0, 0.1, 0])     # move 10 cm along the global y-axis
zrotate3d!(mirror, deg2rad(45))       # rotate 45° about z through the mirror position
```

## Objects with multiple shapes

Objects which consist of more than one shape, e.g. a [`CubeBeamsplitter`](@ref) or a doublet lens, move all of their shapes as one rigid body. Unless the object defines its own reference point, the position and orientation of its **first** shape are used as the kinematic center, see [`BeamletOptics.MultiShape`](@ref).

## Groups of optical elements

For the easier representation of a group of [`BeamletOptics.AbstractObject`](@ref)s that moves as one, the [`ObjectGroup`](@ref) can be used. Refer to the [Lens groups](@ref) example for more information. The group's pivot (its `center`) can be moved without moving its members via [`set_pivot3d!`](@ref).

```@docs; canonical=false
ObjectGroup
set_pivot3d!(::BeamletOptics.AbstractObjectGroup, ::Any)
```

## Assemblies: group or multi-shape object

An assembly such as a telescope, a periscope or an interferometer arm can be modeled in two ways, and the choice follows from whether its parts touch:

| Parts of the assembly | Model | `interact3d` |
|---|---|---|
| separate, with air between them | an [`ObjectGroup`](@ref) or an own subtype of [`BeamletOptics.AbstractObjectGroup`](@ref) | none: each part is traced and interacts on its own |
| touching and interacting as one, e.g. cemented or with a shared face | an [`BeamletOptics.AbstractObject`](@ref) with the [`BeamletOptics.MultiShape`](@ref) trait, like the [`CubeBeamsplitter`](@ref) | required: the solver hands every hit to the object, which forwards it to the part that was hit |

The solver can find out by itself which of several separate parts a ray hits next. Parts that touch need a [`BeamletOptics.Hint`](@ref) across their shared face, which only the component knows. An own subtype of `AbstractObjectGroup` adds named parts and an own constructor to a group and needs no `interact3d` method; its fields and its constructor check are listed in the docstring of [`BeamletOptics.AbstractObjectGroup`](@ref), and the agent skill has a complete example in `templates/08_custom_group.jl`. A `MultiShape` object without an `interact3d` method stops every ray that hits it, with the warning `No interact3d method defined for the MultiShape object`.
