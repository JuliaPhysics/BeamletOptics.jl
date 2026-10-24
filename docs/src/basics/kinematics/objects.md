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
| moved as one unit with a fixed pivot and named fields, e.g. a cemented doublet or a cube beamsplitter | an [`BeamletOptics.AbstractObject`](@ref) with the [`BeamletOptics.MultiShape`](@ref) trait whose parts are objects, like the [`CubeBeamsplitter`](@ref) | none: the solver hands every hit to the part that was hit |

The solver finds by itself which of several parts a ray hits next, also where parts touch: at a shared face it collects the surfaces of all touching parts and derives the media on both sides from them (see [Transitions](@ref)), e.g. a ray passes from the front into the back lens of a cemented doublet, or from one of two touching plates into the other. An own subtype of `AbstractObjectGroup` adds named parts and an own constructor to a group and needs no `interact3d` method; its fields and its constructor check are listed in the docstring of [`BeamletOptics.AbstractObjectGroup`](@ref), and the agent skill has a complete example in `templates/08_custom_group.jl`. A `MultiShape` object of plain shapes (no objects) needs an `interact3d` method; without one, it stops every ray that hits it, with the warning `No interact3d method defined`.
