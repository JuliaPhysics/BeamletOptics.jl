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
set_pivot3d!(::ObjectGroup, ::Any)
```
