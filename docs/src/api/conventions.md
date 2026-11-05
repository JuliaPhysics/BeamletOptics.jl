# Conventions

In order to ensure implicit and explicit compliance with large parts of the API, the following conventions need to be used unless specified otherwise. Below you can find some essential conventions. For specific components or parts of the code base, refer to the respective documentation if certain guidelines need to be followed.

!!! warning
    Failure to comply with the following conventions can lead to spurious effects and silent bugs when using the API of this package!

## Global optical axis

Commonly, the z-axis is used as the principal optical axis when defining equations or alignment of optical systems and models. **This is not the case for BMO**, which uses the global y-axis in positive direction as the "global optical axis" in which effects are described. This is motivated by the plotting axes of [Makie](https://docs.makie.org/stable/), which uses a coordinate system where the x- and y-axis form the horizontal plane and the z-axis is orthogonal (upwards).

!!! info
    For (optical) equations that are depended on a global coordinate system, use a basis where the **global y-axis (`[0,1,0]`)** is the optical axis and the x-z-axes form the transverse plane. Global propagation along this axis is defined with a direction vector of `[0,1,0]`.

## Right-handedness

All **coordinate systems** are or must be defined **right-handed**! All normal vectors are or must be defined right-handed! 

## Counter-clockwise rotation

All **rotations** are or must be performed in a **counter-clockwise** manner for a positive rotation angle ``\theta > 0`` and vice-versa! This is the right-hand rule: a positive ``\theta`` about an axis turns counter-clockwise when seen from the tip of the axis vector towards its origin, e.g. a rotation about `[0,0,1]` turns the x-axis towards the y-axis.

Rotations are **extrinsic**: the `axis` of [`rotate3d!`](@ref) is given in global coordinates, not in the local frame of the rotated object, and the rotation acts about the global axis through the pivot (by default the `position` of the object). A rotation matrix ``R`` is applied from the left, i.e. the new orientation is `orientation = R * orientation`, so that a sequence of rotations is multiplied from right to left, the first rotation being the rightmost factor. For the difference between extrinsic and intrinsic rotations, refer to this [article](https://dominicplein.medium.com/extrinsic-intrinsic-rotation-do-i-multiply-from-right-or-left-357c38c1abfd).

## Normal vector direction

For a closed volume BMO assumes that the surface normal vector points outside of the volume.

The normal of a face of a [`BeamletOptics.Mesh`](@ref) follows from the order of its vertices by the right-hand rule: the normal is the cross product of the edges from the first vertex to the second and to the third. The vertices of each face of a closed mesh must therefore run counter-clockwise when seen from outside of the volume.

!!! info "Normal vector direction definition"
    Equations to calculate optical effects often rely on the normal vector at the ray intersection location to work correctly and point in a specific direction.
    It is important to ensure that this condition is fulfilled when spurious effects occur.