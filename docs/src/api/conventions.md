# Conventions

In order to ensure implicit and explicit compliance with large parts of the API, the following conventions need to be used. Below you can find some essential conventions that need to be followed. For specific components or parts of the code base, refer to the documentation if certain guidelines need to be followed.

!!! warning
    Failure to comply with the following conventions can lead to spurious effects and silent bugs when using the API of this package!

## Global optical axis

Commonly, the z-axis is used as the principal optical axis when defining equations or alignment of optical systems and models. **This is not the case for BMO**, which uses the global y-axis in positive direction as the "global optical axis" in which effects are described. This is motivated by the plotting axes of [Makie](https://docs.makie.org/stable/), which uses a coordinate system where the x- and y-axis form the horizontal plane and the z-axis is orthogonal (upwards).

!!! info
    For (optical) equations that are depended on a global coordinate system, use a basis where the **global y-axis (`[0,1,0]`)** is the optical axis and the x-z-axes form the transverse plane. Global propagation along this axis is defined with a direction vector of `[0,1,0]`.

## Right-handedness

All **coordinate systems** are or must be defined **right-handed**! All normal vectors are or must be defined right-handed! 

## Counter-clockwise rotation

All **rotations** are or must be performed in a **counter-clockwise** manner for a positive rotation angle ``\theta > 0`` and vice-versa!
For a definition of rotation matrix order, refer to this [article](https://dominicplein.medium.com/extrinsic-intrinsic-rotation-do-i-multiply-from-right-or-left-357c38c1abfd).