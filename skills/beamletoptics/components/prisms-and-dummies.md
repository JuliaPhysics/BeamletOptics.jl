# Prisms, retroreflectors & dummies

## Prisms

All prisms are extruded along local z (`height`, centered at z = 0). A face "faces a beam along +y" if its
outward normal is -y. `translate3d!` and the rotations act on the origin given below.

- `RightAnglePrism(leg_length, height, n)`: refractive right-angle prism. The legs have the normals -x and -y,
  the hypotenuse has the normal `(1, 1, 0)/√2`. **The origin is the center of the hypotenuse face**, not the
  right-angle edge, which lies at `(-leg_length/2, -leg_length/2, 0)`.
  - As spawned, the leg with normal -y faces a beam along +y: the beam is totally internally reflected at the
    hypotenuse (needs `n > √2`) and leaves along -x. After `yrotate3d!(prism, π)` it leaves along +x.
  - `zrotate3d!(prism, 5π/4)` (225°) makes the hypotenuse face a beam along +y: the beam is reflected at both legs
    and returned antiparallel, mirrored about the line through the origin.
- `Prism(vertices, height, n)`: right prism over any strictly convex polygon. `vertices` are `(x, y)` points in [m] in the local
  frame (either orientation, validated), extruded by `height` along local z, origin at the polygon coordinates' origin.
  A face faces a beam along +y if its edge's outward normal is -y; otherwise orient it with the kinematic functions.
- `EquilateralPrism(side, height, n)`: dispersing prism, origin at the centroid, apex along +y, base normal -y.
  A symmetric (minimum deviation) ray runs parallel to the base, entering and leaving the faces next to the apex, and is deflected towards the base.
  For a beam along +y, `zrotate3d!(prism, π/2 - θ)` with `θ = asin(n * sin(π/6)) - π/6` (index at the design
  wavelength) sets minimum deviation; the beam is deflected towards +x.
- `DovePrism(length, aperture, height, n)`: trapezoid with 45° ends, long axis along y, TIR base face has normal -x,
  origin at the center of the cross-section's bounding box. As spawned, an end face faces a beam along +y.
  The beam leaves parallel to +y, but in general laterally shifted. Requires `length > 2 * aperture > 0`.
- `Prism(shape, n)`: generic refractive body from an SDF or mesh shape (advanced).

Dispersion needs a λ-dependent index (`SellmeierEquation`). Trace one `Beam` per wavelength.

## Retroreflector

`Retroreflector(scale)`: corner cube; `scale = 1e-3` makes it mm-sized.

## Mechanics and obstacles

| Constructor | Behavior |
|-------------|----------|
| `NonInteractableObject(path_to_stl)` | visual only, rays pass through (mounts, housings) |
| `NonInteractableObject(shape)` | visual only, generic shape |
| `IntersectableObject(path_to_stl)` | hard stop: rays terminate on it (apertures, beam dumps) |

STL files are usually in mm; scale and translate them to match the SI scene (see the package tutorials,
e.g. `docs/src/tutorials/michelson.md`). To move a dummy together with its optic, put both in an `ObjectGroup`.
