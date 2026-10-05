# Prisms, retroreflectors & dummies

## Prisms

- `RightAnglePrism(leg_length, height, n)`: refractive right-angle prism with legs in x/y and height in z.
  `n` may be a number for a constant refractive index. It is not aligned with the y-axis at spawn, so render it once to check its orientation.
- `Prism(vertices, height, n)`: right prism over any strictly convex polygon. `vertices` are `(x, y)` points in [m] in the local
  frame (either orientation, validated), extruded by `height` along local z, origin at the polygon coordinates' origin.
  Not aligned with the y-axis; orient it with the kinematic functions.
- `EquilateralPrism(side, height, n)`: dispersing prism, origin at the centroid, apex along +y, base normal -y.
  A symmetric (minimum deviation) ray runs parallel to the base, entering and leaving the faces next to the apex, and is deflected towards the base.
- `DovePrism(length, aperture, height, n)`: trapezoid with 45° ends, long axis along y, TIR base face has normal -x.
  A ray along its axis exits undeviated. Requires `length > 2 * aperture > 0`.
- `Prism(shape, n)`: generic refractive body from an SDF or mesh shape (advanced).

Dispersion needs a λ-dependent index (`SellmeierEquation`). Trace one `Beam` per wavelength.

## Retroreflector

`Retroreflector(scale)`: corner cube; `scale = 1e-3` makes it mm-sized.

## Mechanics and obstacles

| Constructor | Behavior |
|-------------|----------|
| `MeshDummy(path_to_stl)` | visual only, rays pass through (mounts, housings) |
| `NonInteractableObject(shape)` | visual only, generic shape |
| `IntersectableObject(path_to_stl)` | hard stop: rays terminate on it (apertures, beam dumps) |

STL files are usually in mm; scale and translate them to match the SI scene (see the package tutorials,
e.g. `docs/src/tutorials/michelson.md`). To move a dummy together with its optic, put both in an `ObjectGroup`.
