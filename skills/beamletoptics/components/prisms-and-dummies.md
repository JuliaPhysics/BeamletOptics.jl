# Prisms, retroreflectors & dummies

## Prisms

- `RightAnglePrism(leg_length, height, n)`: refractive right-angle prism with legs in x/y and height in z.
  It is not aligned with the y-axis at spawn, so render it once to check its orientation.
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
