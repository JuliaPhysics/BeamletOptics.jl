# Conventions: units, coordinates, signs

## Units

Everything is SI. The package defines no unit helpers except `BeamletOptics.inch` (= 25.4e-3, not exported).
Define your own at the top of a script:

```julia
const mm = 1e-3
const µm = 1e-6
const nm = 1e-9
```

| Quantity            | Unit        | Example                         |
|---------------------|-------------|---------------------------------|
| lengths, radii      | m           | `25.4mm`, `0.1`                 |
| wavelength `λ`      | m           | `532e-9`, `1064nm`              |
| angles              | rad         | `deg2rad(45)`                   |
| power `P0`          | W           | `P0 = 1e-3`                     |
| field `E0`          | V/m         | `E0 = [1, 0, 0]`                |
| intensity           | W/m²        | output of `intensity(detector)` |

Defaults (`src/Config.jl`): `get_default_wavelength()` = 1e-6 m, `get_default_waist()` = 1e-3 m,
`get_default_power()` = 1e-3 W, `get_default_r_max()` = 100, `get_default_depth_max()` = 100.

## Coordinate system

- Right-handed global frame. **+y is the default optical axis**; x and z span the transverse plane,
  **z points up** (as in Makie).
- Positive rotation angles rotate counter-clockwise (right-hand rule).
- Components spawn **at the origin, aligned with +y**:
  - lenses: front vertex at y = 0, extending to y = `thickness(lens)`
  - plano mirrors: reflecting face at y = 0, normal along y, substrate towards +y
  - `Detector`, `PolarizationFilter`: zero-thickness square/disc in the x-z plane at y = 0
  - `CubeBeamsplitter`: centered at the origin, coating at 45° to the y-axis
  - curved mirrors (`ParabolicMirror`, `ConicMirror`, …): vertex at the origin, opening towards −y
- A typical source starts before the first element, e.g. `Beam([0, -0.05, 0], [0, 1, 0], λ)`.
- `orientation(obj)` returns a 3×3 matrix whose columns are the local x, y, z axes;
  `direction(obj) == orientation(obj)[:, 2]` is the local optical axis.
- Surface normals of closed volumes point outward.

## Kinematics

| Function                               | Semantics                                                    |
|----------------------------------------|--------------------------------------------------------------|
| `translate3d!(x, Δ)`                   | **relative** move by `Δ`                                     |
| `translate_to3d!(x, p)`                | absolute: `position(x) == p` afterwards                      |
| `xrotate3d!(x, θ)` / `y…` / `z…`       | rotate by θ about the global x/y/z axis through `position(x)`|
| `rotate3d!(x, axis, θ[, pivot])`       | rotate about an arbitrary axis, optionally about `pivot`     |
| `rotate3d!(x, R[, pivot])`             | apply rotation matrix `R`                                    |
| `align3d!(x, target_dir)`              | rotate so that `direction(x)` points along `target_dir`      |
| `reset_translation3d!(x)`              | move back to the origin                                      |
| `reset_rotation3d!(x)`                 | back to identity orientation (objects/groups only, not beams)|
| `set_pivot3d!(group, p)`               | move the pivot of an `ObjectGroup`/beam group, not its members |

All rotations are cumulative. To fold a +y beam by 90° with a plano mirror at the origin,
rotate the mirror by ±45°: `zrotate3d!(m, deg2rad(45))` folds into the x-y plane,
`xrotate3d!(m, deg2rad(-45))` folds towards +z.

Moving a source (`Ray`, `Beam`, beamlet, beam group) **resets** its traced path. Only root beams can be
moved, not children created by splitters.

## Sign conventions

- **Lens radii (ISO 10110)**: R > 0 if the center of curvature lies at +y of the surface (to the right),
  R < 0 if at −y, `Inf` for plano. Biconvex: `r1 > 0, r2 < 0`. Plano-convex with the curved side first:
  `r1 > 0, r2 = Inf`.
- `BeamletOptics.lensmakers_eq(R1, R2, n)` returns the thin-lens focal length **f** (not 1/f).
- **Conic mirrors**: `R > 0` is concave, opening towards −y; `k = -1` parabola, `k = 0` sphere.
  Ellipsoid: `s`, `s′` same sign; hyperboloid: opposite signs.
- **Even asphere coefficients** `α` enter as `Σ α[i] r^(2i)`, so `α[1]` is the r² term, `α[2]` the r⁴ term, …,
  in SI units. Convert mm-based catalog values `A_2i` with `A_2i * (1e3)^(2i-1)`:
  `[0, A4*1e3^3, A6*1e3^5, A8*1e3^7, ...]`.
- **Gaussian beam**: `w0` is the waist *radius* (1/e² intensity); `gauss_parameters` returns
  curvature `R` as 1/r (not a radius) and Gouy phase ψ = −atan(z/z_R).
- **Beamsplitter** `reflectance` keyword is the reflected **power** fraction (0.7 → 70:30, default 0.5 → 50:50);
  internally the amplitudes are `r = √reflectance`, `t = √(1 − r²)`. The reflection phase depends on the
  beam type and the incidence side (see `components/beamsplitters.md`). Children are ordered `[transmitted, reflected]`.
- **Detector local coordinates** `(x, z)` form a left-handed frame with the detector normal (the normal
  points against the incoming beam, initially −y).
