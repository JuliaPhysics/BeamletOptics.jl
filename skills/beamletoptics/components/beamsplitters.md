# Beamsplitters & compensators

The `reflectance` keyword is the reflected **power** fraction: `0.7` gives a 70:30 splitter, the default
`0.5` a 50:50 splitter. Internally the field amplitudes are `r = √reflectance` and `t = √(1 − r²)`.

The reflection phase depends on which side of the coating is hit, the same for every beam type that
carries a field: a beam that faces the coating normal (`dot(dir, normal) < 0`) is reflected with an
extra π, one from the other side without (Stokes relations, r' = −r). Only then do the two output
ports of an interferometer sum to its input power.

| Beam type | Reflected field |
|-----------|-----------------|
| `Ray` (`Beam` of rays) | no field is carried, so no phase |
| `PolarizedRay` | Jones matrix `SPBasis(-r, 0, 0, r)` from the back, `SPBasis(r, 0, 0, -r)` when facing the normal |
| `GaussianBeamlet` | scalar amplitude `r` from the back, `-r` when facing the normal |
| `AstigmaticGaussianBeamlet` | chief polarization as for `PolarizedRay` |

Transmitted beams get the amplitude `t` and no phase shift. Behind the coating of a plate beamsplitter the
field of a polarized beam is scaled for the glass (the coating transmits the power fraction `t²`), and it
leaves the uncoated back with the Fresnel loss of that surface. In interferometers, only relative phases
between arms are meaningful. Do not add phase offsets by hand to "correct" for the splitter; if a
fringe offset matters, check which side of the coating each arm reflects from.

Splitting creates a beam tree: `beam.children == [transmitted, reflected]`.

| Constructor | Geometry at spawn |
|-------------|-------------------|
| `CubeBeamsplitter(leg_length, n; reflectance = 0.5)` | centered at the origin, coating at 45° to the y-axis |
| `PolarizingCubeBeamsplitter(leg_length, n; Rs = 1, Rp = 0)` | as built: a +y beam transmits its x-component (p) along +y and reflects its z-component (s) along −x; `Rp > 0` for a finite extinction ratio |
| `RectangularPlateBeamsplitter(width, height, thickness, n; reflectance = 0.5)` | coating centered at the origin, substrate towards +y (from y = 0 to `thickness`) |
| `RoundPlateBeamsplitter(diameter, thickness, n; reflectance = 0.5)` | as above, round |
| `ThinBeamsplitter(width[, height]; reflectance = 0.5)` | zero-thickness coating (testing, composites) |
| `RoundThinBeamsplitter(diameter; reflectance = 0.5)` | round, zero thickness |
| `RectangularCompensatorPlate(width, height, thickness, n)` | uncoated plate, first surface at the origin, aligned with +y |

Plate beamsplitters are usually rotated by 45° (`zrotate3d!(bs, deg2rad(45))`) and paired with a
compensator plate in interferometers. A cube needs no compensator. See `templates/04_michelson_scan.jl`
for a working cube setup (`zrotate3d!(cbs, deg2rad(-90))` puts the reflected arm on +x for a +y input).
