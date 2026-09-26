# Beamsplitters & compensators

The `reflectance` keyword is the reflected **power** fraction: `0.7` gives a 70:30 splitter, the default
`0.5` a 50:50 splitter. Internally the field amplitudes are `r = √reflectance` and `t = √(1 − r²)`.
The reflected beam gets a π phase jump.

Splitting creates a beam tree: `beam.children == [transmitted, reflected]`.

| Constructor | Geometry at spawn |
|-------------|-------------------|
| `CubeBeamsplitter(leg_length, n; reflectance = 0.5)` | centered at the origin, coating at 45° to the y-axis |
| `RectangularPlateBeamsplitter(width, height, thickness, n; reflectance = 0.5)` | coating centered at the origin, substrate towards −y |
| `RoundPlateBeamsplitter(diameter, thickness, n; reflectance = 0.5)` | as above, round |
| `ThinBeamsplitter(width[, height]; reflectance = 0.5)` | zero-thickness coating (testing, composites) |
| `RoundThinBeamsplitter(diameter; reflectance = 0.5)` | round, zero thickness |
| `RectangularCompensatorPlate(width, height, thickness, n)` | uncoated plate, first surface at the origin, aligned with +y |

Plate beamsplitters are usually rotated by 45° (`zrotate3d!(bs, deg2rad(45))`) and paired with a
compensator plate in interferometers. A cube needs no compensator. See `templates/04_michelson_scan.jl`
for a working cube setup (`zrotate3d!(cbs, deg2rad(-90))` puts the reflected arm on +x for a +y input).
