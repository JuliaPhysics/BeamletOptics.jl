# Detectors

```julia
Detector(edge_length, stop = true)
```

- An infinitely thin square screen at the origin in the x-z plane. Its normal points along −y, so beams
  travelling along +y hit it. Orient it like any other object, e.g. after `zrotate3d!(det, deg2rad(90))`
  it catches beams travelling along −x.
- With `stop = false`, beams continue behind the detector.
- Stores *hits* of exactly **one** type (Ray, PolarizedRay, GaussianBeamlet or AstigmaticGaussianBeamlet).
- **Hits are values** of the beam state at the solve. Solving again, or moving a source, mirror or lens, never changes stored hits.
- **Moving the detector empties it:** every kinematic call on it (also a move of a group containing it) clears its hits.
- **Accumulates across solves:** by default the solver does not reset detectors, so hits of consecutive solves sit side by side and all beams of a beam group, and several sources solved one after another, superpose on purpose. The user decides when to evaluate. To start empty, solve with `initialize = true`, or call `initialize!(system)` (all detectors of a system, nested groups included) or `empty!(det)` (one detector).
- `spot_diagram` shows all stored hits. `electric_field`, `intensity` and `optical_power` superpose **all** stored hits **coherently**, as for beams that exist at the same time. For solves that follow each other in time, know what this computes, or evaluate after each solve and reset (`initialize = true`).

## Readout

| Call | Returns |
|------|---------|
| `spot_diagram(det)` | `Vector{Point2}` of local `(x, z)` hit points (beamlets: projected 1/e² contour, kwargs `num_spots`, `crop_factor`) |
| `electric_field(det; kwargs...)` | `(xs, zs, E)`: complex field on an `n×n` grid (hits of `PolarizedRay`s and `AstigmaticGaussianBeamlet`s: a matrix of complex 3-vectors, use `sum(abs2, e)` per element) |
| `intensity(det[, Z]; kwargs...)` | `(xs, zs, I)` with I = \|E\|²/(2Z) in W/m² |
| `optical_power(det; kwargs...)` | integrated power in W |

Grid keywords:

- `n = 100`: grid points per axis
- `crop_factor`: scales the automatic window around the hits; default `1` for rays (box around the hit points) and `3` for beamlets (three beam radii around each hit)
- `center = Centroid()` or `MinMax()`: how the automatic window is centered
- `x_min, x_max, z_min, z_max`: explicit window in local coordinates (m)
- `x0_shift, z0_shift`: shift the window
- `num_spots`: beamlet contour resolution used for the automatic limits

`progress = true` shows a progress bar in the terminal once the readout has run for
`get_progress_threshold()` s (default 5 s).

Local `(x, z)` form a left-handed frame with the normal, as seen by the incoming beam.
For polarized hits, `I = intensity.(E)` converts a matrix of field vectors to intensity.
Ray-based PSFs are not normalized: compare shapes and positions, not absolute Strehl values.

If nothing reaches the detector, readout throws "No hits available on detector."
