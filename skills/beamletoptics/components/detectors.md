# Detectors

```julia
Detector(edge_length, stop = true)
```

- An infinitely thin square screen at the origin in the x-z plane. Its normal points along −y, so beams
  travelling along +y hit it. Orient it like any other object, e.g. after `zrotate3d!(det, deg2rad(90))`
  it catches beams travelling along −x.
- With `stop = false`, beams continue behind the detector.
- Stores *hits* of exactly **one** type (Ray, PolarizedRay, GaussianBeamlet or AstigmaticGaussianBeamlet).
- **Accumulates across solves:** call `empty!(det)` before each `solve_system!` that reuses it.
- Read out all results before moving the detector again.

## Readout

| Call | Returns |
|------|---------|
| `spot_diagram(det)` | `Vector{Point2}` of local `(x, z)` hit points (beamlets: projected 1/e² contour, kwargs `num_spots`, `crop_factor`) |
| `electric_field(det; kwargs...)` | `(xs, zs, E)`: complex field on an `n×n` grid (polarized hits: a matrix of 3-vectors) |
| `intensity(det[, Z]; kwargs...)` | `(xs, zs, I)` with I = \|E\|²/(2Z) in W/m² |
| `optical_power(det; kwargs...)` | integrated power in W |

Grid keywords:

- `n = 100`: grid points per axis
- `crop_factor = 1`: enlarges the automatic window around the hits
- `center = Centroid()` or `MinMax()`: how the automatic window is centered
- `x_min, x_max, z_min, z_max`: explicit window in local coordinates (m)
- `x0_shift, z0_shift`: shift the window
- `num_spots`: beamlet contour resolution used for the automatic limits

Local `(x, z)` form a left-handed frame with the normal, as seen by the incoming beam.
For polarized hits, `I = intensity.(E)` converts a matrix of field vectors to intensity.
Ray-based PSFs are not normalized: compare shapes and positions, not absolute Strehl values.

If nothing reaches the detector, readout throws "No hits available on detector."

## Fourier Beamlet Detector

```julia
FourierBeamletPropagator(shape; stop = true, is_planar = false)
FourierBeamletPropagator(edge_length; stop = true, is_planar = false)
FourierBeamletPropagator(x_edge, z_edge; stop = true, is_planar = true)
FourierBeamletPropagator(dx, dy, dz; stop = true, is_planar = false)
FourierBeamletPropagator(object; stop = true, is_planar = false)
```

- High-performance wave-optical detector using accelerated Non-Uniform Fast Fourier Transforms (`NonuniformFFTs.jl`) to synthesize coherent continuous fields without mesh discretization artifacts.
- Converts incident rays and beamlets (`Ray`, `PolarizedRay`, `GaussianBeamlet`, `AstigmaticGaussianBeamlet`) into spectral Fourier modes (`FourierKPoint`) preserving phase, optical path length (OPL), and 3D polarization vectors.
- Supports planar 2D aperture surfaces (`is_planar = true`), volumetric cuboids (`dx, dy, dz`), and general 3D shapes.
- **Accumulates across solves:** call `empty!(det)` before each `solve_system!` that reuses it.

### Readout & Field Synthesis

| Call | Returns |
|------|---------|
| `synthesize_field(det, grid)` | Complex 3D electric field array `E` on a `SpatialGrid` via NUFFT (V/m) |
| `calculate_observables(det, (xs, zs))` | `DetectorObservables` struct containing `intensity`, `phase`, `stokes`, and `power` |
| `calculate_observables!(obs, det, (xs, zs))` | In-place zero-allocation evaluation of optical observables |
| `generate_fourier_spectrum(det)` | `Vector{FourierKPoint}` containing all discrete plane-wave modes |

