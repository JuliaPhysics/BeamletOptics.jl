# Beams, rays & sources

All positions and directions are 3-vectors in m; directions are normalized for you. λ is in m.

## Rays and ray beams

| Constructor | Notes |
|-------------|-------|
| `Ray(pos, dir, λ = 1e-6)` | single geometric ray, starts in n = 1 |
| `PolarizedRay(pos, dir, λ = 1e-6, E0 = [1, 0, 0])` | `E0` complex 3-vector (V/m), **must be ⟂ `dir`** |
| `Beam(pos, dir, λ = 1e-6)` | beam of `Ray`s: the usual thing to trace |
| `Beam(pos, dir, λ, E0)` | beam of `PolarizedRay`s |
| `Beam(ray)` | wrap an existing ray |

A `Beam` is a tree: `rays(beam)` are its segments, `beam.children` the sub-beams created by
beamsplitters (`[transmitted, reflected]`), `beam.parent` the parent. Only root beams can be moved;
moving resets the trace. `empty!(beam)` resets it manually.

## Gaussian beamlets

```julia
GaussianBeamlet(pos, dir, λ = 1e-6, w0 = 1e-3; M2 = 1, P0 = 1e-3, z0 = 0, support = nothing)
```

- Stigmatic TEM00 beam represented by 3 rays (chief, waist, divergence). `w0` = waist **radius**,
  `z0` = waist offset along the axis, `P0` in W, `M2` beam quality factor.
- Valid for on-axis, untilted, (nearly) aberration-free systems only.
- Queries: `gauss_parameters(g, z)` → `(w, R, ψ, w0)`, `waist_parameters`, `rayleigh_range(g; M2)`,
  `optical_power(g)`, `electric_field(g, r, z)`, `point_on_beam(g, t)`.

```julia
AstigmaticGaussianBeamlet(pos, dir, λ, w0; M2, P0, E0, support, z0)
AstigmaticGaussianBeamlet(pos, dir, λ, w0_x, w0_y; M2_x, M2_y, P0, E0, support, z0, z0_x, z0_y)
```

- General astigmatic beamlet (9 rays, polarized chief ray). Use it for tilted, off-axis, cylindrical
  or polarization-dependent setups. `E0` defaults to a field along the `support` axis scaled to `P0`.
- `solve_system!` checks the optical invariant and stops tracing if it is violated.
- Queries: `gauss_parameters(agb, z)` → `(w1, w2, R1, R2, ψ, w01, w02)`, `rayleigh_range(agb)` →
  `(z_rx, z_ry)`, `intensity(agb, r, z)`, `electric_field`, `BeamletOptics.polarized_field`.

## Beam groups (sources)

All sources are iterable collections of beams (`BeamletOptics.beams(src)`). `solve_system!(system, src)`
traces them multithreaded, and `render!(ax, src; render_every = 5)` draws them.
Kinematics act on the whole group; `set_pivot3d!` moves its pivot.

| Constructor | Samples |
|-------------|---------|
| `CollimatedSource(pos, dir, diameter, λ = 1e-6; num_rings = 10, num_rays = 100num_rings, basis)` | parallel rays on concentric rings (use `num_rays ≥ 20 num_rings`) |
| `UniformDiscSource(pos, dir, diameter, λ = 1e-6; num_rays = 1000, basis)` | parallel rays, equal-area (Fibonacci) sampling. **Use for PSF/intensity.** |
| `PointSource(pos, dir, θ, λ = 1e-6; num_rings = 10, num_rays, basis)` | diverging rays, half-angle `θ` (rad, < π) |
| `UniformPointSource(pos, dir, θ, λ = 1e-6; num_rays = 1000, basis)` | diverging rays, equal solid angle |
| `CollimatedSource(beams, diameter, pos, dir)` / `PointSource(beams, NA, pos, dir)` | wrap your own beams (e.g. polarized ones) |
| `CollimatedGaussianBeamletSource(pos, dir, D, λ, w0s; n_grid = 20, basis, randomize_axes, rng)` | square grid of astigmatic beamlets, `w0s ≈ D/n_grid` |
| `SphericalGaussianBeamletSource(pos, dir, θ, λ; num_rings, num_rays, overlap = 1.2, basis, randomize_axes, rng, P0, E0)` | diverging beamlet fan |
| `EllipticalGaussianBeamletSource(pos, dir, θ_x, θ_y, λ; …)` | elliptical diverging beamlet fan |
| `GaussianBeamletDecomposition(pos, dir, λ, w0; n_grid, overlap, basis, randomize_axes, rng, P0, E0, threshold = 1e-4)` | one large Gaussian beam decomposed into beamlets |
| `WavefrontBeamletDecomposition(x, y, amplitude, phase, dir, λ; threshold, overlap, basis, randomize_axes, rng, E0)` | arbitrary sampled wavefront; the grid is centered at the origin, move it with `translate3d!` |
| `AstigmaticBeamGroup(beams, pos, dir_or_orientation)` | wrap your own astigmatic beamlets |

Trace beamlet groups onto a `Detector` and read them out with `intensity`/`electric_field`, which sum
the beamlets coherently.
