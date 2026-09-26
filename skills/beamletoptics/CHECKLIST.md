# Checklist before handing over a simulation

## Correctness

- [ ] All lengths and wavelengths in **meters** (`532e-9`, not `532`); angles in **radians** (`deg2rad`).
- [ ] Refractive indices are callables of λ (`λ -> 1.5`, `SellmeierEquation`, `DiscreteRefractiveIndex`).
- [ ] `DiscreteRefractiveIndex` contains the exact wavelength of every source (no interpolation).
- [ ] Lens radius signs follow ISO 10110 (biconvex: `r1 > 0`, `r2 < 0`).
- [ ] Geometry assumes +y as the optical axis; rotations are cumulative and relative.
- [ ] Sources start **outside** of all elements and point at them.
- [ ] `empty!(detector)` before every re-solve; results read out before the detector is moved.
- [ ] Beam model fits the question (see below).
- [ ] Script was **run**, and the printed numbers were compared with a paraxial or analytic estimate.
- [ ] For non-trivial geometry, a rendered image was checked.

## Choosing the beam model

| Question                                   | Use |
|--------------------------------------------|-----|
| ray paths, spot diagram, aberrations       | `Beam`, `CollimatedSource`, `PointSource` |
| PSF / intensity from rays                  | `UniformDiscSource`, `UniformPointSource` |
| on-axis laser beam: waist, z_R, power      | `GaussianBeamlet` |
| tilted/off-axis/cylindrical Gaussian beam  | `AstigmaticGaussianBeamlet` |
| extended Gaussian/diverging wavefronts     | `CollimatedGaussianBeamletSource`, `SphericalGaussianBeamletSource`, `GaussianBeamletDecomposition` |
| polarization with rays                     | `Beam(pos, dir, λ, E0)` (E0 ⟂ dir) |

## Common pitfalls

1. **Detectors accumulate** hits across solves. Forgetting `empty!(det)` silently sums results.
2. **One hit type per detector**: don't send `Beam`s and `GaussianBeamlet`s to the same detector.
3. **Blocked polarized rays terminate** at the filter: then `last(rays(beam))` is the *incident* ray.
   Check `length(rays(beam))` or read power at a detector instead.
4. **`GaussianBeamlet` is stigmatic only.** It is wrong for tilted or off-axis optics; use
   `AstigmaticGaussianBeamlet`. Beamlets must be much smaller than the optics and must not be clipped.
5. **`lensmakers_eq` returns f**, not 1/f; it is a thin-lens estimate. For thick lenses the focus is
   near the back focal length (plano-convex, curved side first: `BFL ≈ f − ct/n`).
6. **`SphericalLens` with `l = 0`** builds a thin lens via `ThinLens`. If the center thickness is too small for the
   radii and diameter, the constructor throws ("cylinder section length ≤ 0"): increase `l` or reduce `d`.
7. **`Mirror`s reflect on every face**, including the back and edges.
8. **Polarized refraction** traces only the refracted ray (no ghost reflections); no coatings are modeled.
9. **Doublets/triplets** assume flush cemented contact; only spherical menisci are supported;
   aspheres are experimental.
10. **Resonators** (facing mirrors, cavities) need a bounded `r_max`/`depth_max`.
11. **`ObjectGroup` rotates about its `center`**, which starts at the origin, not the centroid:
    use `set_pivot3d!`.
12. **`reset_rotation3d!` throws for rays/beams**; use `align3d!` there.
13. **Moving a child beam is not allowed**; move the source (this also resets its trace).
14. **`retrace = false`** is needed when a moved element newly blocks an old path segment.
15. **Rendering**: `flen` defaults to 1 m for rays; CairoMakie supports only `Axis3`; camera helpers need
    GLMakie `LScene`.
