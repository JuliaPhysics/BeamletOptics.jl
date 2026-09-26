# Workflow patterns

## Environment

- BeamletOptics requires Julia ≥ 1.12. Work in a project environment:
  `julia --project=. -e 'using Pkg; Pkg.add("BeamletOptics")'` (add `CairoMakie`/`GLMakie` for plots).
- Run scripts non-interactively with `julia --project=. script.jl`. First runs precompile and can take
  a minute; don't mistake that for a hang.
- Use `julia -t auto` when tracing sources with many beams (solving a beam group is multithreaded).
- Put reusable component factories in functions, e.g. `LA1805() = SphericalLens(15.5mm, Inf, 8.6mm, 1inch, NBK7)`.
  Every call creates a new, independent object. Objects are mutable and have exactly one pose, so
  create one instance per physical element instead of reusing a variable.

## Build → place → solve → evaluate

1. Construct every element at the origin.
2. Orient (`x/y/zrotate3d!`, `align3d!`), then position (`translate3d!`/`translate_to3d!`).
3. Put everything, including detectors and `MeshDummy` mounts, into one `System`.
4. Create the source *after* deciding where it starts; `solve_system!`.
5. Read out results; print numbers you can sanity-check.

## Stacking lenses along the axis

```julia
L1 = SphericalLens(...); L2 = SphericalLens(...)
gap = 2mm
translate3d!(L2, [0, thickness(L1) + gap, 0])
```

For long stacks keep a running `y` and `translate3d!(Li, [0, y, 0]); y += thickness(Li) + gap_i`
(see `docs/src/assets/examples/double_gauss.jl` in the package repo).

## Parameter scans (focus search, alignment, fringes)

Mutate the scene inside a loop and re-solve. Retracing is automatic and fast.

```julia
for y in ys
    translate_to3d!(det, [0, y, 0])   # absolute
    empty!(det)                       # REQUIRED: detectors accumulate hits
    solve_system!(system, source)
    push!(metric, f(spot_diagram(det)))
end
```

- Read everything you need from a detector **before** moving it again (hits reference beam state).
- If a moved element blocks a path segment that was previously free, pass `retrace = false`.
- See `templates/01_singlet_spot_diagram.jl` (focus scan) and `templates/04_michelson_scan.jl` (fringes).

## Interferometers and coherent sums

- Use `GaussianBeamlet` (on-axis, paraxial) or `AstigmaticGaussianBeamlet` (tilted/off-axis) sources and
  evaluate `optical_power(det)` or `intensity(det)`; phases add coherently on the detector.
- `DiscreteRefractiveIndex` values must exist for the exact source wavelength.
- A beamsplitter produces a beam tree: `beam.children == [transmitted, reflected]`.

## PSFs from rays

- Use `UniformDiscSource` (equal-area sampling) with thousands of rays and a `Detector` near focus;
  read `intensity(det; n, crop_factor)` or give explicit `x_min/x_max/z_min/z_max` windows.
- For unpolarized PSFs, trace two orthogonal polarizations (polarized beams wrapped in a
  `CollimatedSource(beams, diameter, pos, dir)`) and add the intensities incoherently.

## Sub-assemblies

```julia
arm = ObjectGroup([mirror, mount])    # mount = MeshDummy("mount.stl")
set_pivot3d!(arm, position(mirror))   # rotate about the mirror instead of the origin
translate_to3d!(arm, [0.1, 0, 0])
zrotate3d!(arm, deg2rad(90))
system = System([arm, ...])           # groups are flattened automatically
```

Groups can be nested; `translate_to3d!` moves the group `center` onto the target.

## Visual check

When geometry is non-trivial (folds, tilts, groups), render once (`templates/07_render_system.jl`)
and look at the image before trusting numbers. Axis-aligned top view: `Axis3(...; azimuth = -π/2, elevation = π/2)`.

## Debugging checklist

- Beam stops early: print `length(rays(beam))`, and for each `ray` its `position`/`direction`.
  Common causes are the element being too small, the source starting inside an element, a wrong rotation
  sign, a fully blocking polarizer, or `r_max` being reached.
- Gaussian beamlet stops: all 3 (9) rays must hit the same surfaces; the beamlet may be clipped or the
  optical invariant check failed (`check_invariant = false` to diagnose, not to "fix").
- `KeyError` from `DiscreteRefractiveIndex`: wavelength not in the table.
- `MissingBackendError` on `render!`: load `CairoMakie` or `GLMakie`.
- Detector error "No hits available": nothing reached it; check position/orientation (normal −y by default).
