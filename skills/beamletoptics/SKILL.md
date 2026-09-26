---
name: beamletoptics
description: Build, modify, and debug optical simulations with BeamletOptics.jl (Julia, "BMO"). Use when writing Julia code that uses BeamletOptics, or when the user asks to simulate or ray-trace lenses, mirrors, beamsplitters, polarizers, detectors, Gaussian beams/beamlets, interferometers, PSFs, spot diagrams, or to render an optical setup with Makie.
allowed-tools: Read, Write, Edit, Grep, Glob, Bash
---

You are helping the user build optical simulations with **BeamletOptics.jl**: a Julia package for
non-sequential 3D ray tracing and Gaussian beamlet propagation, with a Makie extension for rendering.

When this Skill is active:

- Use only documented BeamletOptics API. **Do not invent constructors, keyword arguments or
  functions.** If unsure, look it up before writing code:
  - `julia --project=<env> scripts/api_lookup.jl <Name> [<Name>...]` prints docstrings and method signatures.
  - `julia --project=<env> scripts/api_lookup.jl --exports` lists every exported name.
  - Read `src/` of the installed package (`julia -e 'using BeamletOptics; println(pkgdir(BeamletOptics))'`).
- Everything is **SI**: meters, radians, watts, V/m. **Wavelengths are in meters** (`532e-9`).
- The **optical axis is +y**, not z. Components spawn at the origin aligned with +y.
- **Refractive indices are callables** `n(λ)` (λ in meters), never plain numbers:
  `λ -> 1.5`, `DiscreteRefractiveIndex(...)` or `SellmeierEquation(...)`.
- **Detectors accumulate.** Call `empty!(detector)` before every `solve_system!` that reuses it.
- **All kinematics are relative** (`translate3d!`, `xrotate3d!`, …), angles in radians (`deg2rad`).
  Use `translate_to3d!` for absolute positioning.
- Run the script after writing it. A simulation that "looks right" but was never executed is not done.

## Default workflow

1) Clarify the optical problem (if not already given)
- Wavelength(s), source type (collimated/diverging, rays vs. Gaussian beam), beam size, power
- Components with catalog data (radii, center thickness, diameter, glass), distances between them
- What to evaluate: spot size, focus position, beam radius, power, fringe pattern, PSF, polarization, a picture

2) Pick the beam model (see `components/beams-and-sources.md`)
- Geometric rays / spot diagrams / aberrations → `Beam`, `CollimatedSource`, `PointSource`
- PSF / coherent intensity from rays → `UniformDiscSource` (equal-area sampling) + `Detector`
- Paraxial laser beam on-axis (waist, Rayleigh range, interferometer power) → `GaussianBeamlet`
- Tilted, off-axis, cylindrical, or polarization-sensitive Gaussian beam → `AstigmaticGaussianBeamlet`
- Polarization with rays → `Beam(pos, dir, λ, E0)` (a `PolarizedRay` beam)

3) Build components at the origin, then place them
- Construct each element (`components/*.md`), then orient with `xrotate3d!/yrotate3d!/zrotate3d!` and
  place with `translate3d!`/`translate_to3d!`. Rotations pivot about `position(obj)`, so for a single
  object the order of rotate/translate does not matter; pass a `pivot` to `rotate3d!` to rotate about another point.
- Stack lenses along +y using `thickness(lens)`.
- Group rigid sub-assemblies with `ObjectGroup` (see `components/systems-and-groups.md`).

4) Solve
- `system = System([obj1, obj2, ...])` (order does not matter, tracing is non-sequential);
  `StaticSystem` for small fixed systems that are solved many times.
- `solve_system!(system, beam_or_source)`.
- For parameter scans, mutate the scene in a loop, `empty!` detectors, and re-solve (retracing is automatic).

5) Evaluate
- Rays: `rays(beam)`, `position`, `direction`, `length(beam)`, `beam.children` (splitters: transmitted first).
- Detector: `spot_diagram`, `intensity`, `electric_field`, `optical_power` (`components/detectors.md`).
- Gaussian: `gauss_parameters(beam, z)`, `waist_parameters`, `rayleigh_range`, `optical_power`.
- Sanity-check numbers against a paraxial estimate (`BeamletOptics.lensmakers_eq`, Airy radius, Malus, …).

6) Visualize (optional, needs a Makie backend) → `VISUALIZATION.md`
- `render!(ax, system)`, `render!(ax, beam)`; CairoMakie for headless PNGs, GLMakie for interactive 3D.

## Safety and non-goals

- BeamletOptics does not model coatings, scattering or stray light from refractive surfaces (polarized
  refraction traces only the refracted ray). Results are simulations, not certified optical designs;
  point out model limits when they matter for the user's question.
- Do not modify the BeamletOptics package source unless the user asks for it; write user scripts instead.

## Local references bundled with this Skill

- API primer: `API.md`
- Units, coordinates, sign conventions: `CONVENTIONS.md`
- Workflow patterns (scans, focus search, interferometers, groups): `WORKFLOW.md`
- Rendering with Makie: `VISUALIZATION.md`
- Pre-delivery checklist and common pitfalls: `CHECKLIST.md`
- Runnable reference scripts: `templates/`
- Helper scripts: `scripts/` (`api_lookup.jl`, `run_templates.jl`)

## Components

- [Beams, rays & sources](./components/beams-and-sources.md)
- [Refractive index / materials](./components/materials.md)
- [Lenses & surfaces](./components/lenses.md)
- [Mirrors](./components/mirrors.md)
- [Beamsplitters & compensators](./components/beamsplitters.md)
- [Polarizers](./components/polarizers.md)
- [Detectors](./components/detectors.md)
- [Prisms, retroreflectors & dummies](./components/prisms-and-dummies.md)
- [Systems & object groups](./components/systems-and-groups.md)
