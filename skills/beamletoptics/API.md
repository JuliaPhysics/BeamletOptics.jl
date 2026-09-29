# BeamletOptics API primer

BeamletOptics (often aliased `const BMO = BeamletOptics`) follows one pattern:

```julia
using BeamletOptics
const mm = 1e-3

# 1. build components (at the origin, aligned with +y)
lens = SphericalLens(50mm, -50mm, 5mm, 25.4mm, λ -> 1.5)
det  = Detector(20mm)

# 2. place them (relative moves, radians)
translate3d!(det, [0, 55mm, 0])

# 3. collect them in a system (order irrelevant, tracing is non-sequential)
system = System([lens, det])

# 4. create a source and solve
beam = Beam([0, -20mm, 0], [0, 1, 0], 633e-9)
solve_system!(system, beam)

# 5. evaluate
rays(beam)             # traced ray segments
spot_diagram(det)      # hits on the detector
```

Names not exported must be qualified: `BeamletOptics.inch`, `BeamletOptics.lensmakers_eq`, …
Use `scripts/api_lookup.jl NAME` to print the docstring and signatures of any name.

## Exported names (by category)

| Category      | Names |
|---------------|-------|
| Kinematics    | `translate3d!`, `translate_to3d!`, `rotate3d!`, `xrotate3d!`, `yrotate3d!`, `zrotate3d!`, `align3d!`, `reset_translation3d!`, `reset_rotation3d!`, `set_pivot3d!`, `position`, `direction`, `orientation` |
| Rays & beams  | `Ray`, `PolarizedRay`, `Beam`, `GaussianBeamlet`, `AstigmaticGaussianBeamlet`, `rays`, `point_on_beam`, `rayleigh_range`, `normal3d` |
| Sources       | `CollimatedSource`, `UniformDiscSource`, `PointSource`, `UniformPointSource`, `set_num_rays!`, `CollimatedGaussianBeamletSource`, `SphericalGaussianBeamletSource`, `EllipticalGaussianBeamletSource`, `GaussianBeamletDecomposition`, `WavefrontBeamletDecomposition`, `AstigmaticBeamGroup` |
| System        | `System`, `StaticSystem`, `solve_system!`, `ObjectGroup` |
| Inspection    | `properties`, `default_properties` |
| Materials     | `DiscreteRefractiveIndex`, `SellmeierEquation` |
| Lenses        | `Lens`, `ThinLens`, `SphericalLens`, `DoubletLens`, `SphericalDoubletLens`, `TripletLens`, `SphericalTripletLens`, `thickness` |
| Surfaces      | `SphericalSurface`, `CircularFlatSurface`, `RectangularFlatSurface`, `EvenAsphericalSurface`, `CylindricalSurface`, `AcylindricalSurface` |
| Mirrors       | `Mirror`, `RoundPlanoMirror`, `SquarePlanoMirror`, `SquarePlanoMirror2D`, `RectangularPlanoMirror`, `RightAnglePrismMirror`, `SphericalMirror`, `ConicMirror`, `OffAxisConicMirror`, `ParabolicMirror`, `OffAxisParabolicMirror`, `EllipsoidalMirror`, `OffAxisEllipsoidalMirror`, `HyperbolicMirror`, `OffAxisHyperbolicMirror`, `Retroreflector` |
| Splitters     | `ThinBeamsplitter`, `RoundThinBeamsplitter`, `RectangularPlateBeamsplitter`, `RoundPlateBeamsplitter`, `CubeBeamsplitter`, `RectangularCompensatorPlate` |
| Prisms        | `Prism`, `RightAnglePrism` |
| Polarizers    | `PolarizationFilter`, `RoundPolarizationFilter`, `LinearPolarizer`, `RoundLinearPolarizer`, `transmission_axis` |
| Detectors     | `Detector`, `spot_diagram`, `intensity`, `electric_field`, `optical_power`, `gauss_parameters`, `waist_parameters`, `Centroid`, `MinMax` |
| Dummies       | `MeshDummy`, `NonInteractableObject`, `IntersectableObject` |
| Config        | `get_default_wavelength`, `get_default_waist`, `get_default_power`, `get_default_r_max`, `get_default_depth_max`, `get_invariant_threshold`, `set_invariant_threshold!`, `get_internal_reflection_threshold`, `get_line_plane_intersection_threshold`, `get_orthogonality_threshold`, `get_sdf_surface_threshold`, `get_sdf_raymarch_eps`, `get_sdf_inside_step`, `get_progress_threshold`, `set_progress_threshold!` |
| Render (Makie)| `render!`, `get_view`, `set_view`, `hide_axis`, `set_orthographic`, `look_at!`, `arrow!`, `render_lcs!`, `set_render_look`, `studio_lighting!` |
| Live (Makie)  | `live_render!`, `update_render!`, `remove_render!`, `pick_object` |

Useful non-exported helpers: `BeamletOptics.inch`, `lensmakers_eq(R1, R2, n)` (returns f),
`divergence_angle(λ, w0, M2)`, `numerical_aperture(θ, n=1)`, `optical_path_length(beam)`,
`isparaxial(system, beam, θ=π/4)`, `fresnel_coefficients`, `beams(group)`, `objects(system)`,
`list_subtypes(T)`, constants `Z_vacuum`, `c_vacuum`.

## Properties of objects

`properties(x)` lists what an object, shape, beam or source stores as `name => value` pairs, e.g.
`"Thickness [m]" => 0.005`, `"Reflectance" => 0.3` (power ratio), `"Hits" => 12`. Values are SI
numbers; the unit is in brackets at the end of the name, names without brackets are dimensionless.
It lists only stored values: e.g. the
radii and diameter of a `Lens` are not stored and not listed.

A custom `AbstractObject` subtype gets the default list (type, position, optical axis, shape). To
show its own parameters, extend it:

```julia
BeamletOptics.properties(x::MyFilter) =
    [default_properties(x); "Optical density" => x.od; "Center wavelength [m]" => x.λc]
```

## Solving

```julia
solve_system!(system, beam; r_max = 100, retrace = true, depth_max = 100,
              check_invariant = true, threshold = get_invariant_threshold())
solve_system!(system, beam_group; progress = true, kwargs...)   # multithreaded over member beams
```

- `r_max`: max. rays per beam leaf (raise it for resonators, e.g. facing mirrors)
- `depth_max`: max. splitting depth of the beam tree
- `retrace`: after the first solve, try to retrace the previous path sequentially and fall back to a
  full non-sequential solve where it breaks. Pass `retrace = false` if an element was moved *into*
  the existing beam path (e.g. a chopper), which otherwise fails silently.
- Julia threads (`julia -t auto`) speed up solving of sources with many beams.
- `progress`: beam groups and detector readout show a progress bar once they have run for
  `get_progress_threshold()` s (default 5 s, `set_progress_threshold!(Inf)` disables it). It is
  only drawn if `stderr` is a terminal, so scripts with piped output see nothing. In the live
  view of BeamletOpticsGUI, the same loops show a progress window in the 3D view instead.

## Beam queries

| Call                              | Returns |
|-----------------------------------|---------|
| `rays(beam)`                      | vector of traced ray segments of this beam (not children) |
| `position(ray)`, `direction(ray)` | start point and unit direction of a segment |
| `length(beam)`                    | geometric path length up to the last intersection |
| `BeamletOptics.optical_path_length(beam)` | optical path length |
| `point_on_beam(beam, t)`          | `(point, ray_index)` at distance `t` along the beam |
| `beam.children`                   | sub-beams after a split, `[transmitted, reflected]` |
| `last(rays(beam)).E0`             | field vector of a polarized beam's last segment |
| `gauss_parameters(g, z)`          | `(w, R, ψ, w0)` for `GaussianBeamlet`; `z` scalar or vector |
| `gauss_parameters(agb, z)`        | `(w1, w2, R1, R2, ψ, w01, w02)` for `AstigmaticGaussianBeamlet` |
| `waist_parameters(g_or_agb, z)`   | local waist geometry including 3D semi-axis vectors |
| `rayleigh_range(g; M2 = 1)`       | Rayleigh range of the first section (M2 is not stored, pass it again) |
| `optical_power(g_or_agb)`         | beamlet power in W |
| `electric_field(g, r, z)`, `intensity(agb, r, z)` | field / intensity at radius `r`, distance `z` |

`z` for Gaussian queries is the distance along the (unfolded) beam path from its start point.

## Detector readout

```julia
empty!(det)                                   # before each solve that reuses det
solve_system!(system, source)
pts        = spot_diagram(det)                # Vector{Point2}: local (x, z) in m
x, z, E    = electric_field(det; n = 200)     # complex field on an n×n grid
x, z, I    = intensity(det; n = 200)          # W/m²
P          = optical_power(det)               # W
```

Grid keywords for `electric_field`/`intensity`/`optical_power`: `n`, `crop_factor`, `x_min`, `x_max`,
`z_min`, `z_max`, `x0_shift`, `z0_shift`, `center = Centroid() | MinMax()`, and for beamlets `num_spots`.
They also take `progress = true` (see Solving). See `components/detectors.md`.

## Makie rendering

`using GLMakie` before or after `using BeamletOptics` activates the extension.
One generic function draws everything: `render!(ax, system)`, `render!(ax, beam)`,
`render!(ax, source)`, `render!(ax, object)`. `live_render!` updates the plots
of moved components and re-solved beams in place. See `VISUALIZATION.md`.
