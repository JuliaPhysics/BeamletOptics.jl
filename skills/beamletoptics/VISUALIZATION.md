# Visualization (Makie extension)

Rendering lives in the `BeamletOpticsMakieExt` extension and uses **GLMakie**:

```julia
using GLMakie, BeamletOptics      # Axis3 and LScene, camera helpers, save to PNG
```

Without a Makie backend, `render!` throws a `MissingBackendError`. GLMakie needs OpenGL: on a
headless Linux machine, run scripts under a virtual display, e.g. `xvfb-run -a julia --project=. script.jl`.

**Fallback without OpenGL or display** (CI containers, cloud sandboxes, no `xvfb`): CairoMakie also
activates the extension and saves PNG/PDF/SVG files. It supports `Axis3` only: no `LScene` and none
of the `LScene` camera helpers below. Swap the import and keep the rest of the script:

```julia
using CairoMakie, BeamletOptics   # Axis3 only; save("setup.png", fig)
```

Use GLMakie when it works and CairoMakie only when GLMakie cannot open a context.

## Minimal example

```julia
fig = Figure()
ax = Axis3(fig[1, 1], aspect = :data)   # or LScene(fig[1, 1])
render!(ax, system)                     # all objects of the system
render!(ax, beam)                       # a solved Beam (and its children)
save("setup.png", fig)                  # or display(fig) for an interactive window
```

## `render!` methods

| Target                        | Useful keywords (defaults) |
|-------------------------------|----------------------------|
| `System`                      | forwarded to each object |
| optical object                | `material` (e.g. `:refractive`, `:mechanics`), `edges`, `color`, `alpha`, `transparency` (override the material of the look) |
| `Beam` / `Ray`                | `color = :blue`, `linewidth = 1`, `flen = 1.0` (length drawn for an open-ended last ray), `show_pos = false`, `show_polarization = false` |
| beam group (source)           | `render_every = 5` (draw every n-th beam), plus beam keywords |
| `GaussianBeamlet`             | `color = :red`, `flen = 0.1`, `show_beams = false` (show chief/waist/divergence rays), `r_res`, `z_res`, `transparency = true` |
| `AstigmaticGaussianBeamlet`   | as above, plus `show_polarization`, `pol_*` |
| polarizer / filter            | `show_transmission_axis = true`, `axis_color`, `axis_linewidth` |

**Always set `flen`** for beams whose last ray does not end on a surface; the default 1 m for rays
usually dwarfs a millimeter-scale setup.

## Look and lighting

- `set_render_look(:modern)` (default: clear glass, metallic mirrors, neutral mechanics) or
  `set_render_look(:cad)` (saturated materials with feature edges); applies to later `render!` calls.
- Each component class has a material (`:refractive`, `:reflective`, `:coating`, `:polarizer`,
  `:detector`, `:mechanics`); `render!(ax, obj; material = :mechanics, edges = false)` overrides it.
- `studio_lighting!(ax)`: key, fill and rim light relative to the camera.

## Views

- Top view (x-y plane, beam along +y upwards): `Axis3(...; azimuth = -π/2, elevation = π/2)`
- Side view along the axis: `Axis3(...; azimuth = 0.0, elevation = 1e-3)`
- Use `aspect = :data` so optics are not distorted; or give explicit `aspect` and `limits`.
- Remove decorations for figures: `hidedecorations!(ax); hidespines!(ax)`.

`LScene` helpers:

- `hide_axis(ls)`, `set_orthographic(ls)`
- `get_view(ls)`: rotate the scene by hand, call it, paste the printed matrix into the script
- `set_view(ls, M)` or `set_view(ls, eye, lookat, up)`
- `look_at!(ls, target, offset; up = [0, 0, 1])`
- `arrow!(ax, pos, dir; scale)`, `render_lcs!(ax, obj; scale, show_labels)` (draw local frames)
- `view_cube!(ls)`: CAD-style cube in a corner, a click sets the standard view

## Live rendering and interactive window (GLMakie)

`render!` creates new plots on every call. For animations or parameter sweeps, render once and update:

```julia
hsys = live_render!(ax, system)   # geometry generated once
hbeam = live_render!(ax, beam)    # all ray segments in one plot
zrotate3d!(mirror, 1e-3)
solve_system!(system, beam)
update_render!(hsys); update_render!(hbeam)   # remove_render!(h) deletes the plots
```

`gui = live_view(system, beam); display(gui)` opens a window in which components are moved with
the mouse and keyboard (`kinematic_controls!`, `h` shows all controls), the system is re-solved
after each change and detector panels update live. `export_changes(gui)` prints the changed poses
as Julia code. Needs an interactive display; not for headless scripts.

Solves longer than `progress_delay` (kwarg, default 0.5 s) run in the background: the window stays
usable and a small progress window appears next to the source being traced or the detector whose
field is computed. Moving a component or pressing `Esc` cancels the solve.

## Detector data plots

```julia
x, z, I = intensity(det; n = 301, crop_factor = 5)
fig = Figure()
ax = Axis(fig[1, 1], aspect = 1, xlabel = "x [µm]", ylabel = "z [µm]")
heatmap!(ax, x * 1e6, z * 1e6, I)
save("psf.png", fig)

pts = spot_diagram(det)
scatter(first.(pts) * 1e6, last.(pts) * 1e6; markersize = 2)
```

See `templates/07_render_system.jl` for a complete headless example.
