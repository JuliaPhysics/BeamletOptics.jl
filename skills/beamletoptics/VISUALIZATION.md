# Visualization (Makie extension)

Rendering lives in the `BeamletOpticsMakieExt` extension. It activates once a Makie backend is loaded:

```julia
using CairoMakie, BeamletOptics   # headless PNG/PDF/SVG; only Axis3
# or
using GLMakie, BeamletOptics      # interactive window; LScene + Axis3, camera helpers
```

Without a backend, `render!` throws a `MissingBackendError`.

## Minimal example

```julia
fig = Figure()
ax = Axis3(fig[1, 1], aspect = :data)   # or LScene(fig[1, 1]) with GLMakie
render!(ax, system)                     # all objects of the system
render!(ax, beam)                       # a solved Beam (and its children)
save("setup.png", fig)                  # CairoMakie; display(fig) with GLMakie
```

## `render!` methods

| Target                        | Useful keywords (defaults) |
|-------------------------------|----------------------------|
| `System`                      | forwarded to each object |
| optical object                | `color`, `transparency`, `alpha` |
| `Beam` / `Ray`                | `color = :blue`, `linewidth = 1`, `flen = 1.0` (length drawn for an open-ended last ray), `show_pos = false`, `show_polarization = false` |
| beam group (source)           | `render_every = 5` (draw every n-th beam), plus beam keywords |
| `GaussianBeamlet`             | `color = :red`, `flen = 0.1`, `show_beams = false` (show chief/waist/divergence rays), `r_res`, `z_res`, `transparency = true` |
| `AstigmaticGaussianBeamlet`   | as above, plus `show_polarization`, `pol_*` |
| polarizer / filter            | `show_transmission_axis = true`, `axis_color`, `axis_linewidth` |

**Always set `flen`** for beams whose last ray does not end on a surface; the default 1 m for rays
usually dwarfs a millimeter-scale setup.

## Views

- Top view (x-y plane, beam along +y upwards): `Axis3(...; azimuth = -π/2, elevation = π/2)`
- Side view along the axis: `Axis3(...; azimuth = 0.0, elevation = 1e-3)`
- Use `aspect = :data` so optics are not distorted; or give explicit `aspect` and `limits`.
- Remove decorations for figures: `hidedecorations!(ax); hidespines!(ax)`.

GLMakie `LScene` helpers:

- `hide_axis(ls)`, `set_orthographic(ls)`
- `get_view(ls)`: rotate the scene by hand, call it, paste the printed matrix into the script
- `set_view(ls, M)` or `set_view(ls, eye, lookat, up)`
- `look_at!(ls, target, offset; up = [0, 0, 1])`
- `arrow!(ax, pos, dir; scale)`, `render_lcs!(ax, obj; scale, show_labels)` (draw local frames)

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
