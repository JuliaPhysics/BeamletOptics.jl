# Render gallery

`render_gallery.jl` renders every exported component of BeamletOptics.jl, one `ObjectGroup`, a
`MeshDummy` and the beam types (`Beam`, `CollimatedSource`, `GaussianBeamlet`,
`AstigmaticGaussianBeamlet`, each traced through a small lens system) into individual images. It
records render time and mesh statistics and serves as the baseline for work on the renderers.

## Running

From the repository root:

```
julia --project=docs gallery/render_gallery.jl
```

GLMakie renders offscreen. A full run takes roughly 1.5 to 2 minutes (most of it compilation). The
output folder `gallery/output/` is cleared at the start of every run and is not committed.

Before rendering, the script compares the catalogue with the component constructors exported in
`src/Exports.jl` and prints the names that are not covered (`missing: ...`). Add new components to
`CATALOGUE` in the script.

## Output

- `output/<category>/<name>_iso.png`: isometric view (from `+x`, `-y`, `+z`)
- `output/<category>/<name>_side.png`: side view along `x`, perpendicular to the optical axis (`+y`)
- `output/index.md`: one table per category
- `output/overview.png`: all isometric images in a grid (saved with CairoMakie, since GLMakie clips figures larger than the screen)

Both views are framed by the script: the automatic recentering of the `LScene` camera is switched off and the eye is placed such that the bounding sphere of the object fits into the field of view.

## Columns of `index.md`

| Column | Meaning |
|:--|:--|
| name | catalogue entry, usually the exported constructor |
| iso, side | the two images |
| time [ms] | duration of `render!` for the isometric view; a warm-up call before it compiles the methods, so this is the second call |
| plots | number of plots `render!` added to the scene |
| vertices, faces | total vertex and face count of all mesh and surface plots (including child plots); for surfaces the grid points and two triangles per grid cell |
| fallback | **yes** if any shape of the entry is rendered by the generic marching cubes method `render!(ax, ::AbstractSDF)` |
| error | first line of the error if building or rendering the entry failed; the run continues with the next entry |
