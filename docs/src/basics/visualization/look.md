# Look

Rendered objects get a material per component class and visible cemented interfaces of doublet
and triplet lenses. Two looks are available, which are selected via [`set_render_look`](@ref):

- `:modern` (default): a restrained palette of clear, slightly tinted glass with highlights,
  metallic mirrors and neutral mechanics. Only the glass (`:refractive`, `:coating` and
  `:interface`) gets faint silhouettes, i.e. its feature edges at 60 % of the opacity of the
  `:cad` look, such that lenses stay readable in front of a housing
- `:cad`: saturated materials with thin dark lines along the feature edges, like a CAD program

```julia
set_render_look(:cad)
```

| material      | components                                 | `:modern`            | `:cad`             |
|:--------------|:-------------------------------------------|:---------------------|:-------------------|
| `:refractive` | lenses, prisms, plates, windows            | clear glass, 0.3     | light blue, 0.5    |
| `:reflective` | mirrors, retroreflector                    | metallic silver      | silver             |
| `:coating`    | beamsplitter coatings                      | pale violet, 0.28    | magenta, 0.6       |
| `:polarizer`  | polarization filters                       | dark slate, 0.75     | dark teal, 0.8     |
| `:detector`   | detectors                                  | graphite             | dark blue          |
| `:mechanics`  | mechanics, dummies and other objects       | neutral grey         | mid grey           |
| `:interface`  | cemented interfaces of doublets, triplets  | pale amber, 0.15     | amber, 0.25        |

The numbers are the opacity `alpha` of transparent materials.

Besides the color and the opacity, each material sets the `transparency`, `diffuse`, `specular`
and `shininess` attributes of the mesh plot. The parts of composite objects are rendered with
their own class, e.g. the prisms and the coating of a [`CubeBeamsplitter`](@ref). The following
keyword arguments of `render!` change the look of an object:

- `material = nothing`: one of the materials above, overrides the component class for all parts
  of the object
- `color`, `alpha`, `transparency`, ...: override the corresponding attribute of the material,
  for all parts of the object. The cemented interfaces keep their amber look.
- `edges`: draws the feature edges, i.e. the edges where the faces of an object meet at an angle
  of more than 30°, and the boundary of open surfaces such as a `Detector`. By default, the `:cad`
  look draws them for all materials except `:mechanics`, whose detailed meshes, e.g. a housing from
  an STL file, would cover the optics, the `:modern` look for the glass only. The parts of
  composite objects contribute edges according to their own class. `edges = true` or `false`
  overrides the look, for all parts of the object. Shapes rendered via the marching cubes fallback
  have no edges. The opacity of the edges follows the opacity of the object, e.g. a nearly
  transparent housing (`color = (:gray70, 0.05)`) gets correspondingly faint edges, while the
  edges of glass stay visible.

```julia
render!(ax, lens)                       # glass of the active look
render!(ax, mirror; edges = true)       # with feature edges in the :modern look
render!(ax, lens; color = :red)         # red glass, same transparency
render!(ax, mount; material = :mechanics, edges = false)
```

## Lighting

The lighting of the scene is set via [`studio_lighting!`](@ref): an ambient light plus a key
light from the upper right front, a fill light from the left and a rim light from behind, all
relative to the camera. Backends with a single directional light (e.g. CairoMakie) get the
ambient and the key light only.

```julia
fig = Figure()
ax = LScene(fig[1, 1])
studio_lighting!(ax)
render!(ax, system)
```

`studio_lighting!(ax; preset = :none)` keeps the default lights of Makie and
`edges = true` or `false` overrides the edges of the look.
