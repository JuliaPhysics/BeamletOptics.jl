# Visualization

As mentioned in other sections of this documentation, the [Makie](https://docs.makie.org) backend can be used in order to generate 2D/3D renderings of optical systems and results generated with this package. Refer to the extensive `Makie` documentation and the **Examples** and the **Tutorials** sections of this package for a variety of showcases on how to visualize your simulation.

## Rendering elements

The main function provided for visualization purposes is the [`render!`](@ref) function. 

```@docs; canonical=false
render!(::Any, ::BeamletOptics._RenderTypes)
```

If a suitable backend is loaded, additional dispatched `render!` functions will become available. For instance, this allows the plotting of a [`GaussianBeamlet`](@ref).

## Loading the extension

Refer to the following snippet for an example on how the extension loading behaves. When only BMO is loaded, the `render!` function becomes available but will throw an [`BeamletOptics.MissingBackendError`](@ref) when trying to plot something.

```julia
julia> using BeamletOptics

julia> methods(render!)
# 1 method for generic function "render!" from BeamletOptics:
 [1] render!(::Any, ::Union{BeamletOptics.AbstractSystem, BeamletOptics.AbstractBeam, BeamletOptics.AbstractObject, BeamletOptics.AbstractObjectGroup, BeamletOptics.AbstractRay, BeamletOptics.AbstractShape}, kwargs...)
     @ C:\Users\anon\.julia\dev\BeamletOptics\src\Render.jl:56

julia> axis = nothing;

julia> mirror = RoundPlanoMirror(25e-3, 5e-3);

julia> render!(axis, mirror)
ERROR: It appears no suitable Makie backend is loaded in this session.
Stacktrace:
 [1] render!(::Nothing, ::Mirror{Float64, BeamletOptics.PlanoSurfaceSDF{Float64}})
   @ BeamletOptics c:\Users\anon\.julia\dev\BeamletOptics\src\Render.jl:46
 [2] top-level scope
   @ REPL[5]:1
```

Once a backend has been loaded, additional dispatched versions of `render!` become available.

```julia
julia> using GLMakie

julia> methods(render!)
# 21 methods for generic function "render!" from BeamletOptics:
  [1] render!(ax::Union{Axis3, LScene}, s::BeamletOptics.UnionSDF; kwargs...)
     @ BeamletOpticsMakieExt C:\Users\anon\.julia\dev\BeamletOptics\ext\RenderSDF.jl:32
  [2] render!(axis::Union{Axis3, LScene}, css::BeamletOptics.ConcaveSphericalSurfaceSDF; color, kwargs...)
     @ BeamletOpticsMakieExt C:\Users\anon\.julia\dev\BeamletOptics\ext\RenderLenses.jl:1
  [3] render!(axis::Union{Axis3, LScene}, css::BeamletOptics.ConvexSphericalSurfaceSDF; color, kwargs...)
     @ BeamletOpticsMakieExt C:\Users\anon\.julia\dev\BeamletOptics\ext\RenderLenses.jl:31
  [4] render!(axis::Union{Axis3, LScene}, acyl::BeamletOptics.AbstractAcylindricalSurfaceSDF; color, kwargs...)
     @ BeamletOpticsMakieExt C:\Users\anon\.julia\dev\BeamletOptics\ext\RenderCylinderLenses.jl:1
  [5] etc...
```

## Look

Rendered objects get a material per component class and visible cemented interfaces of doublet
and triplet lenses. Two looks are available, which are selected via [`set_render_look`](@ref):

- `:modern` (default): a restrained palette of clear, slightly tinted glass with highlights,
  metallic mirrors and neutral mechanics, without edge lines
- `:cad`: saturated materials with thin dark lines along the feature edges, like a CAD program

```julia
set_render_look(:cad)
```

| material      | components                                 | `:modern`            | `:cad`             |
|:--------------|:-------------------------------------------|:---------------------|:-------------------|
| `:refractive` | lenses, prisms, plates, windows            | clear glass, 0.38    | light blue, 0.5    |
| `:reflective` | mirrors, retroreflector                    | metallic silver      | silver             |
| `:coating`    | beamsplitter coatings                      | light violet, 0.4    | magenta, 0.6       |
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
- `edges`: draws the feature edges (by default in the `:cad` look for all materials except
  `:mechanics`, whose detailed meshes, e.g. a housing from an STL file, would cover the optics),
  i.e. the edges where the faces of an object meet at
  an angle of more than 30°, and the boundary of open surfaces such as a `Detector`. `edges =
  false` switches them off. Shapes rendered via the marching cubes fallback have no edges. The
  opacity of the edges follows the opacity of the object, e.g. a nearly transparent housing
  (`color = (:gray70, 0.05)`) gets correspondingly faint edges, while the edges of glass stay fully
  visible.

```julia
render!(ax, lens)                       # glass of the active look
render!(ax, lens; edges = true)         # with feature edges
render!(ax, lens; color = :red)         # red glass, same transparency
render!(ax, mount; material = :mechanics, edges = false)
```

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

The [`live_view`](@ref) applies the rig, `lighting = :none` keeps the default lights of Makie and
`edges = true` or `false` overrides the edges of the look.

## Camera and scene helpers

Alongside `render!`, a small set of `LScene`-specific helpers is provided for framing and
annotating a 3D scene once a backend is loaded: [`get_view`](@ref), [`set_view`](@ref),
[`set_orthographic`](@ref), [`hide_axis`](@ref), [`look_at!`](@ref), [`arrow!`](@ref) and
[`render_lcs!`](@ref). Like `render!`, each throws a
[`BeamletOptics.MissingBackendError`](@ref) if called before a suitable backend has been
loaded. Refer to the **Reference** page for their full docstrings.

The `get_view`/`set_view(ls, matrix)` pair is meant for interactive use: rotate the scene
by hand, call `get_view(ax)`, and paste the printed matrix back into the script as a
literal passed to `set_view`. This is the pattern used throughout this package's own
tutorials to freeze a camera position found interactively. The `set_view(ls, eye, lookat,
up)` and `look_at!` forms are the reproducible alternative, useful when the viewpoint
should be derived from the scene's own geometry instead of copy-pasted.

## View cube

[`view_cube!`](@ref) adds a small CAD-style cube to a corner of an `LScene`. The cube rotates
with the camera and shows the current orientation of the scene. A left click on a face, an edge or
a corner of the cube moves the camera to the corresponding standard view in a short animation,
while the point the camera looks at and its distance are kept. The region under the cursor is
highlighted.

```julia
using GLMakie, BeamletOptics

fig = Figure()
ax = LScene(fig[1, 1])
render!(ax, system)
cube = view_cube!(ax; size = 110, corner = :top_right)
display(fig)
```

The camera is placed on the side of the clicked face and looks at the system from there:

| Face     | Camera at | Up   |
|:---------|:----------|:-----|
| `Top`    | `+z`      | `+y` |
| `Bottom` | `-z`      | `+y` |
| `Front`  | `-y`      | `+z` |
| `Back`   | `+y`      | `+z` |
| `Right`  | `+x`      | `+z` |
| `Left`   | `-x`      | `+z` |

The 12 edges and 8 corners give the diagonal views between the adjacent faces, e.g. the corner
between `Top`, `Front` and `Right` looks from `(1, -1, 1)`, with `+z` as the up direction. The
transition takes `duration = 0.3` s, `duration = 0` switches the view instantly. The cube is not
affected by zooming or by clip planes. Call `close(cube)` to remove it. The
[live view](@ref "Interactive live view") shows a view cube by default, which is disabled via
`live_view(system, beam; view_cube = false)`.

## Live rendering

`render!` generates new plots on every call. For animations and interactive applications, where
components are moved and the system is re-solved many times per second, use
[`live_render!`](@ref) instead. It returns a handle that re-synchronizes the existing plots with
the current simulation state via [`update_render!`](@ref):

```julia
using GLMakie, BeamletOptics

fig = Figure()
ax = LScene(fig[1, 1])

hsys = live_render!(ax, system)   # geometry is generated once
hbeam = live_render!(ax, beam)    # all ray segments bundled into a single plot

zrotate3d!(mirror, 1e-3)
solve_system!(system, beam)
update_render!(hsys)              # moved components: only their model transformation changes
update_render!(hbeam)             # beam path: point buffer is replaced in place
```

- **Objects and systems**: the geometry is rendered once with `render!`. Kinematic changes
  (`translate3d!`, `rotate3d!`, ...) are then applied as a rigid model transformation of the
  existing plots, which costs microseconds regardless of the mesh resolution.
- **Rays, beams and beam groups**: all segments are drawn by a single `linesegments` plot, so the
  number of plots does not grow with the number of rays. The number of segments may change
  between updates, e.g. when a component is moved out of the beam path.
- **Gaussian beamlets**: the 1/e² envelope of all segments is merged into a single mesh. This
  applies to `GaussianBeamlet`s, `AstigmaticGaussianBeamlet`s and groups of astigmatic beamlets,
  of which every `render_every`-th beamlet is rendered.

Use [`remove_render!`](@ref) to delete the plots of a handle.

## Interactive kinematics

With [`kinematic_controls!`](@ref), the components of a live-rendered system can be grabbed and
moved with the mouse. The `on_change` callback is invoked (at most once per frame) after every
change and is the place to re-solve the system and update dependent plots:

```julia
ctrl = kinematic_controls!(ax, hsys; on_change = obj -> begin
    empty!(detector)
    solve_system!(system, beam)
    update_render!(hbeam)
end)
```

A click selects a component, a drag on the selected component moves or rotates it, and every
other drag rotates the camera as usual, so that rotating the camera never selects or moves a
component by accident. The keyboard controls apply to the selected component and depend on the
mode, which is switched with `m`:

| Input                                   | Move mode                     | Rotate mode                   |
|:----------------------------------------|:------------------------------|:------------------------------|
| Left-click on a component               | Select it                     | Select it                     |
| Left-drag on the selected component     | Move in the horizontal plane  | Rotate around the blue axis   |
| Left-drag elsewhere                     | Rotate the camera             | Rotate the camera             |
| `↑` / `↓`                               | Move along the green arrow    | Rotate around the red ring    |
| `→` / `←`                               | Move along the red arrow      | Rotate around the blue ring   |
| `Page Up` / `Page Down`                 | Move along the blue arrow     | Rotate around the green ring  |

| Input                        | Action                                                          |
|:-----------------------------|:----------------------------------------------------------------|
| `m`                          | Switch between the move and the rotate mode                     |
| Shift (held)                 | Ten times the step size                                         |
| `+` / `-`                    | Increase / decrease the step size along the 1-2-5 sequence      |
| `Backspace`                  | Reset the selected component to its initial pose                |
| `Ctrl`/`Cmd` + `Z`           | Undo the last change                                            |
| `Ctrl`/`Cmd` + `Y`           | Redo, also `Ctrl`/`Cmd` + `Shift` + `Z`                         |
| `Esc`                        | Select the enclosing group, or deselect at the top level        |
| Left-click on empty space    | Deselect                                                        |
| `v`                          | Switch the spectator mode on or off                             |
| `h`                          | Show or hide an overlay of all controls                         |

In the spectator mode, the selection is cleared and all mouse and keyboard input goes to the
camera, such that the system can be viewed without moving a component by accident. Components
whose kinematic trait is `Static` can not be selected.

The selected component is marked by a box and three axes above it: its local y-axis (green), its
local x-axis (red) and the vertical rotation axis (blue). In the move mode the axes are shown as
arrows, in the rotate mode as rings. The first key of each pair moves the component in the
direction of the arrow, or rotates it in the direction of the ring. The current mode and step
size are shown in the hint line at the top of the 3D view.

Clicking a component inside an `ObjectGroup` selects the outermost group first. Clicking the same
component again descends one level into the hierarchy (a subgroup, then the individual object),
so that the group can still be moved as a whole, or a single part can be moved on its own. `Esc`
goes back up one level.

A drag grabs the point under the cursor, which stays under the cursor during the drag. Each drag,
reset and series of steps with the same key (less than 1 s apart) is one entry of the undo history,
which undoes up to 100 changes.

The `constraints` lock axes of individual components, e.g. a mirror in a kinematic mount that can
only be tilted. The axes are named after the gizmo: `:x` (red), `:y` (green) and `:v` (blue, the
`rotation_axis`). Missing fields allow all axes, locked axes are shown faded:

```julia
ctrl = kinematic_controls!(ax, hsys; constraints = Dict(mirror => (; move = (), rotate = (:x, :v))))
```

Call `close(ctrl)` to remove the controls.

## Interactive live view

[`live_view`](@ref) combines [`live_render!`](@ref), [`kinematic_controls!`](@ref), detector
panels and optional sliders into a single ready-to-use window. It is the fastest way to explore
the sensitivity of a system in the REPL: grab a mirror, watch the beam path and the detector
panels update live.

```julia
using GLMakie, BeamletOptics

gui = live_view(system, beam)
display(gui)
```

A [view cube](@ref "View cube") in the top right corner of the 3D view switches to the standard
views with a click, clicks on the cube never select or deselect a component. The live view starts
in the isometric view from the corner between `Top`, `Front` and `Right`, i.e. from `(1, -1, 1)`,
such that the labels of the cube read correctly. The "orthographic" toggle below the 3D view,
next to "auto trace" and "clip beams", switches between perspective and orthographic projection,
`orthographic = true` starts with the latter:

```julia
gui = live_view(system, beam; orthographic = true)
```

More than one `system => beam` pair can be shown in the same 3D view, e.g. the transmitter and
receiver path of a lidar, which are solved with different sources:

```julia
gui = live_view(system_tx => beam_tx, system_rx => source_rx)
```

### Movable sources in the live view

Each source, i.e. the beam or beam group of each `system => beam` pair, is shown with an orange
marker at its position, which points along its direction. A source is selected and moved via its
marker like any component, after which the systems are solved again. For a beam, which only has a
direction, the green axis is its direction. Pass `movable_sources = false` to hide the markers.

### Static context

Additional context that is not part of any `system`, e.g. a housing or an optical table, can be
added directly to `gui.ax` via `render!`:

```julia
render!(gui.ax, housing_mesh; transparency = true, color = (:gray, 0.3))
```

Such geometry is not selectable and does not block clicking on the optics behind it: objects are
picked by intersecting the camera ray with the movable objects of the system, not with everything
drawn in the scene, so a housing mesh in front of a component never gets in the way.

### Clip planes

Clip planes cut the 3D view open, e.g. to look into a housing or behind a component. Only the side
of a plane its normal points to stays visible. `p` adds a plane through the selected component, or
through the point the camera looks at if nothing is selected, with its normal along the view
direction. A plane is selected via the purple handle at its center and moved and rotated like a
component, its normal is the green axis. Moving a plane does not solve the systems.

| Input      | Action                                                  |
|:-----------|:--------------------------------------------------------|
| `p`        | Add a clip plane and select it                          |
| `Delete`   | Remove the selected clip plane                          |
| `c`        | Switch clipping on or off (all planes)                  |
| `Shift+c`  | Flip the selected clip plane, i.e. show the other side  |

Planes can also be given at construction as `point => normal`. The beams are not clipped unless
`clip_beams = true` or the "clip beams" toggle below the 3D view is switched on, the markers of
the sources and planes are never clipped:

```julia
gui = live_view(system, beam; clip_planes = [[0, 0.1, 0] => [0, 1, 0]], clip_beams = true)
```

Makie supports at most 8 clip planes. The selection box of a partly clipped component only covers
its visible part.

### Detector panels

By default (`detectors = :auto`), one panel is shown for every `Detector` of every system,
deduplicated by identity. Each panel shows the spot diagram (`:spot`) for ray-based hits or the
intensity (`:intensity`) for Gaussian beamlet hits, chosen automatically (`:auto`). Pass a vector
to select detectors and modes explicitly, or `[]` to disable the panels:

```julia
gui = live_view(system, beam; detectors = [pd1, pd2 => :spot, pd3 => (:intensity, (; n = 200))])
```

The `kwargs` of the `pd => (mode, kwargs)` form are passed to [`intensity`](@ref); most useful is
a fixed extent via `x_min`, `x_max`, `z_min` and `z_max` (in meters, like the rest of this
package), instead of the automatic crop around the beam.

### Sliders

`sliders` adds custom parameters below the 3D view. Each entry is `"label" => (range, callback)`
(or `(range, callback, startvalue)`); `callback` is called with the current slider value and is
expected to move objects or otherwise change the system:

```julia
gui = live_view(system, beam;
    sliders = ["focus [mm]" => (-1:0.01:4, v -> set_focus(cl, v * 1e-3))])
```

Since this package uses SI units throughout, a slider that is labeled and ranged in millimeters
for convenience must convert its value before applying it, as in `v * 1e-3` above.

### Custom updates

`on_change = (gui, obj) -> ...` is called after every solve, with the moved object or `nothing`
(initial solve, or after a slider change). Use it to plot additional derived quantities into
`gui.fig`. Since `on_change` already runs for the initial solve inside `live_view`, the callback
should only update an `Observable`; the axis and the plot are created once afterwards:

```julia
power = Observable(Point2f[])

function record_power!(gui, obj)
    P = optical_power(pd)
    push!(power[], Point2f(length(power[]) + 1, 1e3 * P))
    notify(power)
    return nothing
end

gui = live_view(system, beam; on_change = record_power!)
# Below a single detector panel, the panels are placed in a grid in gui.fig[1, 2]
power_ax = Axis(gui.fig[1, 2][2, 1]; xlabel = "Update", ylabel = "P [mW]")
lines!(power_ax, power)
on(_ -> autolimits!(power_ax), power)
```

Errors raised inside `on_change` are logged once and do not interrupt the interaction. See the
[Interactive Michelson interferometer](@ref) example for a full callback that tracks the optical
power over time.

### Manual tracing

Solving a large system on every mouse-drag event can be too slow for smooth interaction. With
`auto_trace = false`, `live_view` still updates the 3D view and the sliders immediately, but only
solves the systems (and updates the beams and detector panels) on request: the `Trace (t)` button
below the 3D view, the key `t`, or switching the "auto trace" toggle back on (which solves once if
the state is outdated). While outdated, the beam plots are dimmed and the status line shows a
hint. The initial solve always runs, regardless of `auto_trace`.

With `auto_trace = true`, `live_view` adapts to slow systems as well: if solving takes longer than
`trace_budget` (30 ms by default), the components still follow the mouse immediately, while the
systems are solved once the movement pauses for `idle_delay` (0.2 s). Likewise, detector panels
that take longer than `trace_budget` show a coarse preview while moving, which is refined once the
movement pauses.

### Controls

The 3D view uses the controls of [`kinematic_controls!`](@ref), see
[Interactive kinematics](@ref). In addition, the key `t` solves the systems immediately, see
[Manual tracing](@ref), and `p`, `Delete`, `c` and `Shift+c` control the clip planes, see
[Clip planes](@ref). The keyboard step can be typed into the textbox below the 3D view, e.g.
`250 nm` or `50 µrad`, where the unit selects the move or rotate mode. The status line shows the
pose of the moved component and its change since the window was opened. Names for the status line
and the detector panels are passed via `labels`:

```julia
gui = live_view(system, beam; labels = Dict(m1 => "Mirror 1", pd => "Photodiode"),
    constraints = Dict(m1 => (; move = ())))
```

A complete example, including a custom `on_change` callback, can be found in the
[Interactive Michelson interferometer](@ref) example.
