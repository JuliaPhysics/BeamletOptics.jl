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

With GLMakie, `display(gui)` opens the window without SSAO and with up to 60 fps, since SSAO (e.g.
enabled globally via `GLMakie.activate!(ssao = true)`) multiplies the frame time with large meshes
such as a housing, and GLMakie's default of 30 fps makes rotating the view sluggish. Screen settings
passed to `display(gui; ...)` override these defaults. For static renders with SSAO, pass it to the
display of that figure only, e.g. `display(fig; ssao = true)`.

A [view cube](@ref "View cube") in the top right corner of the 3D view switches to the standard
views with a click, clicks on the cube never select or deselect a component. The live view starts
in the isometric view from the corner between `Top`, `Front` and `Right`, i.e. from `(1, -1, 1)`,
such that the labels of the cube read correctly. The "orthographic" toggle below the 3D view,
next to "auto trace" and "clip beams", switches between perspective and orthographic projection,
`orthographic = true` starts with the latter. Moving along the view direction does not change an
orthographic view, there `W`/`S` zoom like `U`/`O`:

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
direction, the green axis is its direction. If a marker covers small components, the "sources"
toggle below the 3D view or the key `1` hides all markers and shows them again, `show_sources =
false` starts with hidden markers. Pass `movable_sources = false` to omit the markers altogether.

### Extras and static context

Objects that are not part of any `system`, e.g. a housing from a CAD file, are passed as `extras`,
each optionally with the kwargs of its `render!` call:

```julia
housing = MeshDummy("housing.stl")
gui = live_view(system => beam; extras = [housing => (; transparency = true, color = (:gray, 0.3))])
```

Extras are never traced, so they cost nothing in the solves, and moving them does not solve the
systems. Otherwise they act like components: they are selected with a click, in the component menu
or in the object tree (under "Extras"), moved, hidden and exported. They do not block clicking on
the optics behind them: objects are picked by intersecting the camera ray with the optical
components first, which a `NonInteractableObject` never intersects. The card of a
`NonInteractableObject` (e.g. a `MeshDummy`) or an `IntersectableObject` has an "opacity" slider,
which makes it transparent without changing the optics; 0 % hides it. Below 50 % opacity, a click
in the 3D view passes through such an object, i.e. a click into the empty space inside a
transparent housing does not select the housing; select it in the object tree or the component
menu instead. Visible extras count towards the size of the scene, which sets the size of the
source markers and of new clip planes and the view of "fit all" (`g` without a selection).

Plots added directly to `gui.ax` via `render!` are only drawn, i.e. they are neither selectable nor
hideable:

```julia
render!(gui.ax, housing_mesh; transparency = true, color = (:gray, 0.3))
```

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

The subtitle of each panel shows its metrics for alignment, the centroid is marked by a red cross:

- spot diagram: the number of hits `N`, the centroid `c`, the RMS radius `sqrt(mean(|p - c|²))`
  and the geometric radius, i.e. the largest distance from the centroid
- intensity: the power `P`, the peak intensity, the centroid `c` and the 1/e² radii `w` along x
  and z, i.e. twice the standard deviation of the intensity along each axis

The following panel options are not passed to [`intensity`](@ref):

| Option                  | Effect                                                               |
|:------------------------|:---------------------------------------------------------------------|
| `colorscale = :log`     | shows `log10` of the intensity, with a floor of 1e-4 times the maximum |
| `colorrange = (lo, hi)` | fixed color range of the intensity (in `log10` units for `:log`)     |
| `history = true`        | adds an axis with the power (or `N`) and the centroid over the last 300 full solves |
| `profiles = true`       | adds an axis with the intensity along x and z through the centroid   |

```julia
gui = live_view(system, beam; detectors = [pd => (:intensity, (; colorscale = :log, history = true, profiles = true))])
```

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

`on_change = (gui, obj) -> ...` is called after every full solve, with the moved object or
`nothing` (initial solve, or after a slider change). It is not called after the preview solves of
beam groups while moving, see [Manual tracing](@ref), but once the full solve follows. Use it to
record derived quantities, e.g. the optical power on a detector over time, which an own panel
plots, see [Own panels, controls and tools](@ref):

```julia
power = Point2f[]

function record_power!(gui, obj)
    push!(power, Point2f(length(power) + 1, 1e3 * optical_power(pd)))
    return nothing
end

gui = live_view(system, beam; on_change = record_power!)
```

Errors raised inside `on_change` are logged once and do not interrupt the interaction. See the
[Interactive Michelson interferometer](@ref) example for a full callback that tracks the optical
power over time.

### Own panels, controls and tools

Own parts are added to a `gui = live_view(...)` by three functions, which work with both layouts;
the layout decides where the parts go:

| function | `layout = :compact` | `layout = :app` |
|:--|:--|:--|
| [`add_panel!`](@ref) | below the detector panels, right of the 3D view | a tab of the analysis dock |
| [`add_controls!`](@ref) | a row above the status row | a section of the left sidebar, below "Parameters" |
| [`add_tool!`](@ref) | a button (or toggle) in the row below the status line | an icon button (or toggle) in the toolbar, before "Help" |

`add_panel!(f, gui, title)` calls `f(layout)` with the `GridLayout` of the new panel, into which it
builds e.g. an `Axis` with plots. The function returned by `f` is called with the `gui` after each
full solve, like `on_change`, and once right away. In the app layout, only the panel of the shown
tab is updated, hidden panels when their tab is opened, i.e. hidden panels cost nothing. The new tab
stays in the background unless `select = true`:

```julia
add_panel!(gui, "Optical power") do layout
    ax = Axis(layout[1, 1]; xlabel = "Update", ylabel = "P [mW]")
    pts = Observable(copy(power))
    lines!(ax, pts)
    return gui -> (pts[] = copy(power); autolimits!(ax))
end
```

Data that must be recorded after every solve, also while the panel is hidden, is recorded by
`on_change`, as above, and only plotted by the panel. In the compact layout, the panels are placed
in a grid in `gui.fig[1, 2]`; an axis at `gui.fig[1, 2][2, 1]`, below a single detector panel,
works as well, but only `add_panel!` works with both layouts.

`add_controls!(f, gui, title)` calls `f(layout)` to build widgets, e.g. buttons, menus and
textboxes. While a textbox of the controls is focused or a menu is open, the keys of the 3D view and
the camera are ignored, as for the boxes of the live view. A change of the optics from such a
widget calls [`retrace!`](@ref), which solves the systems again like after a slider of the
`sliders` kwarg (or marks the beams as outdated with `auto_trace = false`); the change itself goes
into its first argument, since objects must not change while a solve runs in the background:

```julia
add_controls!(gui, "Mirror") do layout
    box = Textbox(layout[1, 1]; placeholder = "tilt [mrad]", width = 100)
    on(box.stored_string) do s
        θ = tryparse(Float64, s)
        isnothing(θ) || retrace!(() -> zrotate3d!(m1, 1e-3 * θ), gui)
    end
end
```

`add_tool!(f, gui, name)` adds a button that calls `f(gui)`, with `toggle = true` a toggle that
calls `f(gui, active)`. In the app layout, it is an icon of the toolbar (`icon`, e.g. `:measure`,
`:export` or `:chart`, an unknown name lists the valid ones) with a `tooltip`. `key` binds a key of
the 3D view to the tool; keys that the live view, the kinematic controls or the camera use (all
letters) throw an `ArgumentError` naming the binding, e.g. the digits `2`–`9` or the function keys
are free:

```julia
add_tool!(gui, "Center m1"; icon = :fit, key = Keyboard._2) do gui
    retrace!(() -> translate_to3d!(m1, [0, 0.1, 0]), gui)
end
# Shows or hides a reference line along the optical axis
ref = lines!(gui.ax, [Point3f(0, 0, 0), Point3f(0, 0.3, 0)]; color = :gray, visible = false)
add_tool!(gui, "Optical axis"; toggle = true, key = Keyboard._3) do gui, active
    ref.visible[] = active
end
```

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

Beam groups, e.g. a source with thousands of rays, are rendered with `render_every = 5` by
default, i.e. only every fifth beam is drawn. While a component is moved, such groups are only
solved for the rendered beams (preview tracing), the other beams are reset and do not hit the
detectors; the titles of the detector panels end with "(preview)". Once the movement pauses for
`idle_delay`, the full group is solved. The `trace_budget` applies to the preview solve while
moving, such that large groups stay interactive. `preview = false` always solves the full groups:

```julia
gui = live_view(system, source; beam_kwargs = Dict(source => (; render_every = 50)), preview = false)
```

Solves longer than `progress_delay` run in the background with a progress window next to the
source or detector, see "Long solves" in the docstring of [`live_view`](@ref).

### Component card and component menu

Selecting a component, source or clip plane opens a card next to its bounding box in the 3D view,
connected to it by a line. The card follows the camera and the component and stays inside the
view, off the view cube and the other cards. It has the colors of `theme = :light` (default) or
`:dark`, like the app layout, also in the compact layout:

- The head shows the icon of the kind of the component, its label and "hide", which hides it, e.g. a mirror in front of
  the component of interest, and clears the selection. A hidden component can not be selected in
  the 3D view, but stays in the systems, i.e. it is still traced. Selected in the menu, its card
  reads "show", which shows it again. For a clip plane, the head holds "flip" and "remove".
- `x`, `y`, `z` [mm] show the position of the component. Typing a value and pressing `Enter`
  moves the component to this absolute coordinate. The boxes `rx`, `ry` and `rv` [mrad] rotate it
  by the typed angle about the red, green and blue axis of the controls, like the arrow keys in the
  rotate mode, e.g. `rv = 1` equals one key step with a step of 1 mrad. Each input is a step of the
  undo history, the constraints of the component apply. While a box is focused, the keys of the 3D
  view are ignored.
- `step` sets the keyboard step, see below.
- The chevron at the right end of the head collapses the card to its head and expands it again.
- The pin keeps the card with its component when the selection changes, e.g. to watch or type the
  poses of several components; the widgets of a pinned card act on its component. Clicking the
  pin again closes it.
- Dragging the head (icon, label or the free room around them) moves the card to another place in
  the view, e.g. to line up several pinned cards at an edge. It stays there, connected to its
  component by the line, when the camera or the component moves and when it is pinned; the card of
  the selection also keeps its place for the next selected component. The place is kept relative
  to the nearest corner of the view, so a card at an edge stays there when the window is resized.
  A double click on the head places the card next to its component again, as does unpinning it.
- Below the pose, rows of the component type, refreshed after each solve and move:

  | Component | Rows |
  |:----------|:-----|
  | optical components | `beam`: the rays hitting it in the last solve and the angle of incidence of the first one (e.g. the center ray of a ring source) with the range of all, or "not hit" |
  | lenses, prisms | `n` at the wavelength of the hitting beam, `d` the center thickness (`Lens`) |
  | beamsplitters | `split`: R and T (power) of the coating |
  | polarizers | `axis`: the transmission axis about the optical axis, from the horizontal |
  | detectors | `signal`: the power (intensity panel) or the number of rays (spot panel); `panel`: a button that cycles the mode of its detector panel (`auto`, `spot`, `intensity`) and a toggle of the logarithmic color scale, "no panel" for a detector without a panel (see "Detector panels") |
  | ray sources | `λ` and the diameter or NA; sources whose rays can be regenerated (`CollimatedSource`, `PointSource` and their uniform variants, see [`set_num_rays!`](@ref)) add the slider "rays" for their number of rays, which solves again |
  | Gaussian beamlets | `λ`, the waist `w0` and the Rayleigh range `zR` |

The rows and the buttons in the head are declared per type by multiple dispatch, see
[`card_rows`](@ref) and [`card_actions`](@ref): each row is a [`CardRow`](@ref) of texts and
[`CardWidget`](@ref)s, i.e. any `Makie` widget with a function for the value it shows and one for
its input. The card places the widgets, hides them and keeps their clicks and keys from the 3D
view. An own component adds its rows to the pose rows, e.g.

```julia
BeamletOptics.card_rows(l::MyLens) = (pose_card_rows(l)...,
    CardRow("f", CardWidget(Label; value = (gui, l) -> "$(round(1e3 * focal_length(l); digits = 2)) mm")),
    CardRow("n", CardWidget(Slider; range = 1.4:0.01:1.9, solve = true,
        value = (gui, l) -> refractive_index(l), on = (gui, l, n) -> set_index!(l, n))))
```

where `solve = true` solves the systems again after an input, like a move.

With `live_view(...; layout = :app)`, the cards are docked in the "Properties" sidebar instead of
floating next to the components: the card of the selection at the top, with the same rows and
actions, and below it the pinned cards, one below the other, each with its own head (icon, label,
actions, pin and chevron). The pin of the selection pins a card, the pin of a pinned card unpins
it. The sidebar does not scroll: if the cards do not fit, the older pinned cards collapse to their
heads (the one pinned or expanded last stays open) and the property list of the selection is
shortened; a collapsed card is only expanded again by its chevron. The float button in the head of
a pinned card moves it out of the sidebar into the 3D view, where it floats next to its component
as in the compact layout; the dock button in its head moves it back. Only the docked cards take
room in the sidebar. The floating cards and the docked cards are built by the same code from the
same declarations.

The row below the status line holds the component menu and "show all". The menu lists all movable
components and sources by their `labels` (or type), the objects of a group indented after the
group. Selecting an entry selects the component like a click in the 3D view, a click in the 3D view
shows the selected component in the menu. The menu can be searched by typing while it is open.
Clip planes are not listed. "show all" shows all hidden components.

### Beam inspection and measuring

A click on a beam, within 6 px of a rendered segment, marks the point on the beam and shows in the
status line its position [mm], the direction of the beam, the geometric path length and the
optical path length (Σ n·L) from the source [mm] and, for Gaussian beamlets, the radius `w` and
the radius of curvature `R` at this point (see `BeamletOptics.gauss_parameters`). Components take
precedence over beams, i.e. a click on a component still selects it. `esc` or a click elsewhere
removes the marker.

The "measure" toggle in the row of the component menu switches measuring on: two clicks on
components or beams show the distance between the positions of the components or the points of
the beams [mm], its components Δ, and the angle between the optical axes (local y-axes) of two
components, with a dashed line between the points. A third click starts a new measurement,
switching the toggle off clears it.

Both results are also shown on a card at their points in the 3D view, in the style of the
component cards: "Beam" with the rows `at`, `dir`, `path`, `OPL` (and `w`, `R`) and
"Measurement" with `from`, `to`, `d`, `Δ` (and `angle`). The card is replaced by the next
inspection or measurement and removed with it, unless its pin is clicked: a pinned card stays with
its marker or its line until it is unpinned, so several points and distances can be compared. In
the app layout, a pinned result is docked below the inspector like a pinned component.

### Camera tools

| Input       | Action                                                                     |
|:------------|:---------------------------------------------------------------------------|
| `g`         | Zoom to the selected component, or to all systems, keeping the view direction |
| "home"      | Restore the view when the window was shown                                 |
| "views"     | Set one of the saved views                                                 |
| "save view" | Save the current view as `"view n"` and print it as code for `views`       |

Saved views can be passed to the next session via `views`:

```julia
gui = live_view(system, beam;
    views = ["top" => ([0.0, 0.05, 0.5], [0.0, 0.05, 0.0], [0.0, 1.0, 0.0])])
```

### Exporting the changes

The "Export" button next to the status line prints the changed poses as Julia code and copies it
to the clipboard, such that an alignment found interactively can be pasted into the script that
builds the system. [`export_changes`](@ref) returns the same code:

```julia
gui = live_view(system, beam; labels = Dict(m1 => "m1", lens => "lens"))
# move the components, then
code = export_changes(gui)
```

```julia
# Changed poses of the live view, apply to the objects in their initial poses.
# Each rotation is about the position of the object, groups are moved before their objects.

# m1 (Mirror)
rotate3d!(m1, [0.0, 0.0, 1.0], 0.0005)
translate_to3d!(m1, [0.0, 0.1000012, 0.0])
```

Each moved object gets a `rotate3d!` about its own position (only if it was rotated) and a
`translate_to3d!` to its absolute position [m], relative to its pose when the window was opened.
Labels that are valid variable names are used as names, other objects are called `obj1`, `obj2`,
… by their position in the component menu. Clip planes are not exported.

### Controls

The 3D view uses the controls of [`kinematic_controls!`](@ref), see
[Interactive kinematics](@ref). In addition, the key `t` solves the systems immediately, see
[Manual tracing](@ref), `p`, `Delete`, `c` and `Shift+c` control the clip planes, see
[Clip planes](@ref), `1` shows or hides the source markers, see
[Movable sources in the live view](@ref), and `g` zooms to the selection, see [Camera tools](@ref). The keyboard step can be typed into the box `step` of the component card, e.g.
`250 nm` or `50 µrad`, where the unit selects the move or rotate mode, see
[Component card and component menu](@ref). The status line shows the
pose of the moved component and its change since the window was opened. Names for the status line
and the detector panels are passed via `labels`:

```julia
gui = live_view(system, beam; labels = Dict(m1 => "Mirror 1", pd => "Photodiode"),
    constraints = Dict(m1 => (; move = ())))
```

A complete example, including a custom `on_change` callback, can be found in the
[Interactive Michelson interferometer](@ref) example.
