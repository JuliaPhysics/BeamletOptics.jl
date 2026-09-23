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
 [1] render!(::Nothing, ::RoundPlanoMirror{Float64})
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
- **Gaussian beamlets**: the 1/e² envelope of all segments is merged into a single mesh.
  `AstigmaticGaussianBeamlet`s are not supported yet.

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
| `Esc`                        | Select the enclosing group, or deselect at the top level        |
| Left-click on empty space    | Deselect                                                        |
| `h`                          | Show or hide an overlay of all controls                         |

The selected component is marked by a box and three axes above it: its local y-axis (green), its
local x-axis (red) and the vertical rotation axis (blue). In the move mode the axes are shown as
arrows, in the rotate mode as rings. The first key of each pair moves the component in the
direction of the arrow, or rotates it in the direction of the ring. The current mode and step
size are shown in the hint line at the top of the 3D view.

Clicking a component inside an `ObjectGroup` selects the outermost group first. Clicking the same
component again descends one level into the hierarchy (a subgroup, then the individual object),
so that the group can still be moved as a whole, or a single part can be moved on its own. `Esc`
goes back up one level. Call `close(ctrl)` to remove the controls.

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

More than one `system => beam` pair can be shown in the same 3D view, e.g. the transmitter and
receiver path of a lidar, which are solved with different sources:

```julia
gui = live_view(system_tx => beam_tx, system_rx => source_rx)
```

### Static context

Additional context that is not part of any `system`, e.g. a housing or an optical table, can be
added directly to `gui.ax` via `render!`:

```julia
render!(gui.ax, housing_mesh; transparency = true, color = (:gray, 0.3))
```

Such geometry is not selectable and does not block clicking on the optics behind it: objects are
picked by intersecting the camera ray with the movable objects of the system, not with everything
drawn in the scene, so a housing mesh in front of a component never gets in the way.

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

### Controls

The 3D view uses the controls of [`kinematic_controls!`](@ref), see
[Interactive kinematics](@ref). In addition, the key `t` solves the systems immediately, see
[Manual tracing](@ref).

A complete example, including a custom `on_change` callback, can be found in the
[Interactive Michelson interferometer](@ref) example.
