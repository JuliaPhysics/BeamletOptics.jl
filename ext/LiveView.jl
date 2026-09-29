using Makie: Figure, Axis, Label, SliderGrid, GridLayout, DataAspect, Relative, colsize!, colgap!,
             heatmap!, autolimits!, limits!, Button, Toggle, Textbox, Menu, rowgap!, linkxaxes!,
             hidexdecorations!, hidespines!
import InteractiveUtils

"""
    _SolveJob

A solve of the systems (or a computation of the detector fields) of a `LiveView` in a background
task, see `_solve!` and `_run!`.

# Fields

- `task`: runs `_compute` and returns its result
- `done`: notified when `task` ends (or the wait of `_run!` times out)
- `sinks`: the progress outputs of `task`, see `BMO._ProgressSink`: one per source, then one per
  detector panel
- `anchors`: the position of the progress window of each sink, i.e. of its source or detector
- `apply`: shows the result of `task` in the live view, called on the render task
- `obj`: the moved object, or `nothing`
- `timing`: the duration field of the live view (e.g. `:solve_time`) that a cancelled job updates
- `t0`: `time()` at the start
- `shown`: the loop whose progress window is shown, `(; k, t0, t, count)`: the index of its sink,
  its start and the time and count when its window appeared (`k = 0` before the first window)
"""
mutable struct _SolveJob
    task::Task
    done::Base.Event
    sinks::Vector{BMO._ProgressSink}
    anchors::Vector{Point3f}
    apply::Function
    obj::Any
    timing::Symbol
    t0::Float64
    shown::@NamedTuple{k::Int, t0::Float64, t::Float64, count::Int}
end

"""
    AbstractLiveLayout

Layout of a [`LiveView`](@ref), i.e. where its widgets are placed in the figure, see the
`layout` kwarg of [`live_view`](@ref). The logic of the live view (tracing, clipping, measuring,
camera, inspection) is shared by all layouts: it works on the observables of the widgets, e.g.
`clicks` of a button, `active` of a toggle and `stored_string` of a textbox, never on their type
or position. A layout `L` holds its own state (e.g. the slots of its sidebars) and is stored in
`gui.layout` of the `LiveView{L}`.

# Interface

A subtype `L <: AbstractLiveLayout` has the field `theme`, the color tokens of the `theme` kwarg
(see `_APP_THEMES`), in which the floating cards and the progress window are drawn in all layouts,
and implements

- `_build_layout(layout::L, fig, spec) -> NamedTuple`: creates the widgets in the `Figure` `fig`,
  which `_figure(layout, size)` created. `spec` holds the inputs of `live_view`: `specs` (detector
  panels, see `_panel_specs`), `slider_specs`, `labels`, `lighting`, `view_cube`, `auto_trace`,
  `clip_beams`, `orthographic`, `show_sources` and `view_specs`. The result has the fields
  - `ax`: the `LScene` of the 3D view, with `studio_lighting!` applied, and `cube`: its view cube
    or `nothing`
  - `panels`: the `DetectorPanel`s of `spec.specs`, `sliders`: a `SliderGrid` of
    `spec.slider_specs` or `nothing`, `status`: the `Label` of the status line
  - the built-in tools by role, e.g. `trace_button` (anything with `clicks::Observable{Int}`) or
    `auto_trace_toggle` (anything with `active::Observable{Bool}`, initialized from `spec`), created
    by `_build_tools(layout, spec)`, see "Tools" below
  - `info`: the `Label` of the last solve, the number of rays and the projection, or `nothing`

  `live_view` collects the widgets in a `_LayoutWidgets`, whose typed fields reject a missing or
  wrong widget at construction.
- `_build_menus(layout::L, w, options, views_options) -> (; menu, views_menu)`: called with the
  result `w` of `_build_layout` once the movable objects are known; the component `menu` (a `Menu`
  with the `options`, or `nothing`) and the `views_menu` (a `Menu` with the `views_options`)

and optionally, with defaults for any layout,

- `_figure(layout::L, size)`: the `Figure`, and `_default_size(layout::L)`: its size unless given
- `_connect_layout!(gui::LiveView{L})`: connects the widgets that only the layout has, e.g.
  collapsing; its listeners belong in `gui.controls.listeners`
- where the controls of the selected object are shown, i.e. the widgets that
  [`card_actions`](@ref) and [`card_rows`](@ref) declare, built by the same code for any host (see
  `_AbstractCard`): by default on the floating card next to the object (`_ComponentCard`), which
  holds the step box. A layout that shows them elsewhere, e.g. the app layout on a card docked in
  its inspector (`_DockedCard`), returns `false` from `_selection_card_shown(gui)` and refreshes
  them in `_refresh_inspector!(gui; force)`, called with the cards (see `_update_inspector!`);
  `_step_box(layout, w, card)` returns its `Textbox` of the keyboard step, `_layout_boxes(gui)`
  the textboxes that take the keyboard (see `_typing`). Pinned cards float by default; a layout
  that shows them elsewhere, e.g. the app layout below its inspector, implements `_pin!(gui, obj)`,
  `_unpin!(gui, obj)` and `_is_pinned(gui, obj)`, see `_toggle_pin!`; if it also lets a pinned
  card float, `_float!(gui, obj)` and `_dock!(gui, obj)`, whose floating cards get the button that
  docks them (`_card_tools!(layout, c)`).
- `_outside_view(gui)`: `true` while the mouse is over a part of the layout whose clicks must not
  reach the controls of the 3D view, e.g. the sidebars of the app layout (none by default)
- hooks called by the shared logic: `_on_solved!(gui)` after a solve is shown,
  `_on_selected!(gui)` after the selection changed, `_on_clipping!(gui)` after clipping was
  switched, `_on_clip_planes_changed!(gui)` after a clip plane was added or removed,
  `_on_hidden!(gui)` after objects were hidden or shown and `_on_pinned!(gui)` after a card was
  pinned or unpinned
- which detector panels are computed and how their results are shown, for layouts that show only
  some of them: `_computed_panels(gui, preview)`, `_shown_panels(gui)` and
  `_apply_panel!(gui, p, field; coarse, preview)` (all panels by default), with the hooks
  `_on_solve_started!(gui)` and `_on_applied!(gui)` around a solve, see `LiveDock.jl`
- `_show_hint(gui)`: how a hidden object is shown again, for the status line
- colors of the 3D view: by default from the tokens of the `theme` (see `LiveLayout.jl`),
  `_clip_plane_color(layout)`, `_marker_stroke(layout)` (the outline of the handles of sources,
  clip planes and measured points), `_beam_style(layout, beam)` (default kwargs of `live_render!`
  of a source) and `_theme_render!(layout, h)` for the rendered objects
- slots for additional parts: `_add_toolbar_entry!(gui, group)`, `_add_sidebar_section!(gui,
  side, title)` and `_add_dock_panel!(gui, title)`, which return the `GridPosition` or
  `GridLayout` to place widgets in, see `AppLayout`
- the places of the public customization API (see `LiveCustom.jl`), without which it throws for
  the layout: `_add_user_panel!(f, gui, title, select)` for [`add_panel!`](@ref),
  `_controls_slot!(gui, title)` for [`add_controls!`](@ref) and `_tool_widget` (see "Tools")
  for [`add_tool!`](@ref). A layout that
  hides panels implements `_panel_shown(gui, p)` and `_mark_panel_stale!(gui, p)` and updates a
  stale panel once it is shown, see `_UserPanel`

# Tools

The built-in tools (buttons and toggles such as trace, clipping or export) are declared once in
`_BUILTIN_TOOLS` (see `_ToolSpec`). A layout implements `_tool_widget(layout, group, toggle::Val,
label, icon, tooltip, active)`, which creates a tool of a `group` in the place of the layout, e.g.
a labelled button in a row or an icon in a toolbar; `_build_tools(layout, spec)` builds all tools
that `_has_tool(layout, Val(role))` selects (by default those with a field in `_LayoutWidgets`)
this way, and [`add_tool!`](@ref) its tools with the group `:user`.

The layouts of `live_view` are `CompactLayout` and `AppLayout`.
"""
abstract type AbstractLiveLayout end

"""
    _UserPanel

A panel added via [`add_panel!`](@ref): its `title`, the `layout` of its content and the
`update(gui)` returned by the builder (or `nothing`), which runs after each full solve while the
panel is shown, see `_update_user_panels!`. A layout that hides panels (the tabs of the app layout)
marks a hidden panel stale instead and updates it once it is shown, see `_mark_panel_stale!`.
"""
mutable struct _UserPanel
    const title::String
    const layout::GridLayout
    const update::Union{Nothing, Function}
    last_error::Union{Nothing, String}
end

"""
    _UserControls

Controls added via [`add_controls!`](@ref): their `title`, the `layout` of their widgets and the
`update(gui)` returned by the builder (or `nothing`), which runs once when they are added and after
each full solve, like the `update` of a `_UserPanel`, but never stale, since controls are always
shown.
"""
mutable struct _UserControls
    const title::String
    const layout::GridLayout
    const update::Union{Nothing, Function}
    last_error::Union{Nothing, String}
end

"""
    _UserParts()

The parts of a live view added via the customization API (see `LiveCustom.jl`): the `panels`
(`_UserPanel`s), the `controls` (`_UserControls`), the `boxes` (`Textbox`es) and `menus` of added
controls and panels, which take the keyboard like those of the layout (see `_typing`), and the
`keys` of added tools with their names, see [`add_tool!`](@ref).
"""
struct _UserParts
    panels::Vector{_UserPanel}
    controls::Vector{_UserControls}
    boxes::Vector{Textbox}
    menus::Vector{Menu}
    keys::Dict{Keyboard.Button, String}
end

_UserParts() = _UserParts(_UserPanel[], _UserControls[], Textbox[], Menu[], Dict{Keyboard.Button, String}())

"""
    _TraceState

Tracing of a `LiveView`: the systems are solved after each change if `auto[]` (the `active`
observable of the auto trace toggle) is `true`, otherwise only via the trace button, the key `t` or
by switching auto tracing on. `stale` is `true` if the beams and detector panels do not match the
current poses of the objects; the alpha of the beam plots before dimming is stored in
`beam_alphas`. If solving takes longer than `budget` [s], the systems are solved once the movement
pauses for `idle_delay` [s] (`pending`, the moved `pending_obj`, the time of the `last_change`);
`solve_time`, `panel_time` and `preview_time` are the durations of the last solve, panel update
and preview solve [s], `coarse` is `true` while the panels show a preview on a coarse grid.

While moving, beam groups are solved only for their rendered beams if `preview_enabled`; `preview`
is `true` from such a solve (of the moved `preview_obj`) until the full solve. A solve that takes
longer than `budget` continues in the background as `job`, the `progress` window shows its loops
after `progress_delay` [s]. The duration fields are named by `_SolveJob.timing`.
"""
Base.@kwdef mutable struct _TraceState
    auto::Observable{Bool}
    budget::Float64
    idle_delay::Float64
    preview_enabled::Bool
    progress::_ProgressOverlay
    progress_delay::Float64
    stale::Bool = false
    beam_alphas::IdDict{Any, Any} = IdDict{Any, Any}()
    solve_time::Float64 = 0.0
    panel_time::Float64 = 0.0
    preview_time::Float64 = 0.0
    pending::Bool = false
    pending_obj::Any = nothing
    last_change::Float64 = 0.0
    coarse::Bool = false
    preview::Bool = false
    preview_obj::Any = nothing
    job::Union{Nothing, _SolveJob} = nothing
end

"""
    _ClipState

Clip planes of a `LiveView`: the `planes`, applied if `enabled`, to the beams as well if `beams`.
`size` is the edge length of the outline of new planes, fixed at construction, since the bounding
boxes of clipped plots only cover their visible part.
"""
Base.@kwdef mutable struct _ClipState
    planes::Vector{LiveClipPlane} = LiveClipPlane[]
    enabled::Bool = true
    size::Float64
    beams::Bool
end

"""
    _MeasureState

Beam inspection and measuring of a `LiveView`: the inspected point of a beam (see
`_inspect_beam`) and its marker `inspection_plot`; up to two measured `points` `(; point, obj)`,
the `result` `(; distance, angle)` and its `plots`.
"""
Base.@kwdef mutable struct _MeasureState
    inspection::Any = nothing
    inspection_plot::Union{Nothing, AbstractPlot} = nothing
    points::Vector{Any} = Any[]
    result::Any = nothing
    plots::Vector{AbstractPlot} = AbstractPlot[]
end

"""
    _CameraState

Camera tools of a `LiveView`: the `home` view `(eye, lookat, up)`, taken at the first tick
(`home_set`), the saved `views` and the running `animation` of a camera transition.
"""
Base.@kwdef mutable struct _CameraState
    home::NTuple{3, Vector{Float64}} = (zeros(3), zeros(3), zeros(3))
    home_set::Bool = false
    views::Vector{Pair{String, NTuple{3, Vector{Float64}}}}
    animation::Any = nothing
end

"""
    _CardState

Cards of a `LiveView`: the card of the `selection` next to the selected object, which shows the
rows and actions declared for it (see [`card_rows`](@ref)); `all` cards, including the pinned
ones; the listeners that keep the camera from the cards (`shield`, see `_shield_cards!`).
"""
Base.@kwdef mutable struct _CardState
    selection::_ComponentCard
    all::Vector{_ComponentCard} = [selection]
    shield::Vector{Any} = Any[]
end

"""
    _ObjectState

The objects of a `LiveView` as the selection shows them: the movable objects of the component
`menu` (the option `i` selects `menu[i]`), the `hidden` objects (rendered objects, i.e. leaves of
groups, whose plots are invisible, see the action "hide" of the cards), the `opacity` of objects
set via their card (see `_set_opacity!`), and the automatic `names` of objects without a label
with the `counters` of their running indices per type, see `_name_objects!`. `inspected` is shown
on the card of the selection, but not selected for moving: a system or an object that is not
movable, see `_inspect!`; it and `controls.selected[]` exclude each other.
"""
Base.@kwdef mutable struct _ObjectState
    menu::Vector{Any} = Any[]
    hidden::Base.IdSet{Any} = Base.IdSet{Any}()
    opacity::IdDict{Any, Any} = IdDict{Any, Any}()
    names::IdDict{Any, String} = IdDict{Any, String}()
    counters::Dict{String, Int} = Dict{String, Int}()
    inspected::Any = nothing
end

"""
    _LayoutWidgets

The widgets of a `LiveView` that its layout creates (see `_build_layout` and `_build_menus`) and
the shared logic uses. They are typed by what the logic uses, since the layouts use different
widgets: buttons have `clicks::Observable{Int}`, toggles `active::Observable{Bool}`.

- buttons: `trace_button`, `export_button` (see [`export_changes`](@ref)), `home_button`,
  `save_view_button`, `show_all_button` (or `nothing`)
- toggles: `auto_trace_toggle`, `clip_beams_toggle`, `orthographic_toggle`, `sources_toggle` (the
  markers of the movable sources, also the key `1`), `measure_toggle`
- `step_box`: the `Textbox` of the keyboard step (replaced when another card becomes the card of
  the selection, see `_use_card!`); `menu`: the component menu (or `nothing`),
  `views_menu`: the saved views; `view_cube`: the view cube of the 3D view (or `nothing`)
- `info`: the label of the last solve, the number of rays and the projection (or `nothing`)
"""
Base.@kwdef mutable struct _LayoutWidgets
    trace_button::Any
    auto_trace_toggle::Any
    clip_beams_toggle::Any
    orthographic_toggle::Any
    sources_toggle::Any
    measure_toggle::Any
    export_button::Any
    home_button::Any
    save_view_button::Any
    show_all_button::Any
    step_box::Textbox
    menu::Union{Nothing, Menu}
    views_menu::Menu
    view_cube::Union{Nothing, ViewCube}
    info::Union{Nothing, Label} = nothing
end

"""
    LiveView

Interactive window returned by [`live_view`](@ref). The `Figure` is stored in `fig`, the `LScene`
of the 3D view in `ax`, the `KinematicController` in `controls`, the detector `panels`, the
`status` line and the `sliders` (or `nothing`). Use `display` to show the window and `close` to
remove the controls and the view cube.

The state of the shared logic is grouped by concern: `trace` (`_TraceState`), `clip`
(`_ClipState`), `measure` (`_MeasureState`), `camera` (`_CameraState`), `cards` (`_CardState`) and
`objects` (`_ObjectState`); the widgets that the layout creates are in `widgets`
(`_LayoutWidgets`). The export button prints the changed poses as Julia code, see
[`export_changes`](@ref), and copies them to the clipboard if `export_clipboard` is `true`. The
objects of the `extras` kwarg are rendered, selectable and movable, but not part of any system,
see `_live_render_extras!`. Panels, widgets and tool keys added via the customization API are in
`custom`, see `LiveCustom.jl`.

The type parameter `L` is the type of the `layout`, see `AbstractLiveLayout`: `CompactView` and
`AppView` are the `LiveView`s of `live_view(...; layout = :compact)` and `layout = :app`.
"""
Base.@kwdef mutable struct LiveView{L <: AbstractLiveLayout}
    fig::Figure
    ax::LScene
    pairs::Vector{Pair{BMO.AbstractSystem, Any}}
    system_handles::Vector{SystemRenderHandle}
    beam_handles::Vector{AbstractRenderHandle}
    controls::KinematicController
    panels::Vector{Any}
    status::Label
    sliders::Union{Nothing, SliderGrid}
    on_change::Function
    last_error::Union{Nothing, String} = nothing
    labels::IdDict{Any, String}
    export_clipboard::Bool = true
    extras::SystemRenderHandle
    custom::_UserParts = _UserParts()
    trace::_TraceState
    clip::_ClipState
    measure::_MeasureState = _MeasureState()
    camera::_CameraState
    cards::_CardState
    objects::_ObjectState
    widgets::_LayoutWidgets
    # state of the layout, e.g. the slots of the app layout, see `AbstractLiveLayout`
    layout::L
end

# Infers the layout type like the constructor of a non-parametric type
LiveView(fields...) = LiveView{typeof(last(fields))}(fields...)


"""
    display(gui::LiveView; screen_config...)

Displays the window of the `gui`. With GLMakie, the live view gets its own window, rendered without
SSAO and with up to 60 fps: SSAO (e.g. enabled globally via `GLMakie.activate!(ssao = true)`)
multiplies the frame time of a live view with large meshes, and GLMakie's default of 30 fps makes
rotating the view sluggish. The `screen_config` kwargs of the backend override these defaults.

A new window is used instead of reusing the current one, since GLMakie fails to reuse a window
whose screen configuration (e.g. `ssao`) and size both change ("Binding freed Texture").
"""
function Base.display(gui::LiveView; screen_config...)
    if _multi_light_backend() # i.e. GLMakie
        # A figure can only be shown in one screen: an open window of the gui is reused, unless new
        # screen settings are given
        old = Makie.getscreen(gui.fig.scene)
        if !isnothing(old) && isopen(old)
            isempty(screen_config) && return old
            close(old)
        end
        screen = Makie.current_backend().Screen(; ssao = false, framerate = 60.0, screen_config...)
        return display(screen, gui.fig)
    end
    return display(gui.fig; screen_config...)
end

function Base.show(io::IO, gui::LiveView)
    print(io, "LiveView(", length(gui.pairs), " systems, ", length(gui.panels), " detector panels)")
end

function Base.close(gui::LiveView)
    _cancel_solve!(gui)
    close(gui.controls)
    isnothing(gui.widgets.view_cube) || close(gui.widgets.view_cube)
    foreach(_hide_card!, gui.cards.all)
    return nothing
end

"""Logs the error `e` of the `source` once per distinct message, returns the new last message."""
function _log_once(e, last_error, source::String)
    msg = sprint(showerror, e)
    if msg != last_error
        @error "live_view: $source failed" exception = (e, catch_backtrace())
    end
    return msg
end

"""
    live_view(system => beam, ...; kwargs...)
    live_view(system, beam; kwargs...)

Opens a complete interactive window for one or several pairs of `system` and `beam`. All systems
and beams are live-rendered into the same `LScene`, see [`live_render!`](@ref), and can be moved
with the [`kinematic_controls!`](@ref). After each change, all `Detector`s are emptied, all systems
are solved again and the beams and detector panels are updated. Returns a `LiveView` with the
fields `fig`, `ax`, `controls`, `panels`, `status` and `sliders`. Use `display(gui)` to show the
window and `close(gui)` to remove the controls.

Additional context, e.g. an optomechanical assembly from a CAD file, is passed via `extras`: these
objects are rendered, selected, moved, hidden and exported like the components, but they are not
part of any system, i.e. never traced and without cost in the solves; see "Extras and opacity".
Plots added via `render!(gui.ax, ...)` are only drawn. The status line shows the pose of the moved
object and its change since the window was opened.

# Component card, component menu and export

The selected object opens a card next to its bounding box in the 3D view, connected to it by a
line, which follows the camera and the object: right of the box, or left of, below or above it if
there is no room, always inside the view and never over the view cube or another card. It is
hidden without a selection and while a menu is open. "pin" in its head keeps the card with its
object, independent of the selection, e.g. to watch or type the poses of several objects; the
selection then gets a new card. "unpin" closes a pinned card. The head shows the label of the
object and its actions:

- "hide" makes the plots of the object invisible and clears the selection. On a hidden object
  selected in the menu, the button reads "show" and shows it again. Hidden objects can not be
  selected in the 3D view, but are still traced.
- For a clip plane, "flip" and "remove" instead, like `Shift+c` and `Delete`.

Below, `x`, `y`, `z` [mm] show the position of the object, `Enter` in a box moves the object to the
typed absolute coordinate. `rx`, `ry` and `rv` [mrad] rotate the object about the red, green and
blue axis of the controls, like the keys in the rotate mode. Each input is recorded in the undo
history, invalid inputs are reported in the status line. The widgets of a pinned card act on its
object, also if another object is selected. Only the card of the selection has, below its rows,
the `step` box of the keyboard step, e.g. `250 nm` or `50 µrad`, where the unit selects the move
or rotate mode, the "Move"/"Rotate" control, which shows and sets the mode of the controls (it
follows the key `m` and vice versa), and a "Properties" part, collapsed by default, which lists
the properties of the object (see [`properties`](@ref)), the same rows as in the inspector of the
app layout. The state of "Properties" is kept while the card follows the selection. "–" in the
head collapses the card to its head, "+" expands it again. Clicks and drags on the card neither
select objects nor move the camera, and while a box of the card has the focus, the keys of the 3D
view are ignored. Further rows by type, see [`card_rows`](@ref): the number of rays of a source,
the mode (`auto`, `spot`, `intensity`) and the log color scale of the panel of a `Detector` ("no
panel" without one) and the opacity of mechanics. The card, the menus and the buttons take the
colors of the `theme`.

Objects without a `labels` entry are named by their type and a running index, e.g. "Mirror 1" or
"Clip plane 2", in the card, the status line and the menus, like in the object tree of the app
layout. The rows of an own type are added by a method of [`card_rows`](@ref), see the page "Live
view widgets" of the documentation.

Beside the component cards, the card of a system shows the number of its objects, the number of
rays and the duration of the last solve. It is shown, without a gizmo and without a selection
(`controls.selected[]` stays `nothing`), by selecting the system entry ("System 1", ...) in the
component menu. An object that is not movable is shown on its card in the same way instead of
being selected, and its pose boxes reject inputs with a message in the status line. `Esc`, a click
on empty space or the selection of an object closes such a card.

The row below the status line holds a menu of the systems and all movable objects (each system
entry followed by its objects, the objects of a group indented after the group, without clip
planes), which selects an object like a click in the 3D view, and "show all", which shows all
hidden objects. The last cell of the status row is the info label: the duration of the last solve
(or of the preview), the number of rays and the projection.

The "Export" button prints the changed poses as Julia code to `stdout` and copies it to the
clipboard, see [`export_changes`](@ref).

# Extras and opacity

`extras = [housing => (; color = :lightblue), mount]` adds objects without optical function, e.g.
a `MeshDummy` of an STL file, each an `AbstractObject` (or a group of them) with optional kwargs
of its [`render!`](@ref) call. They are not part of any system: never traced, moving them does not
solve the systems. Otherwise they act like the components: a click in the 3D view (where they do
not cover the components behind them, which are picked first) or in the component menu or object
tree selects them, they are moved with the controls if their kinematic trait allows it (static
ones are only rendered and hidden), hidden via "hide" or the eye of the tree, in which they are
listed under "Extras", and their moves are exported by [`export_changes`](@ref). An extra must not
be an object of a system.

The card of mechanics, i.e. a `NonInteractableObject` (e.g. `MeshDummy`) or an
`IntersectableObject`, in a system or among the extras, has an "opacity" slider (0-100 %,
initially the opacity as rendered, e.g. 5 % for a color with alpha 0.05). It scales the alpha of
the plots of the object, including its feature edges, and does not change the optics. Objects below
100 % are drawn with order independent transparency (`transparency = true`), opaque ones without
it, which is cheaper and keeps the correct depth; switching rebuilds the render objects of the
object once (a few ms for a mesh with 1 M triangles). At 0 % the object is hidden like via "hide",
a larger opacity shows it again. The opacity is kept while the object is hidden; showing an object
hidden at 0 % restores its initial opacity.

Mechanics below 50 % opacity (as rendered or set via the slider) are not selected by a click in
the 3D view: a click into the empty space inside or over a transparent housing selects what lies
behind it, e.g. a component, or clears the selection. Select them via the object tree or the
component menu instead. From 50 % on, a click selects them like any component.

# Beam inspection and measuring

A click on a rendered beam (within 6 px) that does not hit a component marks the point on the beam
and shows its position [mm], the direction of the beam, the geometric and optical path length
(Σ n·L) from the source [mm] and, for Gaussian beamlets, the radius `w` and the radius of curvature
`R` in the status line. Components take precedence over beams. `esc` or a click elsewhere removes
the marker.

With the "measure" toggle on, two clicks on components or beams show the distance between the
positions of the components or the points of the beams [mm], and the angle between the optical
axes (local y-axes) of two components, with a dashed line between the points. A third click starts
a new measurement, switching the toggle off clears it.

# Camera tools

The key `g` zooms to the selected object, or to all systems if nothing is selected, while the view
direction is kept. "home" restores the view when the window was shown. The "views" menu sets one
of the `views`, "save view" adds the current view as `"view n"` and prints it as an entry of the
`views` kwarg, e.g. `"view 1" => ([0.1, -0.2, 0.3], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0])`.

# Adaptive tracing

If solving the systems takes longer than `trace_budget`, the objects still follow the mouse and the
keys immediately, while the beams are dimmed. The systems are solved once the movement pauses for
`idle_delay`. Likewise, detector panels that take longer than `trace_budget` show a preview on a
coarse grid while objects are moved, which is refined once the movement pauses.

With `preview = true`, beam groups rendered with `render_every > 1` are solved only for their
rendered beams while objects are moved, and the titles of the detector panels end with
"(preview)". The full beam group is solved once the movement pauses for `idle_delay`. The
`trace_budget` applies to the preview solve while moving. `on_change` is only called after full
solves.

# Long solves

A solve, or the computation of the detector panels, that takes longer than `progress_delay` runs
in the background: the camera can still be moved, the beams are dimmed and the status line shows
"tracing". The loops that show a progress bar in the terminal, i.e. the tracing of a beam group
and the field of a detector panel, show a small progress window in the 3D view next to their source
or detector once they have run for `progress_delay`, with the remaining time; the terminal bar is
not drawn meanwhile. Moving a component or a source, a slider and `Esc` cancel the solve after the
current beam, `t` is ignored until it is done.

# Manual tracing

With `auto_trace = false`, the systems are not solved after each change, which is useful for
systems that take long to solve. Objects and sliders still update the 3D view, while the beams are
dimmed and the status line shows that they are outdated. The systems are solved by the
`Trace (t)` button below the 3D view or the key `t`. The toggle next to the button switches auto
tracing on or off, switching it on solves the systems if they are outdated. The initial solve
always runs.

# Detector panels

By default, one panel per `Detector` of all systems is shown next to the 3D view. The panel shows
the intensity (`:intensity`) for Gaussian beamlet hits and the spot diagram (`:spot`) otherwise,
together with the optical power or the number of hits in its title. The intensity is cropped
automatically around the beam, unless `x_min`, `x_max`, `z_min` and `z_max` are given.

The subtitle of each panel shows its metrics, the centroid is marked by a red cross: the number of
hits, the centroid, the RMS radius and the geometric radius of a spot diagram, or the power, the
centroid, the 1/e² radii along x and z from the second moments and the peak of the intensity. The
following options of the panel kwargs are not passed to `intensity`:

- `colorscale = :linear`: `:log` shows `log10` of the intensity, with a floor of 1e-4 times the
  maximum
- `colorrange = nothing`: fixed color range of the intensity, in `log10` units for `:log`
- `history = false`: adds an axis below the panel with the power (or the number of hits, black)
  and the centroid x (red) and z (blue) over the last 300 full solves
- `profiles = false`: adds an axis below the panel with the intensity along x (red) and z (blue)
  through the centroid

# Clip planes

Clip planes cut the 3D view open, e.g. to look into a housing. Only the side of a plane its normal
points to stays visible. The key `p` adds a plane through the selected object, or through the
`lookat` point of the camera if nothing is selected, with its normal along the view direction.
A plane is selected via the purple handle at its center and moved and rotated like a component,
its normal is the green axis. Moving a plane does not solve the systems. The keys apply as follows:

| key       | action                                   |
|:----------|:-----------------------------------------|
| `p`       | add a clip plane and select it           |
| `Delete`  | remove the selected clip plane           |
| `c`       | switch clipping on or off (all planes)   |
| `Shift+c` | flip the selected clip plane             |

Makie supports at most 8 clip planes. The markers of the sources and planes and the controls are
never clipped, the beams only with `clip_beams = true` or the "clip beams" toggle below the 3D
view. The selection box of a partly clipped component only covers its visible part.

The key `g` zooms to the selection, see "Camera tools".

# App layout

With `layout = :app`, the same live view is arranged as an application window, all keys and mouse
actions in the 3D view are unchanged:

- toolbar (icons with tooltips): trace and auto trace, home, fit (`g`), views, save view,
  orthographic, clipping (`c`), clip beams, sources (`1`), measure, export, the toggles of the
  sidebars and the dock, help (`h`)
- left sidebar: the object tree and, below it, the sliders ("Parameters"). The tree lists each
  system with its objects (groups with their objects, collapsed by default), then the `extras`
  under "Extras", the sources and the clip planes. A click on a name selects the object like a click in the 3D view, and a
  selection in the 3D view highlights its row. The eye hides or shows an object, a group or a
  whole system, the eye in the title of the tree shows all objects again. A click on the name of
  a system row shows the card of the system in the inspector, without selecting anything (the
  expander still folds the row); a click on an object that is not movable shows its card in the
  same way. Objects without a `labels` entry are named by their type and a running index, e.g.
  "Lens 2", also in the status line.
- right sidebar ("Properties"): the card of the selected object, docked instead of floating next
  to it: its name and type, the actions of the card (e.g. "hide", or "flip" and "remove" for a
  clip plane) and a pin, which pins a card to the object; below, the rows of the card (see
  [`card_rows`](@ref), e.g. the pose, the ray slider of a source or the panel options of a
  detector); then the step box, the mode and the properties of the object (see
  [`properties`](@ref)). The pinned cards are docked below, one below the other, each with its own
  head (icon, label, actions, float button, pin and chevron). The sidebar does not scroll: if the
  docked cards do not fit, the older ones collapse to their heads. The float button of a docked
  card moves it into the 3D view, where it floats next to its object like a pinned card of the
  compact layout; the dock button in its head moves it back to the end of the docked cards. A
  card keeps its collapsed state when it moves; pinned again after it was unpinned, it starts
  docked.
- analysis dock below the 3D view: a tab per detector panel, a click on a tab shows its panel.
  Only the panel of the active tab is computed after a solve, the other panels are computed when
  their tab is opened; a collapsed dock computes none. Panels with `history = true` still record
  every full solve. The history and profiles axes are shown beside the panel.
- status bar: the status line and the duration of the last solve, the number of rays and the
  projection

The sidebars and the dock can be collapsed via the toolbar, the 3D view then takes their space.
The component menu of the compact layout is replaced by the tree. In the compact layout, the widgets are at
fixed positions of `gui.fig`, e.g. the detector panels in `gui.fig[1, 2]`, next to which users
may add their own axes.

# Own panels, controls and tools

[`add_panel!`](@ref) adds an own panel (below the detector panels, or a tab of the dock of the app
layout), [`add_controls!`](@ref) own widgets (a row above the status row, or a section of the left
sidebar) and [`add_tool!`](@ref) a button or toggle, optionally with a key (in the tool row, or the
toolbar). [`retrace!`](@ref) solves again after a change from code, e.g. from such a widget.
Widgets of a thing in the scene belong on its card, see [`card_rows`](@ref), own widget types on
cards and in controls use [`card_input`](@ref) and [`card_show!`](@ref).

# Keyword args

- `layout = :compact`: arrangement of the widgets, `:compact` or `:app`, see "App layout"
- `theme = :light`: colors, `:light` or `:dark`, of the whole window (background, cards, buttons,
  menus, status row and progress window) in both layouts. For contrast on a dark 3D view, `:dark`
  draws the rays, the markers and the dark materials of the render look (detectors, polarizers) in
  lighter colors; `:light` keeps the colors of the default look.
- `size`: size of the figure, by default `(1400, 800)` for `:compact` and `(1600, 950)` for `:app`
- `auto_trace = true`: solves the systems after each change, otherwise only on request, see
  "Manual tracing"
- `detectors = :auto`: all `Detector`s of all systems. Alternatively a vector of `pd`,
  `pd => mode` or `pd => (mode, kwargs)`, where `mode` is `:auto`, `:spot` or `:intensity` and
  `kwargs` are passed to `intensity`, e.g. `(; n = 200, x_min = -1e-3, x_max = 1e-3, ...)`,
  except the panel options, e.g. `(; colorscale = :log, history = true)`, see "Detector panels".
  An empty vector disables the panels.
- `on_change = (gui, obj) -> nothing`: called after each full solve with the moved object, or
  `nothing` after a slider change, i.e. not after preview solves, see "Adaptive tracing"
- `sliders = []`: vector of `"label" => (range, callback)` or `"label" => (range, callback, startvalue)`.
  The `callback` is called with the new value, then the systems are solved again.
- `system_kwargs = (;)`: passed to `live_render!` of each system
- `extras = []`: objects that are rendered and can be moved and hidden, but are not traced, a
  vector of `obj` or `obj => render_kwargs`, e.g. `[housing => (; transparency = true, color =
  RGBAf(0.7, 0.8, 0.9, 0.05))]`, see "Extras and opacity"
- `beam_kwargs = Dict()`: `beam => kwargs` passed to `live_render!` of the beam, by default
  `(; render_every = 5)` for beam groups
- `movable_sources = true`: shows an orange marker at each source, i.e. the beam or beam group of
  each pair, with which the source can be selected and moved like the components
- `show_sources = true`: initial visibility of the source markers, which can be switched with the
  "sources" toggle below the 3D view or the key `1`
- `labels = Dict()`: `obj => "name"` for the status line, the titles of the detector panels, the
  component menu and the variable names of [`export_changes`](@ref)
- `trace_budget = 0.03`: [s] duration of a solve or panel update, above which tracing is deferred
  or the panels show a preview, see "Adaptive tracing"
- `idle_delay = 0.2`: [s] pause of the movement after which deferred tracing runs
- `clip_planes = []`: initial clip planes, a vector of `point => normal`, e.g.
  `[[0, 0.1, 0] => [0, 1, 0]]`, see "Clip planes"
- `clip_beams = false`: clips the beams as well, can be switched with the "clip beams" toggle
- `view_cube = true`: shows a view cube in the top right corner of the 3D view, a click on a
  face, edge or corner switches to the corresponding standard view, see [`view_cube!`](@ref)
- `orthographic = false`: starts the 3D view with orthographic instead of perspective projection,
  can be switched with the "orthographic" toggle below the 3D view
- `lighting = :studio`: lighting rig of the 3D view, see [`studio_lighting!`](@ref), `:none`
  keeps the default lights of Makie
- `edges = nothing`: draws the feature edges of the components, by default depending on the look,
  see [`set_render_look`](@ref) and [`render!`](@ref). An `edges` entry of `system_kwargs` takes
  precedence.
- `preview = true`: solves beam groups only for their rendered beams while moving, see
  "Adaptive tracing"
- `views = []`: saved views of the "views" menu, a vector of `"name" => (eye, lookat, up)`, see
  "Camera tools"
- `progress_delay = 0.5`: [s] duration after which a solve continues in the background and a loop
  shows its progress window, see "Long solves"
- all other kwargs are passed to [`kinematic_controls!`](@ref), e.g. `fine_step`, `plane_normal`
  or `rotation_axis`
"""
function live_view(
        pairs::Pair{<:BMO.AbstractSystem}...;
        size = nothing,
        layout::Symbol = :compact,
        theme::Symbol = :light,
        auto_trace::Bool = true,
        detectors = :auto,
        on_change = (gui, obj) -> nothing,
        sliders = [],
        system_kwargs = (;),
        beam_kwargs = Dict(),
        movable_sources = true,
        show_sources::Bool = true,
        labels = Dict(),
        trace_budget = 0.03,
        idle_delay = 0.2,
        clip_planes = [],
        clip_beams::Bool = false,
        view_cube::Bool = true,
        orthographic::Bool = false,
        lighting::Symbol = :studio,
        edges::Union{Nothing, Bool} = nothing,
        preview::Bool = true,
        views = [],
        progress_delay::Real = 0.5,
        extras = [],
        kwargs...
    )
    isempty(pairs) && throw(ArgumentError("live_view requires at least one system => beam pair"))
    ps = Pair{BMO.AbstractSystem, Any}[p for p in pairs]
    # several beams may share a system, which is rendered once
    systems = unique(objectid, first.(ps))
    extra_specs = _extra_specs(extras, systems)
    specs = _panel_specs(detectors, systems)
    slider_specs = [_slider_spec(s) for s in sliders]
    clip_specs = _clip_plane_specs(clip_planes)
    view_specs = _view_specs(views)

    lay = _live_layout(layout, theme)
    fig = _figure(lay, something(size, _default_size(lay)))
    w = _build_layout(lay, fig, (; specs, slider_specs, labels, lighting, view_cube, auto_trace,
        clip_beams, orthographic, show_sources, view_specs))
    ax = w.ax
    # Pose, keyboard step and hide button of the selected object, next to it in the 3D view
    card = _ComponentCard(fig, lay, _card_z(1))

    # `edges` is only passed if given, i.e. custom `render!` methods of user objects do not need to
    # accept it
    sys_kw = isnothing(edges) ? system_kwargs : (; edges, system_kwargs...)
    system_handles = SystemRenderHandle[live_render!(ax, sys; sys_kw...) for sys in systems]
    beam_handles = AbstractRenderHandle[]
    for beam in last.(ps)
        default = beam isa BMO.AbstractBeamGroup ? (; render_every = 5) : (;)
        kw = get(beam_kwargs, beam, default)
        # The planes of the beams are set explicitly by `_apply_clip_planes!`, see `clip_beams`
        push!(beam_handles, live_render!(ax, beam; _beam_style(lay, beam)..., kw...,
            clip_planes = Plane3f[]))
    end

    # A single controller for all systems, otherwise several controllers would compete for events
    handles = reduce(vcat, [h.handles for h in system_handles]; init = ObjectRenderHandle[])
    # The extras are moved and selected like the objects of the systems, but never traced
    extras_handle = _live_render_extras!(ax, extra_specs)
    # Size of the scene (the systems and the visible extras), before any clip plane shrinks the
    # bounding boxes
    extent = _scene_extent((system_handles..., extras_handle))
    if movable_sources
        # Markers of the sources, scaled to the size of the systems
        marker_size = 0.08 * extent
        for src in unique(objectid, last.(ps))
            BMO._is_static(src) || push!(handles,
                _live_render_source!(ax, src; size = marker_size, strokecolor = _marker_stroke(lay)))
        end
    end
    append!(handles, extras_handle.handles)
    parent = IdDict{BMO.AbstractObject, BMO.AbstractObject}()
    foreach(h -> merge!(parent, h.parent), (system_handles..., extras_handle))
    combined = SystemRenderHandle(ax, first(systems), handles, parent)
    # Colors of the render look that the theme of the layout replaces, e.g. of dark detectors
    _theme_render!(lay, combined)
    gui_ref = Ref{LiveView}()
    # Moving a clip plane or an extra does not solve the systems, see `_on_moved!`
    change = function (obj)
        gui = gui_ref[]
        _on_moved!(gui, obj)
        _update_inspector!(gui)
        return nothing
    end
    # Typing into a textbox or the search of the component menu must not trigger the controls
    controls = kinematic_controls!(ax, combined; on_change = change,
        ignore_keys = () -> isassigned(gui_ref) && _typing(gui_ref[]), kwargs...)

    labels_dict = IdDict{Any, String}(labels)
    entries = _menu_entries(controls)
    menus = _build_menus(lay, w, _menu_options(labels_dict, entries), _views_options(view_specs))
    widgets = _LayoutWidgets(; w.trace_button, w.auto_trace_toggle, w.clip_beams_toggle,
        w.orthographic_toggle, w.sources_toggle, w.measure_toggle, w.export_button, w.home_button,
        w.save_view_button, w.show_all_button, step_box = _step_box(lay, w, card), menus.menu,
        menus.views_menu, view_cube = w.cube, w.info)
    trace = _TraceState(; auto = w.auto_trace_toggle.active, budget = trace_budget, idle_delay,
        preview_enabled = preview, progress = _ProgressOverlay(ax, lay.theme), progress_delay)
    gui = LiveView(; fig, ax, pairs = ps, system_handles, beam_handles, controls, w.panels,
        w.status, w.sliders, on_change, labels = labels_dict, extras = extras_handle, trace,
        clip = _ClipState(; size = 1.2 * extent, beams = clip_beams),
        camera = _CameraState(; views = view_specs), cards = _CardState(; selection = card),
        objects = _ObjectState(; menu = Any[first.(entries)...]), widgets, layout = lay)
    gui_ref[] = gui
    # Names of the objects without a label, e.g. for the object tree, see `_name_objects!`
    _name_objects!(gui)
    # Objects must not change while a solve in the background traces them
    controls.before_change = () -> _cancel_solve!(gui)
    for (point, normal) in clip_specs
        _add_clip_plane!(gui, point, normal; select = false)
    end
    isnothing(gui.sliders) || _connect_sliders!(gui, last.(slider_specs))
    _connect_trace!(gui)
    _connect_clip_planes!(gui)
    _connect_tools!(gui)
    _connect_inspection!(gui)
    _connect_camera!(gui)
    _connect_cards!(gui)
    push!(controls.listeners, on(v -> v == gui.clip.beams || _set_clip_beams!(gui, v),
        gui.widgets.clip_beams_toggle.active))
    _connect_projection!(gui, orthographic)
    _connect_sources!(gui)
    _connect_layout!(gui)
    # The info label and the colors of the controls, shared by all layouts
    _connect_theme!(gui)
    _resolve!(gui, nothing)
    # Initial view from the Front-Right-Top corner, in which the labels of the view cube read
    # correctly. Only set once, later changes of the view, e.g. via `set_view`, are kept.
    cam = cameracontrols(ax.scene)
    lookat = Vector{Float64}(cam.lookat[])
    dist = norm(Vector{Float64}(cam.eyeposition[]) .- lookat)
    o, up = _region_view((1, -1, 1))
    set_view(ax, lookat .+ dist .* o, lookat, up)
    # Replaced by the view at the first tick, i.e. when the window is shown
    gui.camera.home = _current_view(gui)
    return gui
end

live_view(system::BMO.AbstractSystem, beam; kwargs...) = live_view(system => beam; kwargs...)
