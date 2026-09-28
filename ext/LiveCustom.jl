#=
Public customization API of the live view: own panels (`add_panel!`), controls (`add_controls!`)
and tools (`add_tool!`), placed by the layout via its slots (see `AbstractLiveLayout`), and
`retrace!`, which solves again after a change from code. Nothing here runs per frame: the panels
are updated after full solves only, the widgets and keys react to their events.
=#

#=
Panels
=#

function add_panel!(f, gui::LiveView, title::AbstractString; select::Bool = false)
    p = _add_user_panel!(f, gui, String(title), select)
    push!(gui.custom.panels, p)
    _register_widgets!(gui, p.layout)
    # The live view has solved already, the panel shows its result right away
    _refresh_user_panel!(gui, p)
    return p.layout
end

"""
    _add_user_panel!(f, gui, title::String, select::Bool) -> _UserPanel

The slot of [`add_panel!`](@ref) in the layout of the `gui`: creates the place of a panel named
`title`, builds its content via `f(layout)` (see `_user_panel`) and returns the panel. `select`
shows it at once in layouts that show one panel at a time.
"""
_add_user_panel!(_, gui::LiveView, ::String, ::Bool) = _slot_error(gui, "place for panels")

"""
Returns the `_UserPanel` named `title` with the content `layout`, built by `f(layout)`: its result
is the `update` of the panel if it is a `Function`, otherwise the panel has no update, e.g. for the
last plot of a `do` block.
"""
_user_panel(f, title::String, layout::GridLayout) =
    _UserPanel(title, layout, _update_function(f(layout)), nothing)

_update_function(update::Function) = update
_update_function(_) = nothing

"""
    _panel_shown(gui, p) -> Bool

Whether the panel `p` of the `gui` is shown, i.e. updated after a solve: always, unless the layout
hides panels, see `_mark_panel_stale!`.
"""
_panel_shown(::LiveView, _) = true

"""
    _mark_panel_stale!(gui, p)

Marks the hidden panel `p` of the `gui` as stale, i.e. to be updated once it is shown; nothing in
layouts that show all panels, see `_panel_shown`.
"""
_mark_panel_stale!(::LiveView, _) = nothing

"""
Whether the detectors and beams of the `gui` hold the result of a full solve, which the `update`
of a panel may read: not while a solve runs in the background or after a preview solve.
"""
_results_valid(gui::LiveView) = !_running(gui) && !gui.preview

"""
    _update_user_panels!(gui)

Calls the `update` of all panels of [`add_panel!`](@ref) that are shown, after a full solve (see
`_apply!`); hidden panels are marked stale instead. Nothing is done without such panels.
"""
function _update_user_panels!(gui::LiveView)
    for p in gui.custom.panels
        _panel_shown(gui, p) ? _update_user_panel!(gui, p) : _mark_panel_stale!(gui, p)
    end
    return nothing
end

"""Updates the panel `p` if it is shown and the last solve is complete, otherwise marks it stale."""
function _refresh_user_panel!(gui::LiveView, p::_UserPanel)
    if _panel_shown(gui, p) && _results_valid(gui)
        _update_user_panel!(gui, p)
    else
        _mark_panel_stale!(gui, p)
    end
    return nothing
end

"""Calls the `update` of the panel `p`, errors are logged once per distinct message."""
_update_user_panel!(gui::LiveView, p::_UserPanel) = _run_update!(gui, p, p.update)
_run_update!(::LiveView, ::_UserPanel, ::Nothing) = nothing
function _run_update!(gui::LiveView, p::_UserPanel, update::Function)
    try
        update(gui)
        p.last_error = nothing
    catch e
        p.last_error = _log_once(e, p.last_error, "update of the panel \"$(p.title)\"")
    end
    return nothing
end

#=
Controls
=#

function add_controls!(f, gui::LiveView, title::AbstractString)
    layout = _controls_slot!(gui, String(title))
    f(layout)
    _on_controls_added!(gui)
    _register_widgets!(gui, layout)
    return layout
end

"""
    _controls_slot!(gui, title::String) -> GridLayout

The slot of [`add_controls!`](@ref) in the layout of the `gui`: the layout of a new, titled group
of widgets. `_on_controls_added!(gui)` is called once the widgets are built.
"""
_controls_slot!(gui::LiveView, ::String) = _slot_error(gui, "place for controls")
_on_controls_added!(::LiveView) = nothing

#=
Widgets that take the keyboard
=#

"""
    _register_widgets!(gui, layout::GridLayout)

Registers the `Textbox`es and `Menu`s in the `layout` (e.g. of [`add_controls!`](@ref)), such that
they take the keyboard like those of the live view: while a box is focused or a menu is open, the
keys of the 3D view are ignored (see `_typing`) and the camera does not take the keys, see
`_keep_keyboard!`.
"""
function _register_widgets!(gui::LiveView, layout::GridLayout)
    foreach(b -> _register_widget!(gui, b), _blocks!(Any[], layout))
    return nothing
end

_register_widget!(::LiveView, _) = nothing
function _register_widget!(gui::LiveView, tb::Textbox)
    push!(gui.custom.boxes, tb)
    push!(gui.controls.listeners, on(_ -> _keep_keyboard!(gui), tb.focused))
    return nothing
end
function _register_widget!(gui::LiveView, m::Menu)
    push!(gui.custom.menus, m)
    push!(gui.controls.listeners, on(_ -> _keep_keyboard!(gui), m.is_open))
    return nothing
end

"""Whether a registered textbox of the `parts` is focused or a registered menu is open."""
_custom_typing(parts::_UserParts) =
    any(tb -> tb.focused[], parts.boxes) || any(m -> m.is_open[], parts.menus)

#=
Tools
=#

function add_tool!(f, gui::LiveView, name::AbstractString; icon::Symbol = :object,
        key::Union{Nothing, Keyboard.Button} = nothing, toggle::Bool = false,
        tooltip::AbstractString = name)
    # An unknown icon throws in all layouts, such that code works with any layout
    _icon(icon)
    _check_key(gui, key, name)
    w = _tool_widget(gui, Val(toggle), _key_label(name, key), icon, _key_label(tooltip, key))
    _connect_tool!(gui, f, w, Val(toggle), String(name))
    _connect_tool_key!(gui, w, key, Val(toggle), String(name))
    return w
end

"""
    _tool_widget(gui, toggle::Val, label, icon, tooltip)

The slot of [`add_tool!`](@ref) in the layout of the `gui`: creates a button (`Val(false)`, with
`clicks`) or a toggle (`Val(true)`, with `active`, initially off), shown with the `label` or the
`icon` and the `tooltip`, which name the key of the tool, if any.
"""
_tool_widget(gui::LiveView, ::Val, _, _, _) = _slot_error(gui, "place for tools")

_key_label(s::AbstractString, ::Nothing) = String(s)
_key_label(s::AbstractString, key::Keyboard.Button) = "$s ($(_key_name(key)))"

"""The name of the `key` as shown in labels, e.g. `2` for `Keyboard._2`."""
_key_name(key::Keyboard.Button) = lstrip(string(key), '_')

"""Calls `f(gui)` after each click of the button `w`, or `f(gui, active)` after the toggle `w`."""
function _connect_tool!(gui::LiveView, f, w, ::Val{false}, name::String)
    last_error = Ref{Union{Nothing, String}}(nothing)
    push!(gui.controls.listeners, on(_ -> _call_tool!(() -> f(gui), last_error, name), w.clicks))
    return nothing
end
function _connect_tool!(gui::LiveView, f, w, ::Val{true}, name::String)
    last_error = Ref{Union{Nothing, String}}(nothing)
    push!(gui.controls.listeners, on(v -> _call_tool!(() -> f(gui, v), last_error, name), w.active))
    return nothing
end

function _call_tool!(call, last_error::Ref, name::String)
    try
        call()
        last_error[] = nothing
    catch e
        last_error[] = _log_once(e, last_error[], "tool \"$name\"")
    end
    return nothing
end

"""
Connects the `key` of the tool `name` with its widget `w`: a press in the 3D view clicks the
button or switches the toggle, unless a textbox or menu takes the keyboard, see `_typing`.
"""
_connect_tool_key!(::LiveView, _, ::Nothing, ::Val, ::String) = nothing
function _connect_tool_key!(gui::LiveView, w, key::Keyboard.Button, toggle::Val, name::String)
    gui.custom.keys[key] = name
    push!(gui.controls.listeners, on(events(gui.ax.scene).keyboardbutton, priority = 200) do event
        (event.action == Keyboard.press && event.key == key) || return Consume(false)
        gui.controls.ignore_keys() && return Consume(false)
        _press!(w, toggle)
        return Consume(true)
    end)
    return nothing
end

_press!(w, ::Val{false}) = (w.clicks[] += 1; nothing)
_press!(w, ::Val{true}) = (w.active[] = !w.active[]; nothing)

"""
Keys of the live view (see `LiveView.jl`, the kinematic controls in `LiveInteraction.jl`) with
their bindings; a tool can not use them, see `_check_key`. The keys of Makie's `Camera3D` are read
from the camera, see `_camera_keys`. `+` and `-` are typed characters, their keys are listed for
the usual layouts.
"""
const _LIVE_VIEW_KEYS = (
    Keyboard.t => "trace (t)",
    Keyboard.escape => "Esc (cancel a solve, clear the inspection, deselect)",
    Keyboard.c => "clipping (c, Shift+c)",
    Keyboard.p => "add a clip plane (p)",
    Keyboard.delete => "remove the clip plane (Delete)",
    Keyboard._1 => "the source markers (1)",
    Keyboard.g => "zoom to the selection (g)",
    Keyboard.v => "the spectator mode (v)",
    Keyboard.h => "the help (h)",
    Keyboard.m => "the move/rotate mode (m)",
    Keyboard.z => "undo (Ctrl+z)",
    Keyboard.y => "redo (Ctrl+y, Ctrl+Shift+z)",
    Keyboard.backspace => "reset the pose (Backspace)",
    Keyboard.up => "the key steps (arrow keys)",
    Keyboard.down => "the key steps (arrow keys)",
    Keyboard.left => "the key steps (arrow keys)",
    Keyboard.right => "the key steps (arrow keys)",
    Keyboard.page_up => "the key steps (Page Up/Down)",
    Keyboard.page_down => "the key steps (Page Up/Down)",
    Keyboard.left_shift => "Shift (large key steps, Shift+c)",
    Keyboard.right_shift => "Shift (large key steps, Shift+c)",
    Keyboard.left_control => "Ctrl (undo, redo)",
    Keyboard.right_control => "Ctrl (undo, redo)",
    Keyboard.left_super => "Cmd (undo, redo)",
    Keyboard.right_super => "Cmd (undo, redo)",
    Keyboard.minus => "the step size (-)",
    Keyboard.equal => "the step size (+)",
    Keyboard.kp_add => "the step size (+)",
    Keyboard.kp_subtract => "the step size (-)",
)

"""The keyboard buttons of a key binding of Makie, e.g. `Keyboard.left_control & Mouse.left`."""
_keyboard_buttons(b::Keyboard.Button) = (b,)
_keyboard_buttons(op::Union{Makie.And, Makie.Or}) =
    (_keyboard_buttons(op.left)..., _keyboard_buttons(op.right)...)
_keyboard_buttons(op::Makie.Exclusively) = Tuple(filter(b -> b isa Keyboard.Button, collect(op.x)))
_keyboard_buttons(bs::Union{Tuple, AbstractVector}) = Tuple(Iterators.flatten(map(_keyboard_buttons, bs)))
# e.g. `Mouse.left`, `Not(...)` or `false` (never pressed)
_keyboard_buttons(_) = ()

"""The keys of the camera of the 3D view of the `gui`, `binding => key`, e.g. `:forward_key`."""
function _camera_keys(gui::LiveView)
    controls = cameracontrols(gui.ax.scene).controls
    return [k => b for k in keys(controls) for b in _keyboard_buttons(Makie.to_value(controls[k]))]
end

"""
    _key_binding(gui, key) -> Union{Nothing, String}

The existing binding of the `key` in the `gui`: of the live view (see `_LIVE_VIEW_KEYS`), of the
`select_modifier` of the controls, of another tool or of the camera; `nothing` if it is free.
"""
function _key_binding(gui::LiveView, key::Keyboard.Button)
    for (k, binding) in _LIVE_VIEW_KEYS
        k == key && return binding
    end
    key in _keyboard_buttons(gui.controls.select_modifier) && return "the select modifier"
    haskey(gui.custom.keys, key) && return "the tool \"$(gui.custom.keys[key])\""
    for (name, k) in _camera_keys(gui)
        k == key && return "Camera3D ($name)"
    end
    return nothing
end

_check_key(::LiveView, ::Nothing, _) = nothing
function _check_key(gui::LiveView, key::Keyboard.Button, name::AbstractString)
    binding = _key_binding(gui, key)
    isnothing(binding) && return nothing
    throw(ArgumentError("the key $(_key_name(key)) of the tool \"$name\" is taken by $binding"))
end

#=
Solving after a change from code
=#

retrace!(gui::LiveView) = retrace!(Returns(nothing), gui)

function retrace!(f, gui::LiveView)
    # Like a slider, see `_connect_sliders!`: `f` may change any object
    _change!(f, gui.controls, nothing)
    update_render!(gui.controls.h)
    _on_change!(gui, nothing)
    return nothing
end

#=
Compact layout: panels below the detector panels in `fig[1, 2]`, controls above the status row,
tools in the tool row
=#

"""
Returns the grid of the panels in `fig[1, 2]` of the compact layout, created without detector
panels: then the rows below the 3D view (sliders, status and tool row, controls) span the new
column, like with detector panels.
"""
_panel_grid!(gui::CompactView) = _panel_grid!(gui, gui.layout.panels)
_panel_grid!(::CompactView, grid::GridLayout) = grid
function _panel_grid!(gui::CompactView, ::Nothing)
    root = gui.fig.layout
    for c in copy(root.content)
        c.span.rows.start > 1 && (root[c.span.rows, 1:2] = c.content)
    end
    grid = gui.layout.panels = GridLayout(root[1, 2])
    colsize!(root, 1, Relative(0.6))
    return grid
end

function _add_user_panel!(f, gui::CompactView, title::String, ::Bool)
    grid = _panel_grid!(gui)
    nrows, ncols = size(grid)
    # Below all panels, over all columns of their grid
    row = isempty(grid.content) ? 1 : nrows + 1
    box = GridLayout(grid[row, 1:ncols])
    Label(box[1, 1], title; font = :bold, halign = :left, tellwidth = false)
    layout = GridLayout(box[2, 1])
    rowgap!(box, 4)
    return _user_panel(f, title, layout)
end

function _controls_slot!(gui::CompactView, title::String)
    root = gui.fig.layout
    row = _GLB.gridcontent(gui.layout.status_row).span.rows.start
    _GLB.insertrows!(root, row, 1)
    box = GridLayout(root[row, 1:_GLB.ncols(root)])
    Label(box[1, 1], title; font = :bold)
    layout = GridLayout(box[1, 2])
    # Keeps the controls left-aligned, like the tool row
    Label(box[1, 3], ""; tellwidth = false)
    colgap!(box, 10)
    return layout
end

"""Returns `n` positions at the end of the tool row of the compact `layout`, before its filler."""
function _tool_slots!(layout::CompactLayout, n::Int)
    row = layout.tool_row
    k = _GLB.ncols(row)
    row[1, k + n] = only(_GLB.contents(row[1, k]))
    return [row[1, k + i] for i in 0:(n - 1)]
end

_tool_widget(gui::CompactView, ::Val{false}, label, _, _) =
    Button(only(_tool_slots!(gui.layout, 1)); label)

function _tool_widget(gui::CompactView, ::Val{true}, label, _, _)
    toggle_pos, label_pos = _tool_slots!(gui.layout, 2)
    toggle = Toggle(toggle_pos; active = false)
    Label(label_pos, label)
    return toggle
end

#=
App layout: panels as tabs of the analysis dock, controls as sections of the left sidebar, tools
in the toolbar group `:user` before "Help"
=#

_mark_panel_stale!(gui::AppView, p::_UserPanel) = (push!(gui.layout.tabs.stale, p); nothing)
_results_valid(gui::AppView) = gui.layout.tabs.hits_valid && !_running(gui) && !gui.preview

"""Hides the blocks of the collapsed `part`, e.g. those added to it since it was collapsed."""
function _hide_collapsed!(part::_AppPart)
    part.shown || foreach(Makie.hide!, _blocks!(Any[], part.grid))
    return nothing
end

"""
Adds the panel as a tab of the dock: the content is built into the shown tab, then the previously
active tab is shown again, unless `select` (or there was none). The panel is updated when its tab
is opened, see `_refresh_tab!`.
"""
function _add_user_panel!(f, gui::AppView, title::String, select::Bool)
    layout = gui.layout
    tabs = layout.tabs
    old = tabs.active
    content = _add_dock_panel!(layout, title; icon = :chart)
    i = tabs.active
    p = try
        _user_panel(f, title, content)
    finally
        (select || old == 0) || _select_tab!(layout, old)
        _hide_collapsed!(layout.dock)
    end
    tabs.panels[i] = p
    return p
end

function _refresh_tab!(gui::AppView, p::_UserPanel)
    _results_valid(gui) || return nothing
    _update_user_panel!(gui, p)
    delete!(gui.layout.tabs.stale, p)
    return nothing
end

_controls_slot!(gui::AppView, title::String) = _add_sidebar_section!(gui, :left, title)
_on_controls_added!(gui::AppView) = _hide_collapsed!(gui.layout.left)

_tool_kwargs(layout::AppLayout) =
    (; _icon_theme(layout.theme)..., size = _TOOLBAR_BUTTON, icon_size = _TOOLBAR_ICON)

_tool_widget(gui::AppView, ::Val{false}, _, icon::Symbol, tooltip::String) =
    _IconButton(_add_toolbar_entry!(gui, :user); icon, tooltip, _tool_kwargs(gui.layout)...)

_tool_widget(gui::AppView, ::Val{true}, _, icon::Symbol, tooltip::String) =
    _IconToggle(_add_toolbar_entry!(gui, :user); icon, tooltip, active = false,
        _tool_kwargs(gui.layout)...)
