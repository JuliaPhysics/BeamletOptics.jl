#=
Cards of the live view: the pose inspector, the card of the selection and the pinned cards,
their widgets and moving them with the mouse; the card widget is in LiveCard.jl
=#

#=
Pose inspector
=#

"""Labels of the boxes of the pose inspector: position [mm] and rotations about the axes [mrad]."""
const _POSE_FIELDS = ("x [mm]", "y [mm]", "z [mm]", "rx [mrad]", "ry [mrad]", "rv [mrad]")

"""Gizmo axes (see `_axis_vectors`) of the rotation boxes of the pose inspector."""
const _POSE_AXES = (:x, :y, :v)

"""Label colors of the pose inspector: the colors of the gizmo axes for the rotations."""
const _POSE_COLORS = (:black, :black, :black, :red, :green, :blue)

"""Returns `true` if a textbox or a menu of the `gui` takes keyboard input."""
function _typing(gui::LiveView)
    any(c -> any(tb -> tb.focused[], _card_boxes(c)), gui.cards.all) && return true
    any(tb -> tb.focused[], _layout_boxes(gui)) && return true
    # the widgets of `add_controls!` and `add_panel!`, see `_register_widgets!`
    _custom_typing(gui.custom) && return true
    return _menu_open(gui)
end

"""
Returns the textboxes of the `gui` outside of the cards, e.g. of the inspector of the app layout,
which take the keyboard like those of the cards, see `_typing`. None by default.
"""
_layout_boxes(::LiveView) = ()

"""Shows `s` in the textbox `tb` without triggering its listeners, `""` shows the placeholder."""
function _set_box!(tb::Textbox, s::String)
    tb.stored_string.val = isempty(s) ? nothing : s
    tb.displayed_string[] == s || (tb.displayed_string[] = s)
    return nothing
end

"""
    _update_inspector!(gui; force = false)

Shows the values of the objects of all cards of the `gui` in their widgets, e.g. the position of
an object in its pose boxes, see `_refresh_card!`, and, in layouts with an inspector, refreshes it,
see `_refresh_inspector!`. Focused textboxes keep the typed text, unless `force`.
"""
function _update_inspector!(gui::LiveView; force::Bool = false)
    foreach(c -> _refresh_card!(gui, c; force), gui.cards.all)
    _refresh_inspector!(gui; force)
    return nothing
end

"""
Refreshes the parts of the layout of the `gui` that show the selected object besides the cards,
e.g. the inspector of the app layout; nothing by default.
"""
_refresh_inspector!(::LiveView; force::Bool = false) = nothing

"""Returns `true` if the constraints of `obj` allow a move by `Δ`, see `kinematic_controls!`."""
function _move_allowed(ctrl::KinematicController, obj, Δ)
    haskey(_constraints_of(ctrl, obj), :move) || return true
    allowed = _allowed_axes(ctrl, obj, :move)
    isempty(allowed) && return iszero(Δ)
    A = hcat(_axis_vectors(ctrl, obj, allowed)...)
    return norm(Δ - A * (pinv(A) * Δ)) <= 1e-12 + 1e-9 * norm(Δ)
end

"""
    _apply_pose_input!(gui, obj, k, s)

Applies the input `s` of the pose box `k` of a card (see `_POSE_FIELDS`) to its object `obj`: for
`k ≤ 3`, moves it to the absolute position x, y or z [mm], otherwise rotates it about the red,
green or blue axis of the controls [mrad], like a key step in the rotate mode. The change is
recorded in the undo history and solved like a key step. An invalid input only shows a message in
the status line, and so does an input for an object that is not movable, e.g. an inspected one (see
`_inspect!`).
"""
_apply_pose_input!(gui::LiveView, ::Nothing, ::Int, _) = _update_inspector!(gui; force = true)
function _apply_pose_input!(gui::LiveView, obj, k::Int, s)
    ctrl = gui.controls
    if !_is_movable(ctrl, obj)
        gui.status.text[] = "$(_label(gui, obj)) is not movable, its pose can not be changed"
        _update_inspector!(gui; force = true)
        return nothing
    end
    x = isnothing(s) ? nothing : tryparse(Float64, strip(s))
    if isnothing(x) || !isfinite(x)
        gui.status.text[] = "invalid input \"$(something(s, ""))\" for $(_POSE_FIELDS[k]), enter a number"
        return nothing
    end
    P0, R0 = _pose(obj)
    if k <= 3
        P = Vector{Float64}(P0)
        P[k] = x / 1e3
        if !_move_allowed(ctrl, obj, P - P0)
            gui.status.text[] = "$(_label(gui, obj)) can not move along $(_POSE_FIELDS[k][1]), see the constraints"
            _update_inspector!(gui; force = true)
            return nothing
        end
        _change!(ctrl, obj) do
            translate_to3d!(obj, P)
            # `translate_to3d!` moves by `P - position`, which may round
            r = P - Vector{Float64}(position(obj))
            iszero(r) || translate3d!(obj, r)
        end
    else
        sym = _POSE_AXES[k - 3]
        if !(sym in _allowed_axes(ctrl, obj, :rotate))
            gui.status.text[] = "$(_label(gui, obj)) can not rotate about this axis, see the constraints"
            _update_inspector!(gui; force = true)
            return nothing
        end
        iszero(x) || _change!(() -> rotate3d!(obj, only(_axis_vectors(ctrl, obj, (sym,))), x / 1e3),
            ctrl, obj)
    end
    P1, R1 = _pose(obj)
    ctrl.last_key_step = nothing
    _push_history!(ctrl, obj, P0, R0, P1, R1)
    _request_update!(ctrl)
    _update_inspector!(gui; force = true)
    return nothing
end

"""
Camera3D takes the keyboard (WASD etc.) only after a click on the background of the 3D view, which
a click on a widget, e.g. the orthographic toggle, undoes. Here the keyboard stays with the camera
of the `gui` unless a textbox or a menu takes the input, see `_typing`.
"""
function _keep_keyboard!(gui::LiveView)
    cam = cameracontrols(gui.ax.scene)
    selected = !_typing(gui)
    cam.selected[] == selected || (cam.selected[] = selected)
    return nothing
end
#=
Component card
=#

"""
Returns the object of the card `c` of the `gui`: the pinned object, the selected (or inspected, see
`_shown_object`) object for the card of the selection (`gui.cards.selection`), `nothing` for a
spare card.
"""
function _card_object(gui::LiveView, c::_ComponentCard)
    c.pinned && return c.obj
    return c === gui.cards.selection ? _shown_object(gui) : nothing
end

"""
    _update_cards!(gui)

Shows the card of the selection (`gui.cards.selection`) next to the selected (or inspected, see
`_inspect!`) object, unless it has a pinned card, and the pinned cards of the `gui` next to their objects; see `_update_card!`. The cards are
placed in this order, each off the view cube and the cards before, so that none covers another:
the cards that the mouse moved to their `spot` first, where they stay (see `_drag_cards!`), then
the card of the selection, at its object; a pinned card without room is collapsed to its head. All
cards are hidden while a menu is open, whose options they would cover. Called every frame, which
moves the cards with the camera and the objects.
"""
function _update_cards!(gui::LiveView)
    menu = _menu_open(gui)
    sel = _shown_object(gui)
    obstacles = _obstacles(gui.widgets.view_cube)
    shown = menu || !_selection_card_shown(gui) || any(c -> c.pinned && c.obj === sel, gui.cards.all) ?
        nothing : sel
    target(c) = c === gui.cards.selection ? shown : c.pinned && !menu ? c.obj : nothing
    order = [gui.cards.selection; filter(c -> c !== gui.cards.selection, gui.cards.all)]
    for moved in (true, false), c in order
        isnothing(c.spot) == moved || _update_card!(gui, c, target(c), obstacles)
    end
    return nothing
end

"""
Returns `true` if the layout of the `gui` shows the card of the selection (`gui.cards.selection`) next to the
selected object; layouts that show the card elsewhere, e.g. docked in the inspector of the app
layout, return `false`. Pinned cards are shown in all layouts.
"""
_selection_card_shown(::LiveView) = true

"""
    _update_card!(gui, c, obj, obstacles)

Shows the card `c` for `obj` (or hides it for `nothing`): builds its declared widgets when it gets
another object (see `_build_content!`) and shows the values of `obj` when it or its pose changes,
then moves the card next to the bounding box of `obj` (see `_card_position`), off the `obstacles`
(see `_avoid`), to which it adds its rectangle, and connects it to `obj` by a line. A pinned card
that finds no room is collapsed to its head (`auto_collapsed`) until there is room again. A card
that the mouse moved stays at its `spot` instead (see `_spot_position`). Only changed values update
the layout.
"""
_update_card!(::LiveView, c::_ComponentCard, ::Nothing, ::Vector{Rect2f}) = _hide_card!(c)
function _update_card!(gui::LiveView, c::_ComponentCard, obj, obstacles::Vector{Rect2f})
    pose = _card_pose(obj)
    if c.pose === nothing || c.pose[1] !== obj
        _build_content!(gui, c, obj)
    end
    # Another object or a new pose, e.g. after a move by a program
    if c.pose === nothing || c.pose[1] !== obj || c.pose[2] != pose
        c.pose = (obj, pose)
        _refresh_card!(gui, c)
    end
    corners = _card_corners(gui, c, obj)
    _update!(c.title.text, _label(gui, obj))
    _show_kind!(c, obj)
    scene = gui.ax.scene
    view = Rect2f(Makie.viewport(scene)[])
    sel = _screen_rect(scene, corners, obj)
    c.auto_collapsed = false
    _fit_properties!(c)
    size = _card_size(c)
    p = isnothing(c.spot) ? _avoid(_card_position(sel, size, view), size, view, obstacles) :
        _spot_position(c.spot, size, view)
    if isnothing(c.spot) && c.pinned && !c.collapsed && _covers(p, size, obstacles)
        c.auto_collapsed = true
        size = _card_size(c)
        p = _avoid(_card_position(sel, size, view), size, view, obstacles)
    end
    _arrange_card!(c, p)
    rect = _card_rect(p, size)
    push!(obstacles, rect)
    _update!(c.link, [_link_anchor(scene, corners, obj), Point2f(minimum(rect) .+ size ./ 2)])
    _update!(c.scene.visible, true)
    return nothing
end

"""
    _build_content!(gui, c, obj)

Builds the widgets of the card `c` of the `gui` for its new object `obj` from the declarations
[`card_actions`](@ref) (in `c.actions`) and [`card_rows`](@ref) (one layout per row in `c.rows`),
see `_add_cell!`. If the declarations have the same layout as those of the widgets on the card,
e.g. for another mirror, the widgets are kept and only take the new declarations. The card is any
host of the declarations, see `_AbstractCard`, e.g. the floating card, whose new widgets come
before the mouse shield of the cards (see `_on_content_built!`), or the docked card of the app.
"""
function _build_content!(gui::LiveView, c::_AbstractCard, obj)
    actions, rows = _declarations(c, obj)
    key = (_layout_key(actions), _layout_key(rows))
    declared = CardWidget[_declared_widgets(actions)..., _declared_widgets(rows)...]
    if key == c.content_key
        c.widgets = [(b, w) for ((b, _), w) in zip(c.widgets, declared)]
        return nothing
    end
    _clear_content!(c)
    for (j, w) in enumerate(actions)
        _add_cell!(gui, c, c.actions[1, j], w)
    end
    for (i, row) in enumerate(rows)
        layout = GridLayout(c.rows[i, 1]; _row_attributes(c)...)
        for (j, cell) in enumerate(row.cells)
            _add_cell!(gui, c, layout[1, j], cell)
        end
    end
    c.content_key = key
    _on_content_built!(gui, c)
    return nothing
end

"""
    _declarations(c, obj) -> (actions, rows)

The declarations of the widgets of `obj` on the card `c`: [`card_actions`](@ref) and
[`card_rows`](@ref), which a host may extend, e.g. the docked card of the app layout.
"""
_declarations(::_AbstractCard, obj) = (card_actions(obj), card_rows(obj))

# The new widgets of a floating card come before the mouse shield of the cards
_on_content_built!(gui::LiveView, ::_ComponentCard) = _shield_cards!(gui)
_on_content_built!(::LiveView, ::_AbstractCard) = nothing

"""
    _add_cell!(gui, c, pos, cell)

Adds a cell of a declaration to the card `c` of the `gui` at the grid position `pos`: a text as a
`Label`, a [`CardWidget`](@ref) as a block of its type with the attributes of the card, see
`_cell_attributes`. The inputs of the block call `on` of the declaration (see `_on_input!`), and a
textbox takes the keyboard like the others of the card.
"""
function _add_cell!(::LiveView, c::_AbstractCard, pos, text::String)
    push!(c.blocks, Label(pos, text; halign = :left, _card_style(c, Label)...))
    return nothing
end
function _add_cell!(gui::LiveView, c::_AbstractCard, pos, w::CardWidget)
    b = w.type(pos; _cell_attributes(c, w)...)
    push!(c.blocks, b)
    push!(c.widgets, (b, w))
    i = length(c.widgets)
    _fix_caret!(c, b)
    _track_textbox!(gui, c, b)
    _listen_input!(gui, c, i, BMO.card_input(b))
    return nothing
end

function _track_textbox!(gui::LiveView, c::_AbstractCard, tb::Textbox)
    push!(c.textboxes, tb)
    push!(c.listeners, on(_ -> _keep_keyboard!(gui), tb.focused))
    return nothing
end
_track_textbox!(::LiveView, ::_AbstractCard, _) = nothing

_listen_input!(::LiveView, ::_AbstractCard, ::Int, ::Nothing) = nothing
function _listen_input!(gui::LiveView, c::_AbstractCard, i::Int, obs::Observable)
    push!(c.listeners, on(v -> _on_input!(gui, c, i, v), obs))
    return nothing
end

"""
    _on_input!(gui, c, i, v)

Applies the input `v` of the declared widget `i` of the card `c` to the object of the card (see
`_card_object`) with `on` of its declaration (and solves again for `solve = true`), then shows the
new values on all cards. Ignored while the card shows new values (`refreshing`).
"""
function _on_input!(gui::LiveView, c::_AbstractCard, i::Int, v)
    c.refreshing && return nothing
    w = c.widgets[i][2]
    _apply_input!(gui, w.on, _card_object(gui, c), v, w.solve)
    _update_inspector!(gui)
    _update_cards!(gui)
    return nothing
end
_apply_input!(gui::LiveView, on, obj, v, solve::Bool) = _apply_on!(gui, on, obj, v, Val(solve))
_apply_input!(::LiveView, ::Nothing, _, _, ::Bool) = nothing
_apply_input!(::LiveView, _, ::Nothing, _, ::Bool) = nothing
_apply_input!(::LiveView, ::Nothing, ::Nothing, _, ::Bool) = nothing
_apply_on!(gui::LiveView, on, obj, v, ::Val{false}) = (on(gui, obj, v); nothing)
# The input changes the optics, see `solve` of `CardWidget`: like a move, via the `on_change` of
# the controls, which callers of the live view may extend
function _apply_on!(gui::LiveView, on, obj, v, ::Val{true})
    _change!(() -> on(gui, obj, v), gui.controls, obj)
    gui.controls.on_change(obj)
    return nothing
end

"""
    _refresh_card!(gui, c; force = false)

Shows the values of the object of the card `c` (see `_card_object`) in its declared widgets, see
`value` of [`CardWidget`](@ref), and in the part of the selection of a floating card, see
`_refresh_selection_part!`; a focused textbox keeps the typed text, unless `force`. Only if the
widgets were built for this object (`c.pose`), see `_update_card!`.
"""
function _refresh_card!(gui::LiveView, c::_AbstractCard; force::Bool = false)
    obj = _card_object(gui, c)
    (isnothing(obj) || c.pose === nothing || c.pose[1] !== obj) && return nothing
    c.refreshing = true
    try
        for (b, w) in c.widgets
            _refresh_widget!(b, w.value, gui, obj; force)
        end
    finally
        c.refreshing = false
    end
    _refresh_selection_part!(gui, c, obj)
    return nothing
end
_refresh_widget!(b, value, gui::LiveView, obj; force::Bool = false) = _show_value!(b, value(gui, obj), force)
_refresh_widget!(_, ::Nothing, ::LiveView, _; force::Bool = false) = nothing

"""
    _card_corners(gui, c, obj)

Returns the corners of the bounding box of the object `obj` of the card `c` of the `gui`: of the
selection box for the card of the selected object, which the controls keep up to date. A pinned
card and the card of an inspected object (see `_inspect!`) take the bounding box of the plots of
`obj` once (see `_card_bbox`) and move it with the pose of `obj` (see `key` and `_card_pose`),
since the plots follow a move only after they are rendered again.
"""
function _card_corners(gui::LiveView, c::_ComponentCard, obj)
    ctrl = gui.controls
    (c.pinned || obj !== ctrl.selected[]) || return ctrl.box_obs[]
    if c.key === nothing || c.key[1] !== obj
        c.corners = _box_corners(_card_bbox(ctrl, obj))
        c.key = (obj, _card_pose(obj))
    end
    return _moved_corners(c.corners, c.key[2], _card_pose(obj))
end

"""
    _card_pose(obj)

The pose of the object `obj` of a card (see `_pose`), by which the card notices a move; `nothing`
for a system, whose card is refreshed after the moves of its objects like any card, see
`_update_inspector!`.
"""
_card_pose(obj) = _pose(obj)
_card_pose(::BMO.AbstractSystem) = nothing

"""
Returns the `corners` of a bounding box taken in the pose `(P0, R0)`, moved to the pose `(P, R)`;
unchanged for an object without a pose, see `_card_pose`.
"""
function _moved_corners(corners, (P0, R0), (P, R))
    T = Matrix{Float64}(R) * Matrix{Float64}(R0)'
    return [Point3f(Vector{Float64}(P) + T * (Vector{Float64}(q) - Vector{Float64}(P0))) for q in corners]
end
_moved_corners(corners, ::Nothing, ::Nothing) = corners

"""
    _card_bbox(ctrl, obj)

The bounding box of the plots of `obj` in the controls `ctrl`, see `_selection_bbox`; for a system,
the union of the boxes of its rendered objects (see `_leaves`), or a box around its first object if
all of them are clipped.
"""
_card_bbox(ctrl::KinematicController, obj) = _selection_bbox(ctrl, obj, _object_plots(ctrl.h, obj))
function _card_bbox(ctrl::KinematicController, sys::BMO.AbstractSystem)
    leaves = _leaves(sys)
    isempty(leaves) && return GeometryBasics.Rect3d(fill(-5e-3, 3), fill(1e-2, 3))
    return _selection_bbox(ctrl, first(leaves), _object_plots(ctrl.h, sys))
end

"""Collapses the card `c` of the `gui` to its head, or expands it again."""
function _toggle_collapsed!(gui::LiveView, c::_ComponentCard)
    c.collapsed = !c.collapsed
    _show_head!(c)
    _update_cards!(gui)
    return nothing
end

"""
    _toggle_pinned!(gui, c)

Pins the card `c` of the selection to the selected object, which keeps the card next to the object
(or at the spot to which the mouse moved it, see `_drag_cards!`) independent of the selection; the
selection gets another card. Unpins a pinned card, which hides it and forgets its spot.
"""
function _toggle_pinned!(gui::LiveView, c::_ComponentCard)
    # The card of an inspected point or a measurement stays, see `_keep_info!`
    c.transient && return _keep_info!(gui, c)
    if c.pinned
        obj = c.obj
        c.pinned, c.obj, c.spot = false, nothing, nothing
        _hide_card!(c)
        # the plots of an inspected point or a measurement, unless it is still pinned elsewhere,
        # e.g. after it was docked, see `_forget!`
        _is_pinned(gui, obj) || _forget!(gui, obj)
    elseif c === gui.cards.selection && !isnothing(_shown_object(gui))
        c.pinned, c.obj, c.key = true, _shown_object(gui), nothing
        _use_card!(gui, _spare_card!(gui))
    end
    # Also resets the pin toggle of a card that cannot be pinned
    _show_head!(c)
    _update_cards!(gui)
    _on_pinned!(gui)
    return nothing
end

"""Unpins the cards of the `gui` that are pinned to `obj`, e.g. a removed clip plane."""
_unpin!(gui::LiveView, obj) = foreach(c -> _toggle_pinned!(gui, c), _floating_cards(gui, obj))

"""Returns `true` if a card of the `gui` is pinned to `obj`: the pinned cards float by default."""
_is_pinned(gui::LiveView, obj) = _is_floating(gui, obj)

"""Returns the pinned cards of the `gui` that float next to `obj` in the 3D view."""
_floating_cards(gui::LiveView, obj) = filter(c -> c.pinned && c.obj === obj, gui.cards.all)

"""Returns `true` if a pinned card of the `gui` floats next to `obj` in the 3D view."""
_is_floating(gui::LiveView, obj) = any(c -> c.pinned && c.obj === obj, gui.cards.all)

"""
    _float!(gui, obj)
    _dock!(gui, obj)

Moves the card pinned to `obj` from the sidebar into the 3D view, next to `obj`, or back. Only a
layout that docks the pinned cards, e.g. the app layout, implements them; by default, e.g. in the
compact layout, which has no sidebar, the pinned cards always float and these do nothing.
"""
_float!(::LiveView, _) = nothing
_dock!(::LiveView, _) = nothing

"""
    _toggle_pin!(gui, obj)

Pins a card to `obj`, independent of the card of the selection, or unpins the cards pinned to
`obj`. Where a pinned card is shown depends on the layout (see `_pin!`, `_unpin!` and
`_is_pinned`): floating next to its object by default, docked below the inspector in the app layout.
"""
function _toggle_pin!(gui::LiveView, obj)
    _is_pinned(gui, obj) ? _unpin!(gui, obj) : _pin!(gui, obj)
    _on_pinned!(gui)
    return nothing
end
_toggle_pin!(::LiveView, ::Nothing) = nothing

"""
    _pin!(gui, obj)
    _pin!(gui, c::_ComponentCard, obj)

Pins a floating card to `obj`: the card `c`, by default a spare card (see `_spare_card!`), see
`_toggle_pin!`. The card is placed next to `obj`, also if the mouse moved it before.
"""
_pin!(gui::LiveView, obj) =_pin!(gui, _spare_card!(gui), obj)
function _pin!(gui::LiveView, c::_ComponentCard, obj)
    c.pinned, c.obj, c.key, c.pose, c.spot = true, obj, nothing, nothing, nothing
    _show_head!(c)
    _update_cards!(gui)
    return nothing
end

# Called after a card was pinned or unpinned, see `AbstractLiveLayout`
_on_pinned!(::LiveView) = nothing

"""Returns a card of the `gui` that is neither pinned nor the card of the selection, or a new one."""
function _spare_card!(gui::LiveView)
    i = findfirst(c -> !c.pinned && c !== gui.cards.selection, gui.cards.all)
    isnothing(i) || return gui.cards.all[i]
    c = _ComponentCard(gui.fig, gui.layout, _card_z(length(gui.cards.all) + 1))
    push!(gui.cards.all, c)
    _connect_card!(gui, c)
    # The listeners of the new widgets come after the mouse shield of the cards, which must come last
    _shield_cards!(gui)
    return c
end

"""
    _ComponentCard(fig, layout::AbstractLiveLayout, z = _CARD_Z)

A floating card in the theme of the `layout`, with the tools of the `layout`, see `_card_tools!`.
"""
_ComponentCard(fig::Figure, layout::AbstractLiveLayout, z::Real = _CARD_Z) =
    _card_tools!(layout, _ComponentCard(fig, layout.theme, z))

"""
    _card_tools!(layout, c::_ComponentCard) -> c

Adds the tools of the `layout` to the head of the floating card `c`: none by default, the
`dock_button` in a layout that docks pinned cards (see `_dock!`).
"""
_card_tools!(::AbstractLiveLayout, c::_ComponentCard) = c

"""Makes `c` the card of the selection of the `gui`, whose `step_box` sets the keyboard step."""
function _use_card!(gui::LiveView, c::_ComponentCard)
    gui.cards.selection = c
    gui.widgets.step_box = c.step_box
    c.key, c.pose, c.spot = nothing, nothing, nothing
    return nothing
end

"""
    _connect_card!(gui, c)

Connects the head and the part of the selection (step, mode and the disclosure of the properties,
see `_connect_selection_part!`) of the card `c` of the `gui`; the declared widgets are connected
when they are built, see `_build_content!`.
"""
function _connect_card!(gui::LiveView, c::_ComponentCard)
    listeners = gui.controls.listeners
    push!(listeners, on(_ -> _toggle_collapsed!(gui, c), c.collapse_button.clicks))
    # The toggle switches itself, `_toggle_pinned!` sets it to the state of the card
    push!(listeners, on(v -> v == _pin_state(c) || _toggle_pinned!(gui, c), c.pin_button.active))
    push!(listeners, on(_ -> _toggle_properties!(gui, c), c.properties_button.clicks))
    _connect_selection_part!(gui, c)
    _connect_dock_button!(gui, c, c.dock_button)
    return nothing
end

#=
Part of the selection, see `_selection_part!`: shared by the floating card of the selection and
the inspector of the app layout, which both have the fields `step_box`, `mode` and `list`
=#

"""
    _connect_selection_part!(gui, part)

Connects the part of the selection `part` (the floating card of the selection or the inspector of
the app layout, see `_selection_part!`) to the `gui`: its step box sets the keyboard step (see
`_set_step!`) and takes the keyboard while focused, its mode control follows the mode of the
controls in both directions, see `_bind_mode!`.
"""
function _connect_selection_part!(gui::LiveView, part)
    listeners = gui.controls.listeners
    push!(listeners, on(s -> _set_step!(gui, s), part.step_box.stored_string))
    push!(listeners, on(_ -> _keep_keyboard!(gui), part.step_box.focused))
    _bind_mode!(gui, part.mode)
    return nothing
end

"""
    _bind_mode!(gui, mode::_Segmented)

Binds the Move/Rotate control `mode` to the mode of the controls of the `gui` in both directions: a
click sets the mode like the key `m` (see `_set_mode!`), the key `m` and a step in another unit (see
`_set_step!`) select the other option.
"""
function _bind_mode!(gui::LiveView, mode::_Segmented)
    ctrl = gui.controls
    sel = mode.selected
    push!(ctrl.listeners, on(m -> _set_mode!(gui, m), sel))
    push!(ctrl.listeners, on(m -> (sel[] == m || (sel[] = m)), ctrl.mode; update = true))
    return nothing
end

"""Sets the mode of the controls of the `gui`, like the key `m`."""
function _set_mode!(gui::LiveView, mode::Symbol)
    ctrl = gui.controls
    ctrl.mode[] == mode && return nothing
    ctrl.mode[] = mode
    _update_selection_box!(ctrl)
    _update_help!(ctrl)
    gui.status.text[] = "$mode mode, step: $(_step_string(mode, ctrl.fine_step, ctrl.fine_angle))"
    return nothing
end

# Names of `properties` that the part of the selection does not list: the header of the card and
# its pose rows show them
const _INSPECTOR_SKIPPED = ("Type", "Position [m]")

"""
    _show_properties!(gui, part, obj)

Shows the properties of `obj` (see `_inspector_rows`) in the property list of the part of the
selection `part`, see `_selection_part!`.
"""
_show_properties!(gui::LiveView, part, obj) = (_set_rows!(part.list, _inspector_rows(gui, obj)); nothing)

"""Returns the rows of the property list for `obj`, see `BeamletOptics.properties`."""
function _inspector_rows(::LiveView, obj)
    props = try
        BMO.properties(obj)
    catch e
        # a failing `properties` method of a user type must not break the live view
        Pair{String, Any}["Error" => sprint(showerror, e)]
    end
    return [_property_row(name, value) for (name, value) in props if !(name in _INSPECTOR_SKIPPED)]
end

"""
The rows of the property list of an inspected system (see `_inspect!`): its properties (see
`BeamletOptics.properties`), then the number of its objects and of its sources.
"""
function _inspector_rows(gui::LiveView, sys::BMO.AbstractSystem)
    rows = invoke(_inspector_rows, Tuple{LiveView, Any}, gui, sys)
    sources = count(p -> p.first === sys, gui.pairs)
    return Tuple{String, String}[rows; ("Objects", string(length(_leaves(sys)))); ("Sources", string(sources))]
end

"""Summary of the live view, shown without a selection, e.g. in the inspector of the app layout."""
function _inspector_rows(gui::LiveView, ::Nothing)
    objects = sum(h -> length(h.handles), gui.system_handles; init = 0)
    rows = Tuple{String, String}[
        ("Systems", string(length(gui.system_handles))), ("Objects", string(objects)),
        ("Sources", string(length(_sources(gui)))), ("Detector panels", string(length(gui.panels))),
        ("Clip planes", string(length(gui.clip.planes))),
        ("Last trace", gui.trace.solve_time > 0 ? _ms_string(gui.trace.solve_time) : "–")]
    return rows
end

"""
    _refresh_selection_part!(gui, c, obj)

Shows the properties of `obj` on the floating card `c` of the selection while they are expanded,
see `_toggle_properties!`; pinned cards have no part of the selection. Nothing for other hosts, e.g.
the inspector of the app layout refreshes its part itself, see `_refresh_inspector!`.
"""
_refresh_selection_part!(::LiveView, ::_AbstractCard, _) = nothing
function _refresh_selection_part!(gui::LiveView, c::_ComponentCard, obj)
    (c.pinned || !c.properties_shown) && return nothing
    _show_properties!(gui, c, obj)
    return nothing
end
_refresh_selection_part!(::LiveView, ::_ComponentCard, ::Nothing) = nothing

"""
Expands the properties of the floating card `c` of the selection of the `gui` below its disclosure
row, or collapses them again, see `_ComponentCard`. The state stays with the card when another
object is selected.
"""
function _toggle_properties!(gui::LiveView, c::_ComponentCard)
    c.properties_shown = !c.properties_shown
    _show_properties_state!(c.properties_button, c.properties_shown)
    _refresh_selection_part!(gui, c, _card_object(gui, c))
    _update_cards!(gui)
    return nothing
end

# The button that docks the pinned card `c` in the sidebar, if the layout has one, see `_card_tools!`
_connect_dock_button!(::LiveView, ::_ComponentCard, ::Nothing) = nothing
_connect_dock_button!(gui::LiveView, c::_ComponentCard, b::_IconButton) =
    (push!(gui.controls.listeners, on(_ -> _dock!(gui, c.obj), b.clicks)); nothing)

"""
    _shield_cards!(gui)

(Re)adds the listeners that keep the presses and the scrolling over the cards of the `gui` from the
camera: after the widgets of all cards (Textbox 70, Button 1), whose presses they would take
otherwise, and before the camera (0).
"""
function _shield_cards!(gui::LiveView)
    ev = events(gui.ax.scene)
    listeners = gui.controls.listeners
    foreach(off, gui.cards.shield)
    filter!(l -> !any(s -> s === l, gui.cards.shield), listeners)
    over = () -> any(c -> _over_card(c, ev), gui.cards.all)
    gui.cards.shield = Any[on(event -> Consume(event.action == Mouse.press && over()), ev.mousebutton; priority = 1),
        on(_ -> Consume(over()), ev.scroll; priority = 1)]
    append!(listeners, gui.cards.shield)
    return nothing
end

"""
    _outside_view(gui) -> Bool

Returns `true` if the mouse is over a part of the layout of the `gui` whose presses are not clicks
into the 3D view, see `ignore_mouse` of the controls: the controls neither select nor deselect on
them. By default none, see `AbstractLiveLayout`.
"""
_outside_view(::LiveView) = false

"""
    _connect_cards!(gui)

Connects the cards of the `gui`: their widgets, their update every frame and the mouse. Presses on
a card reach its widgets only: the controls ignore them (`ignore_mouse`) and the camera does not
get them, see `_shield_cards!`; a drag at its head moves the card, see `_drag_cards!`. A press
elsewhere ends the input into the textboxes of the cards, also if the controls consume it.
"""
function _connect_cards!(gui::LiveView)
    ctrl = gui.controls
    ev = events(gui.ax.scene)
    over = () -> any(c -> _over_card(c, ev), gui.cards.all)
    ctrl.ignore_mouse = () -> over() || _outside_view(gui)
    foreach(c -> _connect_card!(gui, c), gui.cards.all)
    push!(ctrl.listeners, on(_ -> _update_cards!(gui), ev.tick))
    # Before the controls (200)
    push!(ctrl.listeners, on(ev.mousebutton, priority = 250) do event
        event.action == Mouse.press && !over() && foreach(_defocus_card!, gui.cards.all)
        return Consume(false)
    end)
    _drag_cards!(gui)
    _shield_cards!(gui)
    _update_cards!(gui)
    return nothing
end

# Distance [px] that the mouse must move with the button down until a press at the head of a card
# moves it, and the longest pause [s] between the two clicks of a double click there
const _CARD_DRAG_MIN = 3.0f0
const _CARD_DOUBLE_CLICK = 0.4

"""
    _drag_cards!(gui)

Connects the mouse to the handles of the floating cards of the `gui`, i.e. their heads outside the
buttons (see `_over_handle`): a drag moves the card, which then stays at its new `spot` (see
`_card_spot`) with its line to the object, while the camera and the objects move, when it is pinned
and, on the card of the selection, when another object is selected. The cards without a spot keep
off it, see `_update_cards!`. A double click on the handle places the card next to its object again;
so do unpinning it and pinning another card, see `_pin!`. The listeners come before the widgets of
the cards (Button 1) and the mouse shield (1), after the controls (200), which ignore the presses on
the cards.
"""
function _drag_cards!(gui::LiveView)
    ev = events(gui.ax.scene)
    listeners = gui.controls.listeners
    # The card, the mouse at the press, the top left corner of the card at the press and whether it moved
    drag = Ref{Any}(nothing)
    # The card and the time of the last click on a handle, for the double click
    click = Ref{Any}((nothing, 0.0))
    push!(listeners, on(ev.mousebutton, priority = 2) do event
        event.button == Mouse.left || return Consume(false)
        if event.action == Mouse.press
            i = findlast(c -> _over_handle(c, ev), gui.cards.all)
            isnothing(i) && return Consume(false)
            c = gui.cards.all[i]
            r = c.background.layoutobservables.suggestedbbox[]
            drag[] = (c, Point2f(ev.mouseposition[]), Point2f(minimum(r)[1], maximum(r)[2]), false)
            return Consume(true)
        elseif event.action == Mouse.release && !isnothing(drag[])
            c, _, _, moved = drag[]
            drag[] = nothing
            moved && return Consume(true)
            prev, t = click[]
            if prev === c && time() - t < _CARD_DOUBLE_CLICK
                c.spot = nothing
                click[] = (nothing, 0.0)
                _update_cards!(gui)
            else
                click[] = (c, time())
            end
            return Consume(true)
        end
        return Consume(false)
    end)
    push!(listeners, on(ev.mouseposition, priority = 2) do xy
        isnothing(drag[]) && return Consume(false)
        c, start, p0, moved = drag[]
        d = Point2f(xy) - start
        (moved || maximum(abs, d) >= _CARD_DRAG_MIN) || return Consume(true)
        drag[] = (c, start, p0, true)
        size, view = _card_size(c), Rect2f(Makie.viewport(gui.ax.scene)[])
        # The spot of the card inside the view
        c.spot = _card_spot(_spot_position(_card_spot(p0 + d, size, view), size, view), size, view)
        _update_cards!(gui)
        return Consume(true)
    end)
    return nothing
end

