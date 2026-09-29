#=
Selection of the live view: labels, the component menu, hiding and showing, the source markers
=#

"""Formats the position of `obj` in mm."""
function _position_string(obj)
    p = round.(1e3 .* collect(Float64, position(obj)), digits = 6)
    return "(" * join(p, ", ") * ") mm"
end

"""
    _unit_string(x, units)

Formats `x` with 3 significant digits in the first of the `units` (`factor => unit`, from small to
large) in which it is below 1000, e.g. 1000 µm as 1 mm.
"""
function _unit_string(x, units)
    for (i, (factor, unit)) in enumerate(units)
        v = round(x / factor, sigdigits = 3)
        (v < 1000 || i == length(units)) && return "$(_fmt_sigdigits(v)) $unit"
    end
end

_length_string(x) = _unit_string(x, (1e-9 => "nm", 1e-6 => "µm", 1e-3 => "mm"))
_angle_string(x) = _unit_string(x, (1e-6 => "µrad", 1e-3 => "mrad", deg2rad(1) => "°"))

"""
    _label(gui, obj)

The name of `obj` in the `gui`: its entry in `labels`, else its automatic name, e.g. "Lens 2", see
`_name_objects!`, else the name of its type. The same in all layouts.
"""
_label(gui::LiveView, obj) = get(() -> get(gui.objects.names, obj, string(nameof(typeof(obj)))),
    gui.labels, obj)

"""Describes the pose of `obj`, including the change since the controls were enabled."""
function _pose_string(gui, obj)
    s = "$(_label(gui, obj)) at $(_position_string(obj))"
    haskey(gui.controls.init_poses, obj) || return s
    P0, R0 = gui.controls.init_poses[obj]
    P, R = _pose(obj)
    _, angle = _axis_angle_from_rotmatrix(R * R0')
    return s * ", moved by $(_length_string(norm(P - P0))), rotated by $(_angle_string(angle))"
end

"""
    _menu_entries(ctrl)

Returns `(obj, depth)` of all movable objects of the controls `ctrl` in the order of the component
menu: the top-level objects in the order of the handles, the objects of a group after the group,
with the nesting `depth` of the group. Clip planes are not listed.
"""
function _menu_entries(ctrl::KinematicController)
    entries = Tuple{Any, Int}[]
    function add!(obj, depth)
        push!(entries, (obj, depth))
        foreach(c -> add!(c, depth + 1), _children(obj))
        return nothing
    end
    for obj in ctrl.movable
        obj isa LiveClipPlane || add!(obj, 0)
    end
    return entries
end

"""
    _menu_entries(gui)

Returns `(key, depth)` of the entries of the component menu of the `gui`: per system its
`SystemRenderHandle` (depth 0), which inspects the system (see `_inspect!`), followed by its movable
objects (see `_menu_entries(ctrl)`) one level deeper, then the other movable objects, e.g. the
sources and the extras.
"""
function _menu_entries(gui::LiveView)
    entries = _menu_entries(gui.controls)
    out = Tuple{Any, Int}[]
    placed = falses(length(entries))
    for h in gui.system_handles
        push!(out, (h, 0))
        tops = Base.IdSet{Any}(_top_levels(h))
        inside = false
        for (i, (obj, depth)) in enumerate(entries)
            # the objects of a group follow the group
            depth == 0 && (inside = obj in tops)
            if inside && !placed[i]
                push!(out, (obj, depth + 1))
                placed[i] = true
            end
        end
    end
    append!(out, (e for (i, e) in enumerate(entries) if !placed[i]))
    return out
end

#=
Component menu, hide and show
=#

"""
    _menu_options(labels, entries)
    _menu_options(gui, entries)

Returns the options of the component menu of the `entries`, see `_menu_entries`, named after the
`labels` (before the `gui` exists) or after `_label` of the `gui`.
"""
function _menu_options(labels, entries)
    isempty(entries) && return [("no components", 0)]
    return [("  "^depth * get(labels, obj, string(nameof(typeof(obj)))), i)
            for (i, (obj, depth)) in enumerate(entries)]
end
function _menu_options(gui::LiveView, entries)
    isempty(entries) && return [("no components", 0)]
    return [("  "^depth * _label(gui, obj), i) for (i, (obj, depth)) in enumerate(entries)]
end

"""
Shows the entries of `_menu_entries(gui)` with the automatic names of `_label` in the component
menu of the `gui`, if it has one: the systems and their objects, see `gui.objects.menu`.
"""
_refresh_menu_options!(::LiveView, ::Nothing) = nothing
function _refresh_menu_options!(gui::LiveView, menu::Menu)
    entries = _menu_entries(gui)
    gui.objects.menu = Any[first.(entries)...]
    options = _menu_options(gui, entries)
    menu.options[] == options || (menu.options[] = options)
    return nothing
end

"""
    _on_menu_select!(gui, i)

Selects the object of the option `i` of the component menu, like a click in the 3D view, or
inspects the system of the option, see `_select!`. The option `0`, i.e. no selection, is set by
`_on_select!` and ignored.
"""
function _on_menu_select!(gui::LiveView, i)
    1 <= i <= length(gui.objects.menu) || return nothing
    _select!(gui, gui.objects.menu[i])
    return nothing
end

"""
    _select!(gui, obj)
    _select!(gui, h::SystemRenderHandle)

Selects the movable `obj` like a click in the 3D view, e.g. from the component menu or the object
tree, and shows its pose in the status line. Nothing is selected in the spectator mode. The entry
`h` of a system inspects the system instead, see `_inspect!`.
"""
_select!(gui::LiveView, h::SystemRenderHandle) = _inspect!(gui, h)
function _select!(gui::LiveView, obj)
    ctrl = gui.controls
    ctrl.selected[] === obj && return nothing
    if ctrl.spectator[]
        gui.status.text[] = "spectator mode, press v to select components"
        _show_menu_selection!(gui.widgets.menu, 0)
        return nothing
    end
    ctrl.selected[] = obj
    _update_selection_box!(ctrl)
    gui.status.text[] = _pose_string(gui, obj)
    return nothing
end

"""
Shows the selected object of the controls in the component menu and in the component card. Any
change of the selection ends the inspection, see `_inspect!`.
"""
function _on_select!(gui::LiveView)
    gui.objects.inspected = nothing
    _on_shown!(gui)
    return nothing
end

"""
    _on_shown!(gui)

Shows the object of the card of the selection (see `_shown_object`) in the component menu, on the
cards, in the inspector and in the object tree, after the selection or the inspection changed.
"""
function _on_shown!(gui::LiveView)
    obj = _shown_object(gui)
    key = _row_key(gui, obj)
    i = isnothing(obj) ? nothing : findfirst(o -> o === key, gui.objects.menu)
    _show_menu_selection!(gui.widgets.menu, something(i, 0))
    _update_inspector!(gui)
    _update_cards!(gui)
    _on_selected!(gui)
    _show_inspected!(gui, gui.objects.inspected)
    return nothing
end

#=
Inspection: an object that is shown on the card of the selection without being selected for moving
=#

"""
    _shown_object(gui)

The object that the card of the selection of the `gui` shows: the selected object of the controls,
else the inspected one (see `_inspect!`), else `nothing`.
"""
function _shown_object(gui::LiveView)
    sel = gui.controls.selected[]
    return isnothing(sel) ? gui.objects.inspected : sel
end

"""
    _inspect!(gui, obj)
    _inspect!(gui, h::SystemRenderHandle)

Shows `obj` on the card of the selection of the `gui` (floating next to it, or in the inspector of
the app layout) without selecting it for moving, i.e. without gizmo and selection box: a system
(also given by its handle `h`, e.g. from the component menu or the object tree) or an object that
is not movable. The selection of the controls is cleared, since `gui.objects.inspected` and
`controls.selected[]` exclude each other. The inspection ends with Esc, a click on the empty space
of the 3D view or a new selection, see `_end_inspection!` and `_on_select!`. The pose boxes of an
inspected object that is not movable reject inputs, see `_apply_pose_input!`.
"""
_inspect!(gui::LiveView, h::SystemRenderHandle) = _inspect!(gui, h.sys)
function _inspect!(gui::LiveView, obj)
    gui.objects.inspected === obj && return nothing
    ctrl = gui.controls
    if !isnothing(ctrl.selected[])
        ctrl.selected[] = nothing
        _update_selection_box!(ctrl)
    end
    gui.objects.inspected = obj
    # the bounding box of the new object, see `_card_corners`
    gui.cards.selection.key = nothing
    gui.status.text[] = "$(_label(gui, obj)) inspected, not selected for moving, Esc closes its card"
    _on_shown!(gui)
    return nothing
end

"""Ends the inspection of the `gui`, see `_inspect!`; nothing without one."""
function _end_inspection!(gui::LiveView)
    isnothing(gui.objects.inspected) && return nothing
    gui.objects.inspected = nothing
    _on_shown!(gui)
    return nothing
end

"""
    _row_key(gui, obj)

The key of `obj` in the component menu and in the object tree: the object itself, the
`SystemRenderHandle` of a system (or of the extras).
"""
_row_key(::LiveView, obj) = obj
function _row_key(gui::LiveView, sys::BMO.AbstractSystem)
    for h in (gui.system_handles..., gui.extras)
        h.sys === sys && return h
    end
    return sys
end

"""
Shows the inspected object of the `gui` besides the cards, e.g. its row in the object tree of the app
layout; nothing by default.
"""
_show_inspected!(::LiveView, _) = nothing

_show_menu_selection!(::Nothing, _) = nothing
function _show_menu_selection!(menu::Menu, i)
    menu.i_selected[] == i || (menu.i_selected[] = i)
    return nothing
end

#=
Optional widgets: a layout may leave out the component menu and the hide buttons (`nothing`), their
observables are then `nothing` and not listened to
=#

_clicks(::Nothing) = nothing
_clicks(button) = button.clicks
_menu_selection(::Nothing) = nothing
_menu_selection(menu::Menu) = menu.i_selected
_open_observable(::Nothing) = nothing
_open_observable(menu::Menu) = menu.is_open
_is_open(::Nothing) = false
_is_open(menu::Menu) = menu.is_open[]

"""Adds a listener `f` of the observable `obs` to the `listeners`, unless `obs` is `nothing`."""
_listen!(_, _, ::Nothing) = nothing
function _listen!(listeners, f, obs::Observable)
    push!(listeners, on(f, obs))
    return nothing
end

"""
Sets the `visible` attribute of all plots of the rendered objects (leaves) of `obj`. A leaf shown
again whose opacity was set to 0 gets its initial opacity back, see `_set_opacity!`.
"""
function _set_hidden!(gui::LiveView, obj, hide::Bool)
    for leaf in _leaves(obj)
        hide ? push!(gui.objects.hidden, leaf) : delete!(gui.objects.hidden, leaf)
        i = findfirst(oh -> oh.obj === leaf, gui.controls.h.handles)
        isnothing(i) && continue
        visible = !hide && (gui.widgets.sources_toggle.active[] || !_is_source(leaf))
        for plot in gui.controls.h.handles[i].plots
            plot.visible[] == visible || (plot.visible[] = visible)
        end
        hide || _restore_opacity!(gui, leaf, get(gui.objects.opacity, leaf, nothing))
    end
    return nothing
end

_is_source(obj) = obj isa Union{BMO.AbstractBeam, BMO.AbstractBeamGroup}

"""
    _set_show_sources!(gui, show)

Shows or hides the markers of all movable sources of the `gui`, e.g. if a marker covers small
components. Hidden markers can not be selected in the 3D view, a selected source is deselected.
Sources hidden via the "hide" button stay hidden.
"""
function _set_show_sources!(gui::LiveView, show::Bool)
    ctrl = gui.controls
    for oh in ctrl.h.handles
        _is_source(oh.obj) || continue
        visible = show && !(oh.obj in gui.objects.hidden)
        for plot in oh.plots
            plot.visible[] == visible || (plot.visible[] = visible)
        end
    end
    if !show && _is_source(ctrl.selected[])
        ctrl.selected[] = nothing
        _update_selection_box!(ctrl)
    end
    gui.status.text[] = show ? "sources shown" : "sources hidden, press 1 to show them"
    return nothing
end

"""
Connects the "sources" toggle and the key `1` of the `gui`, see `_set_show_sources!`. Not `s`,
which moves the camera backwards (WASD keys of Makie's `Camera3D`).
"""
function _connect_sources!(gui::LiveView)
    listeners = gui.controls.listeners
    push!(listeners, on(v -> _set_show_sources!(gui, v), gui.widgets.sources_toggle.active))
    push!(listeners, on(events(gui.ax.scene).keyboardbutton, priority = 200) do event
        (event.action == Keyboard.press && event.key == Keyboard._1) || return Consume(false)
        gui.controls.ignore_keys() && return Consume(false)
        gui.widgets.sources_toggle.active[] = !gui.widgets.sources_toggle.active[]
        return Consume(true)
    end)
    gui.widgets.sources_toggle.active[] || _set_show_sources!(gui, false)
    return nothing
end

"""
    _toggle_hidden!(gui, obj)

Hides the object `obj` of the `gui` (of a card or a row of the object tree), i.e. makes the plots
of its rendered objects invisible and clears the selection if it is hidden now, e.g. `obj` or an
object of a hidden group or system, or shows it again if all of them are hidden. A hidden object
can not be selected in the 3D view, but it stays in the systems.
"""
function _toggle_hidden!(gui::LiveView, obj)
    ctrl = gui.controls
    hide = !_all_hidden(gui, obj)
    _set_hidden!(gui, obj, hide)
    if hide
        sel = ctrl.selected[]
        if sel === obj || (!isnothing(sel) && _all_hidden(gui, sel))
            ctrl.selected[] = nothing
            _update_selection_box!(ctrl)
        end
        gui.status.text[] = "$(_label(gui, obj)) hidden, $(_show_hint(gui))"
    else
        gui.status.text[] = "$(_label(gui, obj)) shown"
    end
    _update_inspector!(gui)
    _update_cards!(gui)
    _on_hidden!(gui)
    return nothing
end
_toggle_hidden!(gui::LiveView, ::Nothing) = (gui.status.text[] = "select a component to hide it"; nothing)
function _toggle_hidden!(gui::LiveView, ::LiveClipPlane)
    gui.status.text[] = "clip planes can not be hidden, press c to switch clipping off"
    return nothing
end

"""Returns `true` if all rendered objects (leaves) of `obj` are hidden in the `gui`."""
_all_hidden(gui::LiveView, obj) = all(leaf -> leaf in gui.objects.hidden, _leaves(obj))

"""Returns how a hidden object of the `gui` is shown again, for the status line."""
_show_hint(::LiveView) = "select it in the menu to show it again"

"""Shows all hidden objects of the `gui`."""
function _show_all!(gui::LiveView)
    foreach(leaf -> _set_hidden!(gui, leaf, false), collect(gui.objects.hidden))
    gui.status.text[] = "all components shown"
    _update_inspector!(gui)
    _update_cards!(gui)
    _on_hidden!(gui)
    return nothing
end


"""
Connects the export button, the component menu, "show all" and Esc, which ends the inspection (see
`_inspect!`); the cards connect their widgets.
"""
function _connect_tools!(gui::LiveView)
    listeners = gui.controls.listeners
    push!(listeners, on(_ -> _export!(gui), gui.widgets.export_button.clicks))
    _listen!(listeners, i -> _on_menu_select!(gui, i), _menu_selection(gui.widgets.menu))
    push!(listeners, on(_ -> _on_select!(gui), gui.controls.selected))
    _listen!(listeners, _ -> _show_all!(gui), _clicks(gui.widgets.show_all_button))
    # Before the controls, which consume Esc to deselect; passed on
    push!(listeners, on(events(gui.ax.scene).keyboardbutton, priority = 201) do event
        (event.action == Keyboard.press && event.key == Keyboard.escape) || return Consume(false)
        gui.controls.ignore_keys() || _end_inspection!(gui)
        return Consume(false)
    end)
    return nothing
end

"""Returns `true` if the component menu (if any) or the views menu of the `gui` is open."""
_menu_open(gui::LiveView) = _is_open(gui.widgets.menu) || _is_open(gui.widgets.views_menu)

"""Clip planes are numbered, e.g. "Clip plane 2", see `_label`."""
function _clip_plane_label(gui::LiveView)
    n = get(gui.objects.counters, "Clip plane", 0) + 1
    gui.objects.counters["Clip plane"] = n
    return "Clip plane $n"
end
