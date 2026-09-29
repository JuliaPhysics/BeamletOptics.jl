#=
Clip planes of the live view
=#

# Makie ignores all clip planes of a plot beyond the 8th
const _MAX_CLIP_PLANES = 8

const _LIVE_VIEW_HELP = "p: add clip plane, del: remove, c: clipping on/off, shift+c: flip\n1: show/hide sources"

"""Validates the `clip_planes` kwarg of `live_view` and returns a vector of `point => normal`."""
function _clip_plane_specs(clip_planes)
    specs = collect(clip_planes)
    length(specs) > _MAX_CLIP_PLANES &&
        throw(ArgumentError("at most $_MAX_CLIP_PLANES clip planes are supported by Makie"))
    for s in specs
        s isa Pair || throw(ArgumentError("invalid clip plane $s, use point => normal"))
        norm(s.second) < 1e-12 && throw(ArgumentError("the normal of a clip plane must not be zero, got $(s.second)"))
    end
    return specs
end

"""
    _apply_clip_planes!(gui::LiveView)

Applies the clip planes of the `gui` (none if `clipping` is off) to the scene and to all its plots
that did not set `clip_planes` explicitly, i.e. except the markers and the controls. The beams are
clipped only if `clip_beams` is set, which can be switched at runtime. Nothing is written as long
as there are no planes to apply or reset.
"""
function _apply_clip_planes!(gui::LiveView)
    scene = gui.ax.scene
    planes = gui.clip.enabled ? Plane3f[_plane3f(p) for p in gui.clip.planes] : Plane3f[]
    # Plots added later inherit the planes of the scene when they are created
    isempty(planes) && isempty(scene.theme.clip_planes[]) && return nothing
    scene.theme.clip_planes[] = planes
    beam_plots = Base.IdSet{Any}(p for h in gui.beam_handles for p in _beam_plots(h))
    for plot in scene.plots
        (haskey(plot.kw, :clip_planes) || plot in beam_plots) && continue
        # Nested plots inherit the planes of their parent
        plot.clip_planes[] == planes || (plot.clip_planes = planes)
    end
    # The beams are created with explicit `clip_planes`, hence they are set here in both cases
    beam_planes = gui.clip.beams ? planes : Plane3f[]
    for plot in beam_plots
        plot.clip_planes[] == beam_planes || (plot.clip_planes = beam_planes)
    end
    return nothing
end

"""Switches the clipping of the beams of the `gui` on or off, see `clip_beams`."""
function _set_clip_beams!(gui::LiveView, on::Bool)
    gui.clip.beams = on
    gui.widgets.clip_beams_toggle.active[] == on || (gui.widgets.clip_beams_toggle.active[] = on)
    _apply_clip_planes!(gui)
    return nothing
end

"""Switches the clip planes of the `gui` on or off, see `clipping`."""
function _set_clipping!(gui::LiveView, on::Bool)
    gui.clip.enabled = on
    _apply_clip_planes!(gui)
    gui.status.text[] = on ? "clipping on" : "clipping off"
    _on_clipping!(gui)
    return nothing
end

"""Re-applies the clip planes after the clip `plane` has been moved, without solving the systems."""
function _on_clip_change!(gui::LiveView, plane::LiveClipPlane)
    _apply_clip_planes!(gui)
    gui.status.text[] = _pose_string(gui, plane)
    return nothing
end

"""
    _add_clip_plane!(gui, point, normal; select = true)

Adds a clip plane through `point` with the `normal` to the `gui`: renders its marker, registers
it as a movable object of the controls and applies the planes. Selects the plane if `select`.
"""
function _add_clip_plane!(gui::LiveView, point, normal; select::Bool = true)
    ctrl = gui.controls
    plane = LiveClipPlane(point, normal, gui.clip.size)
    push!(ctrl.h.handles, _live_render_clip_plane!(gui.ax, plane;
        color = _clip_plane_color(gui.layout), strokecolor = _marker_stroke(gui.layout)))
    push!(ctrl.movable, plane)
    ctrl.init_poses[plane] = _pose(plane)
    push!(gui.clip.planes, plane)
    haskey(gui.labels, plane) || (gui.labels[plane] = _clip_plane_label(gui))
    _apply_clip_planes!(gui)
    _on_clip_planes_changed!(gui)
    if select
        ctrl.selected[] = plane
        _update_selection_box!(ctrl)
        gui.status.text[] = _pose_string(gui, plane)
    end
    return plane
end

"""
    _remove_clip_plane!(gui, plane)

Removes the clip `plane` from the `gui`: deletes its marker, forgets it in the controls, including
its entries of the undo history, deselects it and applies the remaining planes.
"""
function _remove_clip_plane!(gui::LiveView, plane::LiveClipPlane)
    ctrl = gui.controls
    i = findfirst(oh -> oh.obj === plane, ctrl.h.handles)
    if !isnothing(i)
        remove_render!(ctrl.h.handles[i])
        deleteat!(ctrl.h.handles, i)
    end
    filter!(o -> o !== plane, ctrl.movable)
    delete!(ctrl.init_poses, plane)
    delete!(ctrl.constraints, plane)
    filter!(e -> e.obj !== plane, ctrl.undo_stack)
    filter!(e -> e.obj !== plane, ctrl.redo_stack)
    !isnothing(ctrl.last_key_step) && ctrl.last_key_step.obj === plane && (ctrl.last_key_step = nothing)
    if ctrl.selected[] === plane
        ctrl.dragging = false
        ctrl.drag_start = nothing
        ctrl.press_kind = :none
        ctrl.press_leaf = nothing
        ctrl.selected[] = nothing
        _update_selection_box!(ctrl)
    end
    filter!(p -> p !== plane, gui.clip.planes)
    delete!(gui.labels, plane)
    _apply_clip_planes!(gui)
    _unpin!(gui, plane)
    gui.status.text[] = "clip plane removed"
    _on_clip_planes_changed!(gui)
    return nothing
end

# Other selected objects, e.g. after a click on a button of the card that is being hidden
_remove_clip_plane!(::LiveView, _) = nothing
_flip_clip_plane!(::LiveView, _) = nothing

"""Rotates the clip `plane` by π about its local x-axis, such that the other side is visible."""
function _flip_clip_plane!(gui::LiveView, plane::LiveClipPlane)
    ctrl = gui.controls
    P0, R0 = _pose(plane)
    _change!(() -> rotate3d!(plane, plane.dir[:, 1], π), ctrl, plane)
    P1, R1 = _pose(plane)
    ctrl.last_key_step = nothing
    _push_history!(ctrl, plane, P0, R0, P1, R1)
    update_render!(ctrl.h)
    _update_selection_box!(ctrl)
    _on_clip_change!(gui, plane)
    return nothing
end

"""
    _clip_key!(gui::LiveView, key)

Handles the clip plane keys of the `gui`: `p` adds a plane through the selected object (or the
camera `lookat`) along the view direction, `Delete` removes the selected plane, `c` switches
clipping on and off and `Shift+c` flips the selected plane. Returns whether the key was handled.
"""
function _clip_key!(gui::LiveView, key)
    ctrl = gui.controls
    scene = gui.ax.scene
    sel = ctrl.selected[]
    if key == Keyboard.c
        if _shift_pressed(scene)
            sel isa LiveClipPlane || return false
            _flip_clip_plane!(gui, sel)
            return true
        end
        isempty(gui.clip.planes) && return false
        _set_clipping!(gui, !gui.clip.enabled)
        return true
    elseif key == Keyboard.p
        # Nothing can be selected in the spectator mode
        ctrl.spectator[] && return false
        if length(gui.clip.planes) >= _MAX_CLIP_PLANES
            gui.status.text[] = "at most $_MAX_CLIP_PLANES clip planes"
            return true
        end
        cam = cameracontrols(scene)
        lookat, eye = Vector{Float64}(cam.lookat[]), Vector{Float64}(cam.eyeposition[])
        point = isnothing(sel) ? lookat : Vector{Float64}(position(sel))
        _add_clip_plane!(gui, point, lookat - eye)
        return true
    elseif key == Keyboard.delete
        sel isa LiveClipPlane || return false
        _remove_clip_plane!(gui, sel)
        return true
    end
    return false
end

"""Connects the clip plane keys of the `gui`, see `_clip_key!`."""
function _connect_clip_planes!(gui::LiveView)
    push!(gui.controls.listeners, on(events(gui.ax.scene).keyboardbutton, priority = 200) do event
        event.action == Keyboard.press || return Consume(false)
        gui.controls.ignore_keys() && return Consume(false)
        return Consume(_clip_key!(gui, event.key))
    end)
    gui.controls.help_extra = _LIVE_VIEW_HELP
    _update_help!(gui.controls)
    return nothing
end

