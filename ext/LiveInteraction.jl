import Makie
using Makie: events, on, off, Consume, Mouse, Keyboard, mouseposition_px, Point2f, Vec3f
using LinearAlgebra: tr, diag, I, pinv

"""
    _ray_plane_intersect(origin, dir, plane_point, normal)

Returns the intersection of the ray `origin + t * dir` with the plane through `plane_point`, or
`nothing` if the ray is parallel to the plane or points away from it.
"""
function _ray_plane_intersect(origin, dir, plane_point, normal)
    denom = dot(dir, normal)
    abs(denom) < 1e-10 && return nothing
    t = dot(plane_point .- origin, normal) / denom
    t < 0 && return nothing
    return origin .+ t .* dir
end

"""
    _axis_angle_from_rotmatrix(R)

Returns the `axis` and `angle` of the rotation matrix `R` such that `rotate3d(axis, angle) ≈ R`.
"""
function _axis_angle_from_rotmatrix(R::AbstractMatrix{T}) where {T}
    angle = acos(clamp((tr(R) - 1) / 2, -one(T), one(T)))
    angle < 1e-8 && return [0.0, 0.0, 1.0], 0.0
    if angle > π - 1e-6
        # sin(angle) ≈ 0, use the symmetric part of R instead
        S = (R + I) / 2
        d = diag(S)
        i = argmax(d)
        axis = zeros(3)
        axis[i] = sqrt(max(d[i], 0.0))
        for j in 1:3
            j == i && continue
            axis[j] = S[i, j] / axis[i]
        end
        return normalize(axis), angle
    end
    axis = [R[3, 2] - R[2, 3], R[1, 3] - R[3, 1], R[2, 1] - R[1, 2]] ./ (2 * sin(angle))
    return normalize(axis), angle
end

"""Returns the 12 edges of the box `bb` as input for `linesegments`."""
function _bbox_wireframe(bb)
    lo, hi = minimum(bb), maximum(bb)
    c = [Point3f(x, y, z) for x in (lo[1], hi[1]), y in (lo[2], hi[2]), z in (lo[3], hi[3])]
    pts = Point3f[]
    for j in 1:2, k in 1:2
        push!(pts, c[1, j, k], c[2, j, k])
    end
    for i in 1:2, k in 1:2
        push!(pts, c[i, 1, k], c[i, 2, k])
    end
    for i in 1:2, j in 1:2
        push!(pts, c[i, j, 1], c[i, j, 2])
    end
    return pts
end

_shift_pressed(scene) = Keyboard.left_shift in events(scene).keyboardstate ||
                        Keyboard.right_shift in events(scene).keyboardstate

"""Returns whether the `modifier` (a `Keyboard.Button`, or a tuple/vector of them meaning "any of")
is currently held in the `scene`."""
_modifier_held(scene, modifier::Keyboard.Button) = modifier in events(scene).keyboardstate
_modifier_held(scene, modifiers) = any(m -> m in events(scene).keyboardstate, modifiers)

"""Returns a short display name for a `select_modifier` (a `Keyboard.Button`, or a tuple/vector of
them), e.g. `"left_shift"`."""
_modifier_name(btn::Keyboard.Button) = string(btn)
_modifier_name(btns) = join(string.(btns), "/")

# Ctrl (Windows/Linux) or Cmd (macOS) held, used for the undo/redo shortcuts
const _CTRL_OR_CMD_KEYS = (Keyboard.left_control, Keyboard.right_control,
    Keyboard.left_super, Keyboard.right_super)

_px(scene) = Tuple(Float64.(mouseposition_px(scene)))

_other_mode(mode::Symbol) = mode == :move ? :rotate : :move

"""Formats `x` with 3 significant digits, without a trailing `.0` for integer values."""
function _fmt_sigdigits(x)
    v = round(x, sigdigits = 3)
    return isinteger(v) && abs(v) < 1e15 ? string(Int(v)) : string(v)
end

"""Returns the current keyboard step of the `mode`, in nm for the move and µrad for the rotate mode."""
function _step_string(mode::Symbol, fine_step, fine_angle)
    return mode == :move ? "$(_fmt_sigdigits(fine_step * 1e9)) nm" :
           "$(_fmt_sigdigits(fine_angle * 1e6)) µrad"
end

const _SPECTATOR_HINT = "spectator mode, v: edit, h: show controls"

function _help_hint(mode::Symbol, fine_step, fine_angle)
    step = _step_string(mode, fine_step, fine_angle)
    return "$mode mode, step $step, +/-: step, m: switch mode, v: spectator, h: show controls"
end

const _SPECTATOR_HELP = """
    spectator mode, v: switch to edit mode
    all clicks and drags: camera
    components can not be selected or moved
    h: hide controls"""

function _help_text(mode::Symbol, fine_step, fine_angle, select_modifier = nothing)
    step = _step_string(mode, fine_step, fine_angle)
    if mode == :move
        verb = "move along"
        up, left, page, drag = "green arrow", "red arrow", "blue arrow", "move in plane"
    else
        verb = "rotate around"
        up, left, page, drag = "red ring", "blue ring", "green ring", "rotate around blue ring"
    end
    click = isnothing(select_modifier) ? "click" : "$(_modifier_name(select_modifier))+click"
    return """
    $mode mode, m: switch to $(_other_mode(mode)) mode
    $click: select, again: part of a group
    drag selection: $drag, other drags: camera
    ↑/↓: $verb $up
    ←/→: $verb $left
    page up/down: $verb $page
    step: $step, +/-: change step, shift: 10× step
    backspace: reset, esc: enclosing group or deselect
    ctrl/cmd+z: undo, ctrl/cmd+y or ctrl/cmd+shift+z: redo
    v: spectator mode, h: hide controls"""
end

const _STEP_MANTISSAS = (1.0, 2.0, 5.0)
# Bounds of the keyboard steps, changed via + and -
const _FINE_STEP_LIMITS = (1e-12, 1.0)
const _FINE_ANGLE_LIMITS = (1e-9, π / 4)

"""
    _next_step(x, dir::Int)

Returns the next value after `x > 0` in the 1-2-5 sequence (…, 1, 2, 5, 10, 20, …) in the direction
`dir` (`1`: up, `-1`: down). If `x` is not a value of the sequence, the next value of the sequence
in the direction `dir` is returned.
"""
function _next_step(x::Real, dir::Int)
    x > 0 || throw(ArgumentError("step must be positive, got $x"))
    dir in (-1, 1) || throw(ArgumentError("dir must be -1 or 1, got $dir"))
    e = floor(log10(x))
    # Adjacent decades as well, in case of rounding errors of log10
    grid = [k * 10.0^d for d in (e - 1):(e + 2) for k in _STEP_MANTISSAS]
    i = findfirst(g -> isapprox(x, g; rtol = 1e-6), grid)
    isnothing(i) || return grid[i + dir]
    return dir > 0 ? grid[findfirst(>(x), grid)] : grid[findlast(<(x), grid)]
end

# Green, red and blue axes of the controls: local y-axis, local x-axis and rotation axis
const _AXES_COLORS = [:green, :red, :blue]
# Symbols of the gizmo axes in the same order as _AXES_COLORS, and as accepted by `constraints`
const _AXES_SYMS = (:y, :x, :v)
const _GIZMO_AXES = (:x, :y, :v)
const _GIZMO_FADE_ALPHA = 0.15
const _RING_RES = 32

"""
    _gizmo(mode, origin, axes, l)

Returns the arrows, labels and rings that visualize the controls `axes` at the `origin`. In the
move mode, the arrows of length `l` point along the `axes`. In the rotate mode, the arrows mark the
positive direction of rotation on rings with radius `l` around the `axes`.
"""
function _gizmo(mode::Symbol, origin, axes, l)
    r = l
    arrow_pos = Point3f[]
    arrow_dir = Vec3f[]
    label_pos = Point3f[]
    ring_pts = Point3f[]
    for a in axes
        e1 = normalize(cross(a, abs(a[1]) < 0.9 ? [1, 0, 0] : [0, 1, 0]))
        e2 = cross(a, e1)
        ts = LinRange(0, 1.5π, _RING_RES + 1)
        for i in 1:_RING_RES
            push!(ring_pts, Point3f(origin + r * (cos(ts[i]) * e1 + sin(ts[i]) * e2)))
            push!(ring_pts, Point3f(origin + r * (cos(ts[i + 1]) * e1 + sin(ts[i + 1]) * e2)))
        end
        if mode == :move
            push!(arrow_pos, Point3f(origin))
            push!(arrow_dir, Vec3f(l * a))
            push!(label_pos, Point3f(origin + 1.2 * l * a))
        else
            # Arrow tip at the end of the ring, label on the axis of the ring
            t = last(ts)
            push!(arrow_pos, Point3f(origin + r * (cos(t) * e1 + sin(t) * e2)))
            push!(arrow_dir, Vec3f(0.35 * r * (-sin(t) * e1 + cos(t) * e2)))
            push!(label_pos, Point3f(origin + 1.3 * r * a))
        end
    end
    return arrow_pos, arrow_dir, label_pos, ring_pts
end

"""One entry of the undo history of [`kinematic_controls!`](@ref): the pose of `obj` before and
after one gesture (a mouse drag, a reset or one or several merged keyboard steps)."""
struct _HistoryEntry
    obj::_LiveMovable
    P0::Point3{Float64}
    R0::Matrix{Float64}
    P1::Point3{Float64}
    R1::Matrix{Float64}
end

# Last key step, see `_record_key_step!`
const _KeyStepInfo = NamedTuple{(:obj, :key, :time), Tuple{_LiveMovable, Keyboard.Button, Float64}}

"""
    KinematicController

Returned by [`kinematic_controls!`](@ref). The currently selected object is stored in the
`selected` `Observable`, the current mode (`:move` or `:rotate`) in the `mode` `Observable`. Use
`close` to remove the controls.
"""
mutable struct KinematicController{H <: SystemRenderHandle}
    ax::_RenderEnv
    h::H
    movable::Vector{_LiveMovable}
    init_poses::IdDict{_LiveMovable, Tuple{Point3{Float64}, Matrix{Float64}}}
    selected::Observable{Union{Nothing, _LiveMovable}}
    mode::Observable{Symbol}
    on_change::Function
    plane_normal::Vector{Float64}
    rotation_axis::Vector{Float64}
    rotate_speed::Float64
    fine_step::Float64
    fine_angle::Float64
    throttle::Bool
    # all input goes to the camera, toggled via v
    spectator::Observable{Bool}
    select_modifier::Any
    drag_threshold::Float64
    # locked move/rotate axes per object, see the `constraints` kwarg
    constraints::IdDict{Any, NamedTuple}
    # screen-space pick radius of sources [px]
    source_pick_radius::Float64
    # disables all keyboard handling while true, e.g. a focused text field elsewhere in the figure
    ignore_keys::Function
    # interaction state
    dirty::Bool
    dragging::Bool
    plane_point::Vector{Float64}
    grab_offset::Vector{Float64}
    last_mouse::NTuple{2, Float64}
    press_pos::Union{Nothing, NTuple{2, Float64}}
    press_leaf::Union{Nothing, _LiveMovable}
    press_kind::Symbol
    # pose of the selected object at the start of the current drag, for the undo history
    drag_start::Union{Nothing, Tuple{Point3{Float64}, Matrix{Float64}}}
    # undo/redo history, one entry per gesture, see `_push_history!` / `_record_key_step!`
    undo_stack::Vector{_HistoryEntry}
    redo_stack::Vector{_HistoryEntry}
    last_key_step::Union{Nothing, _KeyStepInfo}
    # selection box and gizmo of the keyboard controls
    box_obs::Observable{Vector{Point3f}}
    arrow_pos::Observable{Vector{Point3f}}
    arrow_dir::Observable{Vector{Vec3f}}
    label_pos::Observable{Vector{Point3f}}
    ring_pts::Observable{Vector{Point3f}}
    arrow_color::Observable{Vector{Makie.RGBAf}}
    label_color::Observable{Vector{Makie.RGBAf}}
    ring_color::Observable{Vector{Makie.RGBAf}}
    gizmo_size::Observable{Float64}
    gizmo_visible::Observable{Bool}
    # controls overlay, toggled via h
    help_obs::Observable{String}
    help_shown::Bool
    # additional lines of the overlay, e.g. the keys of `live_view`
    help_extra::String
    plots::Vector{AbstractPlot}
    listeners::Vector{Any}
    # last error of on_change, logged only once
    last_error::Union{Nothing, String}
end

function Base.show(io::IO, ctrl::KinematicController)
    obj = ctrl.selected[]
    sel = isnothing(obj) ? "none" : string(nameof(typeof(obj)))
    print(io, "KinematicController(", length(ctrl.movable), " movable, selected = ", sel,
        ", mode = ", ctrl.mode[], ")")
end

"""Returns the objects of the group `obj`, or an empty vector if `obj` is not a group."""
_children(obj) = obj isa BMO.AbstractObjectGroup ? collect(_LiveMovable, BMO.shape(obj)) :
                 _LiveMovable[]

"""Returns all objects of the group `obj` that are not groups themselves (recursively), or `[obj]`."""
function _leaves(obj)
    obj isa BMO.AbstractObjectGroup || return _LiveMovable[obj]
    return reduce(vcat, (_leaves(c) for c in BMO.shape(obj)); init = _LiveMovable[])
end

"""Returns the `obj` and all its nested objects and subgroups (recursively)."""
function _descendants(obj)
    out = _LiveMovable[obj]
    for c in _children(obj)
        append!(out, _descendants(c))
    end
    return out
end

"""Returns the chain `[leaf, parent of leaf, …, top-level object]` of the hierarchy of `ctrl.h`."""
function _chain(ctrl, leaf)
    chain = _LiveMovable[leaf]
    while haskey(ctrl.h.parent, last(chain))
        push!(chain, ctrl.h.parent[last(chain)])
    end
    return chain
end

"""An object is movable if its top-level object is one of the movable objects of the `ctrl`."""
function _is_movable(ctrl::KinematicController, obj)
    top = _top_level(ctrl.h, obj)
    return any(o -> o === top, ctrl.movable)
end

function _object_plots(h::SystemRenderHandle, obj)
    plots = AbstractPlot[]
    for leaf in _leaves(obj)
        i = findfirst(oh -> oh.obj === leaf, h.handles)
        isnothing(i) || append!(plots, _pickable_plots(h.handles[i]))
    end
    return plots
end

"""
    _drill_select(ctrl, leaf)

Returns the new selection after a click on the `leaf`: the top-level object of the `leaf` on the
first click, then one level further down the hierarchy towards the `leaf` on each further click.
"""
function _drill_select(ctrl::KinematicController, leaf)
    chain = _chain(ctrl, leaf)
    sel = ctrl.selected[]
    i = isnothing(sel) ? nothing : findfirst(o -> o === sel, chain)
    isnothing(i) && return last(chain)
    return chain[max(i - 1, 1)]
end

"""Axes of the keyboard controls: local y-axis, local x-axis and rotation axis of the `obj`."""
function _control_axes(ctrl::KinematicController, obj)
    R = _pose(obj)[2]
    return (Vector{Float64}(R[:, 2]), Vector{Float64}(R[:, 1]), ctrl.rotation_axis)
end

"""Returns the vectors of the gizmo axes `syms` (a subset of `:x`, `:y`, `:v`) of `obj`."""
function _axis_vectors(ctrl::KinematicController, obj, syms)
    y, x, v = _control_axes(ctrl, obj)
    lookup = (y = y, x = x, v = v)
    return [lookup[s] for s in syms]
end

"""Returns the constraints `NamedTuple` of `obj` (`(; move, rotate)`, either field possibly
missing), or `(;)` if `obj` has no entry in `ctrl.constraints`."""
_constraints_of(ctrl::KinematicController, obj) = get(ctrl.constraints, obj, (;))

"""Returns the allowed axes (a tuple of `:x`, `:y`, `:v`) of `obj` for `kind` (`:move` or
`:rotate`), all axes by default."""
function _allowed_axes(ctrl::KinematicController, obj, kind::Symbol)
    return get(_constraints_of(ctrl, obj), kind, _GIZMO_AXES)
end

"""Validates the `constraints` kwarg of [`kinematic_controls!`](@ref)."""
function _validate_constraints(constraints)
    for (obj, c) in constraints
        for k in keys(c)
            k in (:move, :rotate) ||
                throw(ArgumentError("invalid constraints field :$k, use :move or :rotate"))
        end
        for kind in (:move, :rotate)
            haskey(c, kind) || continue
            for a in c[kind]
                a in _GIZMO_AXES ||
                    throw(ArgumentError("invalid constraint axis :$a for :$kind, use :x, :y or :v"))
            end
        end
    end
    return nothing
end

"""Returns the colors of the gizmo axes `[y, x, v]` (green, red, blue) of `obj`'s `kind` controls
(`:move` or `:rotate`), faded to indicate axes locked by the `constraints`."""
function _gizmo_colors(ctrl::KinematicController, obj, kind::Symbol)
    allowed = _allowed_axes(ctrl, obj, kind)
    colors = Makie.RGBAf[]
    for (sym, c) in zip(_AXES_SYMS, _AXES_COLORS)
        rgba = Makie.RGBAf(Makie.to_color(c))
        push!(colors, sym in allowed ? rgba : Makie.RGBAf(rgba.r, rgba.g, rgba.b, _GIZMO_FADE_ALPHA))
    end
    return colors
end

function _update_help!(ctrl::KinematicController)
    if ctrl.spectator[]
        ctrl.help_obs[] = ctrl.help_shown ? _SPECTATOR_HELP : _SPECTATOR_HINT
        return nothing
    end
    help = _help_text(ctrl.mode[], ctrl.fine_step, ctrl.fine_angle, ctrl.select_modifier)
    isempty(ctrl.help_extra) || (help *= "\n" * ctrl.help_extra)
    ctrl.help_obs[] = ctrl.help_shown ? help :
                      _help_hint(ctrl.mode[], ctrl.fine_step, ctrl.fine_angle)
    return nothing
end

"""Changes the keyboard step of the current mode of the `ctrl` by one value of the 1-2-5 sequence."""
function _change_step!(ctrl::KinematicController, dir::Int)
    if ctrl.mode[] == :move
        ctrl.fine_step = clamp(_next_step(ctrl.fine_step, dir), _FINE_STEP_LIMITS...)
    else
        ctrl.fine_angle = clamp(_next_step(ctrl.fine_angle, dir), _FINE_ANGLE_LIMITS...)
    end
    _update_help!(ctrl)
    return nothing
end

_is_finite_box(bb) = all(isfinite, minimum(bb)) && all(isfinite, maximum(bb))

"""
    _selection_bbox(ctrl, obj, plots)

Returns the bounding box of the `plots` of the selected `obj`. With clip planes, see
[`live_view`](@ref), the bounding box of a plot only covers its visible part and is empty if the
plot is clipped entirely. Such plots are skipped; if all are clipped, a box around the `position`
of `obj` is returned, with the edge length of a source marker (8 % of the visible scene, or 1 cm).
"""
function _selection_bbox(ctrl::KinematicController, obj, plots)
    bbs = filter(_is_finite_box, [Makie.boundingbox(p) for p in plots])
    isempty(bbs) || return reduce(GeometryBasics.union, bbs)
    all_plots = reduce(vcat, (oh.plots for oh in ctrl.h.handles); init = AbstractPlot[])
    scene_bbs = filter(_is_finite_box, [Makie.boundingbox(p) for p in all_plots])
    w = isempty(scene_bbs) ? 1e-2 :
        0.08 * maximum(GeometryBasics.widths(reduce(GeometryBasics.union, scene_bbs)))
    return GeometryBasics.Rect3d(Vector{Float64}(position(obj)) .- w / 2, fill(w, 3))
end

function _update_selection_box!(ctrl::KinematicController)
    obj = ctrl.selected[]
    plots = isnothing(obj) ? AbstractPlot[] : _object_plots(ctrl.h, obj)
    if isempty(plots)
        isempty(ctrl.box_obs[]) || (empty!(ctrl.box_obs[]); notify(ctrl.box_obs))
        ctrl.gizmo_visible[] && (ctrl.gizmo_visible[] = false)
        return nothing
    end
    bb = _selection_bbox(ctrl, obj, plots)
    ctrl.box_obs[] = _bbox_wireframe(bb)
    # Place the gizmo above the object, where it is not covered by beams through the object
    w = GeometryBasics.widths(bb)
    l = (ctrl.mode[] == :move ? 1.2 : 0.8) * maximum(w)
    v = ctrl.rotation_axis
    offset = ctrl.mode[] == :move ? 0.3 * l : 1.4 * l
    origin = Vector{Float64}(position(obj)) + (dot(abs.(v), w) / 2 + offset) * v
    arrow_pos, arrow_dir, label_pos, ring_pts = _gizmo(ctrl.mode[], origin, _control_axes(ctrl, obj), l)
    ctrl.arrow_pos.val = arrow_pos
    ctrl.arrow_dir.val = arrow_dir
    ctrl.label_pos.val = label_pos
    ctrl.ring_pts.val = ring_pts
    ctrl.gizmo_size.val = l
    colors = _gizmo_colors(ctrl, obj, ctrl.mode[])
    ctrl.arrow_color.val = colors
    ctrl.label_color.val = colors
    ctrl.ring_color.val = repeat(colors; inner = 2 * _RING_RES)
    foreach(notify, (ctrl.gizmo_size, ctrl.arrow_pos, ctrl.arrow_dir, ctrl.label_pos, ctrl.ring_pts,
        ctrl.arrow_color, ctrl.label_color, ctrl.ring_color))
    ctrl.gizmo_visible[] || (ctrl.gizmo_visible[] = true)
    return nothing
end

function _request_update!(ctrl::KinematicController)
    if ctrl.throttle
        ctrl.dirty = true
    else
        _apply_update!(ctrl)
    end
    return nothing
end

function _apply_update!(ctrl::KinematicController)
    ctrl.dirty = false
    update_render!(ctrl.h)
    _update_selection_box!(ctrl)
    obj = ctrl.selected[]
    isnothing(obj) && return nothing
    # Errors must not propagate into the render loop
    try
        ctrl.on_change(obj)
        ctrl.last_error = nothing
    catch e
        msg = sprint(showerror, e)
        if msg != ctrl.last_error
            @error "kinematic_controls!: `on_change` callback failed" exception = (e, catch_backtrace())
            ctrl.last_error = msg
        end
    end
    return nothing
end

"""
    _set_pose!(obj, P, R)

Sets the pose of `obj` to position `P` and orientation `R`, by rotating around the axis-angle of
`R * orientation(obj)'` and then translating to `P`.
"""
function _set_pose!(obj, P, R)
    axis, angle = _axis_angle_from_rotmatrix(R * _pose(obj)[2]')
    angle > 1e-12 && rotate3d!(obj, axis, angle)
    translate_to3d!(obj, P)
    return nothing
end

function _reset_pose!(ctrl::KinematicController, obj)
    P0, R0 = ctrl.init_poses[obj]
    _set_pose!(obj, P0, R0)
    return nothing
end

const _HISTORY_LIMIT = 100

"""
    _push_history!(ctrl, obj, P0, R0, P1, R1)

Pushes a new undo entry moving `obj` from the pose `(P0, R0)` to `(P1, R1)`, unless nothing
changed. Starts a new gesture: clears the redo stack and drops the oldest entry once the history
exceeds `_HISTORY_LIMIT`.
"""
function _push_history!(ctrl::KinematicController, obj, P0, R0, P1, R1)
    (P0 == P1 && R0 == R1) && return nothing
    push!(ctrl.undo_stack, _HistoryEntry(obj, P0, R0, P1, R1))
    length(ctrl.undo_stack) > _HISTORY_LIMIT && popfirst!(ctrl.undo_stack)
    empty!(ctrl.redo_stack)
    return nothing
end

"""
    _record_key_step!(ctrl, obj, key, P0, R0, P1, R1)

Records a keyboard step of `obj` from `(P0, R0)` to `(P1, R1)` in the undo history. A step of the
same `key` on the same `obj` as the previous one, no more than 1 s apart, is merged into the same
entry instead of pushing a new one.
"""
function _record_key_step!(ctrl::KinematicController, obj, key, P0, R0, P1, R1)
    (P0 == P1 && R0 == R1) && return nothing
    now = time()
    m = ctrl.last_key_step
    if !isnothing(m) && m.obj === obj && m.key === key && (now - m.time) <= 1.0 &&
       !isempty(ctrl.undo_stack)
        e = pop!(ctrl.undo_stack)
        push!(ctrl.undo_stack, _HistoryEntry(obj, e.P0, e.R0, P1, R1))
    else
        _push_history!(ctrl, obj, P0, R0, P1, R1)
    end
    ctrl.last_key_step = (obj = obj, key = key, time = now)
    return nothing
end

"""
    _undo!(ctrl::KinematicController)

Undoes the last recorded gesture (a mouse drag, a reset or one or several merged keyboard steps),
if any, and selects its object. Returns whether an entry was undone.
"""
function _undo!(ctrl::KinematicController)
    isempty(ctrl.undo_stack) && return false
    e = pop!(ctrl.undo_stack)
    push!(ctrl.redo_stack, e)
    ctrl.last_key_step = nothing
    _set_pose!(e.obj, e.P0, e.R0)
    ctrl.selected[] === e.obj || (ctrl.selected[] = e.obj)
    _update_selection_box!(ctrl)
    _request_update!(ctrl)
    return true
end

"""
    _redo!(ctrl::KinematicController)

Redoes the last undone gesture, if any, and selects its object. Returns whether an entry was
redone.
"""
function _redo!(ctrl::KinematicController)
    isempty(ctrl.redo_stack) && return false
    e = pop!(ctrl.redo_stack)
    push!(ctrl.undo_stack, e)
    ctrl.last_key_step = nothing
    _set_pose!(e.obj, e.P1, e.R1)
    ctrl.selected[] === e.obj || (ctrl.selected[] = e.obj)
    _update_selection_box!(ctrl)
    _request_update!(ctrl)
    return true
end

"""
    _key_step!(ctrl, obj, key, factor)

Moves or rotates the `obj` depending on the mode of the `ctrl`. Returns `false` if the `key` is
not a control key. If the corresponding axis is locked by the `constraints` of `obj`, the key is
still consumed (returns `true`) but nothing moves.
"""
function _key_step!(ctrl::KinematicController, obj, key, factor)
    y, x, v = _control_axes(ctrl, obj)
    # Axis and sign of the step for each key, see _help_text
    axis_sym, axis, sign = if ctrl.mode[] == :move
        key == Keyboard.up ? (:y, y, 1) : key == Keyboard.down ? (:y, y, -1) :
        key == Keyboard.right ? (:x, x, 1) : key == Keyboard.left ? (:x, x, -1) :
        key == Keyboard.page_up ? (:v, v, 1) : key == Keyboard.page_down ? (:v, v, -1) :
        (nothing, nothing, 0)
    else
        key == Keyboard.up ? (:x, x, 1) : key == Keyboard.down ? (:x, x, -1) :
        key == Keyboard.left ? (:v, v, 1) : key == Keyboard.right ? (:v, v, -1) :
        key == Keyboard.page_up ? (:y, y, 1) : key == Keyboard.page_down ? (:y, y, -1) :
        (nothing, nothing, 0)
    end
    isnothing(axis) && return false
    axis_sym in _allowed_axes(ctrl, obj, ctrl.mode[]) || return true # locked: consumed, no motion
    if ctrl.mode[] == :move
        translate3d!(obj, (sign * factor * ctrl.fine_step) .* axis)
    else
        rotate3d!(obj, axis, sign * factor * ctrl.fine_angle)
    end
    return true
end

"""
    _set_spectator!(ctrl, on::Bool)

Switches the spectator mode of the `ctrl` on or off. In the spectator mode, the selection is
cleared and all mouse and keyboard input goes to the camera.
"""
function _set_spectator!(ctrl::KinematicController, on::Bool)
    ctrl.spectator[] = on
    if on
        ctrl.dragging = false
        ctrl.press_kind = :none
        ctrl.drag_start = nothing
        ctrl.last_key_step = nothing
        ctrl.selected[] = nothing
        _update_selection_box!(ctrl)
    end
    _update_help!(ctrl)
    return nothing
end

"""Returns the drag plane intersection of the ray through the mouse position."""
function _mouse_plane_hit(scene, ctrl::KinematicController)
    ray = Makie.ray_at_cursor(scene)
    return _ray_plane_intersect(Vector{Float64}(ray.origin), Vector{Float64}(ray.direction),
        ctrl.plane_point, ctrl.plane_normal)
end

_default_pick(ax) = Makie.pick(Makie.get_scene(ax))

"""
    _ray_box(origin, dir, bb)

Returns the distance along the ray from `origin` along `dir` to the box `bb`, `0` if `origin` is
inside of the box, or `nothing` if the ray misses the box.
"""
function _ray_box(origin, dir, bb)
    lo, hi = minimum(bb), maximum(bb)
    tmin, tmax = -Inf, Inf
    for i in 1:3
        if abs(dir[i]) < 1e-12
            lo[i] <= origin[i] <= hi[i] || return nothing
        else
            t1, t2 = (lo[i] - origin[i]) / dir[i], (hi[i] - origin[i]) / dir[i]
            tmin, tmax = max(tmin, min(t1, t2)), min(tmax, max(t1, t2))
        end
    end
    return tmax < max(tmin, 0) ? nothing : max(tmin, 0)
end

"""
    _ray_pick(objects, origin, dir, box = obj -> nothing)

Returns `(obj, t)`: the object among `objects` with the nearest hit of the ray from `origin` along
`dir`, and the hit distance `t`, or `(nothing, nothing)` if the ray misses all `objects`. Objects
are hit via `BMO.intersect3d`, i.e. `ObjectGroup`s and `MultiShape` objects are picked as a whole.
Objects for which `intersect3d` errors (e.g. custom objects without a geometry implementation) or
returns `nothing` (e.g. `NonInteractableObject`) are skipped. Sources, which have no geometry, are
hit via the bounding box returned by `box`.
"""
function _ray_pick(objects, origin, dir, box = obj -> nothing)
    ray = BMO.Ray(Vector{Float64}(origin), Vector{Float64}(dir))
    best, tmin = nothing, Inf
    for obj in objects
        t = if obj isa BMO.AbstractObject
            isect = try
                BMO.intersect3d(obj, ray)
            catch e
                e isa InterruptException && rethrow()
                nothing
            end
            isnothing(isect) ? nothing : BMO.length(isect)
        else
            bb = box(obj)
            isnothing(bb) ? nothing : _ray_box(origin, dir, bb)
        end
        isnothing(t) && continue
        if 0 <= t < tmin
            best, tmin = obj, t
        end
    end
    return best, isnothing(best) ? nothing : tmin
end

"""
    _ray_pick(ctrl::KinematicController, scene)

Returns `(obj, t)`: the object hit by the camera ray at the current cursor position of `scene`
among all objects of the movable objects of `ctrl`, i.e. objects of groups are returned instead of
the groups, and the hit distance `t`; or `(nothing, nothing)` if the ray misses everything. Sources
(non-`AbstractObject` leaves, i.e. beams and beam groups) are hit via the bounding box of their
marker as usual, and additionally whenever the screen-space projection of their `position` is
within `ctrl.source_pick_radius` pixels of the cursor, using as `t` the distance along the ray to
the point of the ray closest to the source position (so that the nearest candidate still wins).
"""
function _ray_pick(ctrl::KinematicController, scene)
    r = Makie.ray_at_cursor(scene)
    origin, dir = Vector{Float64}(r.origin), Vector{Float64}(r.direction)
    leaves = reduce(vcat, (_leaves(o) for o in ctrl.movable); init = _LiveMovable[])
    box = function (obj)
        plots = _object_plots(ctrl.h, obj)
        return isempty(plots) ? nothing : mapreduce(Makie.boundingbox, GeometryBasics.union, plots)
    end
    best, tmin = _ray_pick(leaves, origin, dir, box)
    tmin = isnothing(tmin) ? Inf : tmin
    cursor = _px(scene)
    for obj in leaves
        obj isa BMO.AbstractObject && continue # sources only, already covered by the bbox hit above
        p = Vector{Float64}(position(obj))
        px = Makie.project(scene, :data, :pixel, Point3(p))
        hypot(px[1] - cursor[1], px[2] - cursor[2]) <= ctrl.source_pick_radius || continue
        t = dot(p .- origin, dir)
        if 0 <= t < tmin
            best, tmin = obj, t
        end
    end
    return best, isnothing(best) ? nothing : tmin
end

"""
    kinematic_controls!(ax, h::SystemRenderHandle; kwargs...)

Enables mouse and keyboard controls for the objects of the live-rendered system `h`, see
[`live_render!`](@ref). Objects can be grabbed with the mouse and moved or rotated, while the
[`update_render!`](@ref) of `h` and the `on_change` callback are called after each change. The camera
can be used as usual as long as no object is grabbed. Returns a `KinematicController`.

The controls have a move and a rotate mode, which are switched with the key `m`. The selected object
is marked by a box and three axes above the object: its local y-axis (green), its local x-axis
(red) and the `rotation_axis` (blue). In the move mode the axes are shown as arrows, in the rotate
mode as rings. The key `h` shows or hides an overlay of all controls.

The key `v` switches the spectator mode on or off, in which the selection is cleared and all
mouse and keyboard input goes to the camera, such that nothing can be moved by accident.

Objects with a `Static` kinematic trait, see `BeamletOptics.kinematic_trait_of`, can not be
selected. Sources can be moved via their marker, see [`live_view`](@ref).

# Mouse controls

A click selects, a drag moves only what is already selected, every other drag rotates the camera:

- left-click on an object: selects it, so that a camera drag never moves or rotates a component by
  accident
- left-drag on the selected object: moves it within the plane through the grabbed point
  (`plane_normal`), such that the point under the cursor at the start of the drag stays under the
  cursor; or rotates it around the `rotation_axis` in the rotate mode. If the exact point under the
  cursor is not known (a custom `pick` function, or the `Makie.pick` fallback), the object's
  `position` is used instead
- left-drag elsewhere (background or an unselected object): rotates the camera as usual
- left-click on empty space: deselects the current object

Object groups (e.g. `ObjectGroup`) are selected as a whole by the first click. Each further click on
the selected group selects the next level of the hierarchy towards the object under the cursor,
i.e. a subgroup or a single object, which is then moved on its own. A group is moved and rotated
around its `position`, the group center.

Sources (beams and beam groups) are additionally picked when their `position`, projected to pixel
space, is within `source_pick_radius` of the cursor, in addition to the bounding box of their
marker, see [`live_view`](@ref).

# Keyboard controls

The following keys apply to the selected object, pressing shift multiplies the step size by 10:

| key                   | move mode (`fine_step`) | rotate mode (`fine_angle`) |
|:----------------------|:------------------------|:---------------------------|
| `↑`/`↓`               | along the green arrow   | around the red ring        |
| `→`/`←`               | along the red arrow     | around the blue ring       |
| `page up`/`page down` | along the blue arrow    | around the green ring      |

The first key moves the object in the direction of the arrow, or rotates it in the direction of
the ring. In addition, `Backspace` resets the object to its initial pose and `esc` deselects it. If
the object is part of a group, `esc` selects the enclosing group instead.

The keys `+` and `-` increase or decrease the step size of the current mode (`fine_step` or
`fine_angle`) along the 1-2-5 sequence, e.g. 10 nm → 20 nm → 50 nm → 100 nm. They also work
without a selected object. The current step size is shown in the hint line.

# Undo/redo

`Ctrl`+`Z` (`Cmd`+`Z` on macOS) undoes the last gesture, `Ctrl`+`Y` or `Ctrl`+`Shift`+`Z`
(`Cmd`+`Y`/`Cmd`+`Shift`+`Z`) redoes it. A gesture is a mouse drag, a reset (`Backspace`) or a run
of keyboard steps of the same key on the same object less than 1 s apart, each recorded as a single
history entry (up to 100 entries). Undoing or redoing selects the affected object and re-triggers
`on_change`. A new gesture clears the redo history. Not available in the spectator mode.

# Constraints

`constraints` locks the move and/or rotate axes of specific objects, e.g. a mirror in a kinematic
mount that may only tilt. A locked axis does not move via the mouse or the keyboard (the key is
still consumed) and is shown faded on the gizmo. For an object inside a group, the constraints of
the object (or subgroup) that is currently selected apply.

# Keyword args

- `objects = nothing`: the movable top-level objects, all top-level objects of `h` by default. The
  objects of a movable group are movable as well.
- `on_change = obj -> nothing`: called with the moved object after each change, e.g. to solve the system
- `plane_normal = [0, 0, 1]`: normal of the plane for mouse translation
- `rotation_axis = [0, 0, 1]`: rotation axis for mouse rotation, blue axis of the keyboard controls
- `rotate_speed = deg2rad(0.5)`: mouse rotation angle per pixel [rad]
- `fine_step = 10e-9`: keyboard translation step [m]
- `fine_angle = 10e-6`: keyboard rotation step [rad]
- `mode = :move`: initial mode, `:move` or `:rotate`
- `throttle = true`: limits updates to one per frame
- `show_help = false`: shows the controls overlay initially, otherwise only a hint
- `pick = nothing`: objects are selected by intersecting the camera ray with the movable objects, so
  meshes that are not part of the system (e.g. housings) do not block the selection; otherwise a
  function `ax -> (plot, index)`. The grab point of a drag falls back to the object's `position`
  when a custom `pick` is used, since it does not provide a hit point.
- `select_modifier = nothing`: if set to a `Keyboard.Button` (or a tuple/vector of them, meaning
  "any of"), clicks and drags only react on objects while the modifier is held; otherwise they
  always go to the camera. Suggested: `Keyboard.left_shift` (`Ctrl`+click is taken by Makie's
  `Camera3D` reset).
- `drag_threshold = 3`: [px] mouse movement between press and release that turns a click into a
  drag
- `spectator = false`: starts in the spectator mode
- `constraints = Dict()`: maps an object to a `NamedTuple` `(; move = axes, rotate = axes)`, where
  `axes` is a tuple of the allowed gizmo axes `:x` (red, local x), `:y` (green, local y) and `:v`
  (blue, `rotation_axis`); a missing field means all axes are allowed, `()` means none, see
  "Constraints"
- `source_pick_radius = 15`: [px] screen-space pick radius of source markers, see "Mouse controls"
- `ignore_keys = () -> false`: predicate checked before any keyboard handling (including `v`, `h`,
  `m` and undo/redo); while `true`, keys are not handled and passed on, e.g. while a text field
  elsewhere in the figure is focused
"""
function kinematic_controls!(
        ax::_RenderEnv,
        h::SystemRenderHandle;
        objects = nothing,
        on_change = obj -> nothing,
        plane_normal = [0, 0, 1],
        rotation_axis = [0, 0, 1],
        rotate_speed = deg2rad(0.5),
        fine_step = 10e-9,
        fine_angle = 10e-6,
        mode::Symbol = :move,
        throttle::Bool = true,
        show_help::Bool = false,
        pick = nothing,
        select_modifier = nothing,
        drag_threshold = 3,
        spectator::Bool = false,
        constraints = Dict(),
        source_pick_radius = 15,
        ignore_keys = () -> false
    )
    mode in (:move, :rotate) || throw(ArgumentError("mode must be :move or :rotate, got :$mode"))
    # Objects are compared by identity
    constraints_dict = IdDict{Any, NamedTuple}(constraints)
    _validate_constraints(constraints_dict)
    scene = Makie.get_scene(ax)
    movable = _LiveMovable[]
    if isnothing(objects)
        # Top-level objects, since the handles of groups belong to the objects of the group
        for oh in h.handles
            top = _top_level(h, oh.obj)
            BMO._is_static(top) && continue
            any(o -> o === top, movable) || push!(movable, top)
        end
    else
        append!(movable, objects)
    end
    # Initial poses of all levels, such that each object and subgroup can be reset
    init_poses = IdDict{_LiveMovable, Tuple{Point3{Float64}, Matrix{Float64}}}(
        obj => _pose(obj) for top in movable for obj in _descendants(top))

    # Selection box and gizmo, updated via Observables. The hidden gizmo is placed in the center of
    # the system with a negligible size, since it counts towards the limits of the scene.
    box_obs = Observable(Point3f[])
    obj_plots = reduce(vcat, (oh.plots for oh in h.handles); init = AbstractPlot[])
    bb = isempty(obj_plots) ? GeometryBasics.Rect3d(zeros(3), ones(3)) :
         mapreduce(Makie.boundingbox, GeometryBasics.union, obj_plots)
    center = Vector{Float64}(minimum(bb) + GeometryBasics.widths(bb) / 2)
    arrow_pos, arrow_dir, label_pos, ring_pts = _gizmo(mode, center,
        ([0, 1, 0], [1, 0, 0], [0, 0, 1]), 1e-3 * maximum(GeometryBasics.widths(bb)))
    arrow_pos, arrow_dir = Observable(arrow_pos), Observable(arrow_dir)
    label_pos, ring_pts = Observable(label_pos), Observable(ring_pts)
    default_colors = Makie.RGBAf[Makie.RGBAf(Makie.to_color(c)) for c in _AXES_COLORS]
    arrow_color = Observable(default_colors)
    label_color = Observable(default_colors)
    ring_color = Observable(repeat(default_colors; inner = 2 * _RING_RES))
    gizmo_size = Observable(1.0)
    gizmo_visible = Observable(false)
    mode_obs = Observable(mode)
    labels = Makie.lift(m -> m == :move ? ["↑", "→", "page up"] : ["page up", "↑", "←"], mode_obs)
    plots = AbstractPlot[
        linesegments!(ax, box_obs; color = :yellow, linewidth = 2, clip_planes = Plane3f[]),
        # Arrow dimensions relative to the gizmo size, such that ring arrows are mostly tip
        arrows3d!(ax, arrow_pos, arrow_dir; color = arrow_color, visible = gizmo_visible,
            markerscale = gizmo_size, shaftradius = 0.025, tipradius = 0.1, tiplength = 0.3,
            overdraw = true, clip_planes = Plane3f[]),
        linesegments!(ax, ring_pts; color = ring_color, linewidth = 3, overdraw = true,
            visible = Makie.lift((v, m) -> v && m == :rotate, gizmo_visible, mode_obs),
            clip_planes = Plane3f[]),
        text!(ax, label_pos; text = labels, color = label_color, visible = gizmo_visible,
            fontsize = 20, align = (:center, :center), overdraw = true, clip_planes = Plane3f[])
    ]
    help_obs = Observable("")
    # Drawn in the 2D scene of the axis, such that it does not count towards the limits of the scene
    help_pos = Makie.lift(vp -> Point2f(minimum(vp)[1] + 10, maximum(vp)[2] - 10), ax.scene.viewport)
    push!(plots, text!(ax.blockscene, help_pos; text = help_obs, space = :pixel,
        align = (:left, :top), fontsize = 14, color = :gray40))

    ctrl = KinematicController(
        ax, h, movable, init_poses, Observable{Union{Nothing, _LiveMovable}}(nothing),
        mode_obs, on_change, normalize(Float64.(plane_normal)), normalize(Float64.(rotation_axis)),
        Float64(rotate_speed), Float64(fine_step), Float64(fine_angle), throttle,
        Observable(spectator), select_modifier, Float64(drag_threshold),
        constraints_dict, Float64(source_pick_radius), ignore_keys,
        false, false, zeros(3), zeros(3), (0.0, 0.0), nothing, nothing, :none,
        nothing, _HistoryEntry[], _HistoryEntry[], nothing,
        box_obs, arrow_pos, arrow_dir, label_pos, ring_pts, arrow_color, label_color, ring_color,
        gizmo_size, gizmo_visible, help_obs, show_help, "", plots, Any[], nothing
    )

    # High priority, so that the camera does not receive events while an object is dragged
    l1 = on(events(scene).mousebutton, priority = 200) do event
        (event.button == Mouse.left && !ctrl.spectator[]) || return Consume(false)
        if event.action == Mouse.press
            if !isnothing(ctrl.select_modifier) && !_modifier_held(scene, ctrl.select_modifier)
                # Modifier not held: every click and drag goes to the camera, no state change
                return Consume(false)
            end
            local t
            if pick === nothing
                leaf, t = _ray_pick(ctrl, scene)
                if isnothing(leaf)
                    plot, _ = _default_pick(ax)
                    leaf = isnothing(plot) ? nothing : _pick_leaf(h, plot)
                    t = nothing
                end
            else
                plot, _ = pick(ax)
                leaf = isnothing(plot) ? nothing : _pick_leaf(h, plot)
                t = nothing
            end
            ctrl.press_pos = _px(scene)
            if isnothing(leaf) || !_is_movable(ctrl, leaf)
                ctrl.press_leaf = nothing
                ctrl.press_kind = :background
                return Consume(false)
            end
            ctrl.press_leaf = leaf
            sel = ctrl.selected[]
            if !isnothing(sel) && any(o -> o === sel, _chain(ctrl, leaf))
                # Dragging the already selected object: block the camera immediately
                ctrl.press_kind = :pending_drag
                ctrl.last_mouse = ctrl.press_pos
                ctrl.plane_point = Vector{Float64}(position(sel))
                # Grab at the exact point under the cursor, so that point stays under the cursor
                # during the drag; falls back to the pivot if the hit point is not known
                hit_point = nothing
                if !isnothing(t)
                    r = Makie.ray_at_cursor(scene)
                    hit_point = Vector{Float64}(r.origin) .+ t .* Vector{Float64}(r.direction)
                    ctrl.plane_point = hit_point
                end
                if !isnothing(hit_point)
                    ctrl.grab_offset = Vector{Float64}(position(sel)) .- hit_point
                else
                    hit = _mouse_plane_hit(scene, ctrl)
                    ctrl.grab_offset = isnothing(hit) ? zeros(3) : ctrl.plane_point .- hit
                end
                return Consume(true)
            else
                # A press on an unselected object might turn into a camera drag: let it through
                ctrl.press_kind = :pending_select
                return Consume(false)
            end
        elseif event.action == Mouse.release
            cur = _px(scene)
            moved = isnothing(ctrl.press_pos) ? Inf : hypot((cur .- ctrl.press_pos)...)
            kind = ctrl.press_kind
            consume = false
            if ctrl.dragging
                ctrl.dragging = false
                consume = true
                if !isnothing(ctrl.drag_start)
                    obj = ctrl.selected[]
                    P0, R0 = ctrl.drag_start
                    P1, R1 = _pose(obj)
                    ctrl.drag_start = nothing
                    ctrl.last_key_step = nothing
                    _push_history!(ctrl, obj, P0, R0, P1, R1)
                end
            elseif kind == :pending_drag
                # Released before crossing the threshold: a click, not a drag
                if moved < ctrl.drag_threshold
                    ctrl.selected[] = _drill_select(ctrl, ctrl.press_leaf)
                    _update_selection_box!(ctrl)
                end
                consume = true
            elseif kind == :pending_select
                # Only select if this was a click, not a camera rotation
                if moved < ctrl.drag_threshold
                    ctrl.selected[] = _drill_select(ctrl, ctrl.press_leaf)
                    _update_selection_box!(ctrl)
                end
            elseif kind == :background
                if moved < ctrl.drag_threshold
                    ctrl.selected[] = nothing
                    _update_selection_box!(ctrl)
                end
            end
            ctrl.press_kind = :none
            ctrl.press_leaf = nothing
            ctrl.press_pos = nothing
            return Consume(consume)
        end
        return Consume(false)
    end

    l2 = on(events(scene).mouseposition, priority = 200) do _
        if !ctrl.dragging
            ctrl.press_kind == :pending_drag || return Consume(false)
            cur = _px(scene)
            moved = hypot((cur .- ctrl.press_pos)...)
            moved >= ctrl.drag_threshold || return Consume(true)
            ctrl.dragging = true
            ctrl.drag_start = _pose(ctrl.selected[])
        end
        obj = ctrl.selected[]
        if ctrl.mode[] == :move
            hit = _mouse_plane_hit(scene, ctrl)
            if !isnothing(hit)
                target = hit .+ ctrl.grab_offset
                allowed = _allowed_axes(ctrl, obj, :move)
                if !isempty(allowed)
                    Δ = target .- Vector{Float64}(position(obj))
                    A = hcat(_axis_vectors(ctrl, obj, allowed)...)
                    # Projection onto the allowed axes, which may be linearly dependent, e.g. if
                    # a local axis is parallel to the rotation axis
                    translate3d!(obj, A * (pinv(A) * Δ))
                    _request_update!(ctrl)
                end
            end
        else
            mp = _px(scene)
            dx = mp[1] - ctrl.last_mouse[1]
            ctrl.last_mouse = mp
            if dx != 0 && :v in _allowed_axes(ctrl, obj, :rotate)
                rotate3d!(obj, ctrl.rotation_axis, ctrl.rotate_speed * dx)
                _request_update!(ctrl)
            end
        end
        return Consume(true)
    end

    l3 = on(events(scene).keyboardbutton, priority = 200) do event
        ctrl.ignore_keys() && return Consume(false)
        event.action in (Keyboard.press, Keyboard.repeat) || return Consume(false)
        if event.action == Keyboard.press && event.key == Keyboard.v
            _set_spectator!(ctrl, !ctrl.spectator[])
            return Consume(true)
        end
        # Only the overlay can be toggled in the spectator mode
        ctrl.spectator[] && event.key != Keyboard.h && return Consume(false)
        if event.action == Keyboard.press && event.key in (Keyboard.h, Keyboard.m)
            if event.key == Keyboard.h
                ctrl.help_shown = !ctrl.help_shown
            else
                ctrl.mode[] = _other_mode(ctrl.mode[])
                _update_selection_box!(ctrl)
            end
            _update_help!(ctrl)
            return Consume(true)
        end
        if event.key in (Keyboard.z, Keyboard.y) && _modifier_held(scene, _CTRL_OR_CMD_KEYS)
            if event.key == Keyboard.y || _shift_pressed(scene)
                _redo!(ctrl)
            else
                _undo!(ctrl)
            end
            return Consume(true)
        end
        obj = ctrl.selected[]
        isnothing(obj) && return Consume(false)
        if event.key == Keyboard.escape
            # One level up in the hierarchy of groups, deselect at the top level
            ctrl.selected[] = get(ctrl.h.parent, obj, nothing)
            _update_selection_box!(ctrl)
        elseif event.key == Keyboard.backspace
            P0, R0 = _pose(obj)
            _reset_pose!(ctrl, obj)
            P1, R1 = _pose(obj)
            ctrl.last_key_step = nothing
            _push_history!(ctrl, obj, P0, R0, P1, R1)
            _request_update!(ctrl)
        else
            P0, R0 = _pose(obj)
            if _key_step!(ctrl, obj, event.key, _shift_pressed(scene) ? 10 : 1)
                P1, R1 = _pose(obj)
                _record_key_step!(ctrl, obj, event.key, P0, R0, P1, R1)
                _request_update!(ctrl)
            else
                return Consume(false)
            end
        end
        return Consume(true)
    end
    # Typed characters instead of keys, since + and - depend on the keyboard layout
    l4 = on(events(scene).unicode_input, priority = 200) do char
        ctrl.ignore_keys() && return Consume(false)
        (char in ('+', '-') && !ctrl.spectator[]) || return Consume(false)
        _change_step!(ctrl, char == '+' ? 1 : -1)
        return Consume(true)
    end
    push!(ctrl.listeners, l1, l2, l3, l4)
    _update_help!(ctrl)

    if throttle
        push!(ctrl.listeners, on(_ -> ctrl.dirty && _apply_update!(ctrl), events(scene).tick))
    end
    return ctrl
end

function Base.close(ctrl::KinematicController)
    foreach(off, ctrl.listeners)
    empty!(ctrl.listeners)
    foreach(p -> delete!(p.parent, p), ctrl.plots)
    empty!(ctrl.plots)
    return nothing
end
