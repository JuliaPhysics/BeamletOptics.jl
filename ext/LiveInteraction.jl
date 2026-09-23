import Makie
using Makie: events, on, off, Consume, Mouse, Keyboard, mouseposition_px, Point2f, Vec3f
using LinearAlgebra: tr, diag, I

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

function _help_hint(mode::Symbol, fine_step, fine_angle)
    step = _step_string(mode, fine_step, fine_angle)
    return "$mode mode, step $step, +/-: step, m: switch mode, h: show controls"
end

function _help_text(mode::Symbol, fine_step, fine_angle)
    step = _step_string(mode, fine_step, fine_angle)
    if mode == :move
        verb = "move along"
        up, left, page, drag = "green arrow", "red arrow", "blue arrow", "move in plane"
    else
        verb = "rotate around"
        up, left, page, drag = "red ring", "blue ring", "green ring", "rotate around blue ring"
    end
    return """
    $mode mode, m: switch to $(_other_mode(mode)) mode
    left-drag: $drag
    ↑/↓: $verb $up
    ←/→: $verb $left
    page up/down: $verb $page
    step: $step, shift: 10× step
    +/-: change step
    r: reset, esc: deselect
    click again: select part of a group, esc: up one level
    h: hide controls"""
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

"""
    KinematicController

Returned by [`kinematic_controls!`](@ref). The currently selected object is stored in the
`selected` `Observable`, the current mode (`:move` or `:rotate`) in the `mode` `Observable`. Use
`close` to remove the controls.
"""
mutable struct KinematicController{H <: SystemRenderHandle}
    ax::_RenderEnv
    h::H
    movable::Vector{BMO.AbstractObject}
    init_poses::IdDict{BMO.AbstractObject, Tuple{Point3{Float64}, Matrix{Float64}}}
    selected::Observable{Union{Nothing, BMO.AbstractObject}}
    mode::Observable{Symbol}
    on_change::Function
    plane_normal::Vector{Float64}
    rotation_axis::Vector{Float64}
    rotate_speed::Float64
    fine_step::Float64
    fine_angle::Float64
    throttle::Bool
    # interaction state
    dirty::Bool
    dragging::Bool
    plane_point::Vector{Float64}
    grab_offset::Vector{Float64}
    last_mouse::NTuple{2, Float64}
    bg_press_pos::Union{Nothing, NTuple{2, Float64}}
    # selection box and gizmo of the keyboard controls
    box_obs::Observable{Vector{Point3f}}
    arrow_pos::Observable{Vector{Point3f}}
    arrow_dir::Observable{Vector{Vec3f}}
    label_pos::Observable{Vector{Point3f}}
    ring_pts::Observable{Vector{Point3f}}
    gizmo_size::Observable{Float64}
    gizmo_visible::Observable{Bool}
    # controls overlay, toggled via h
    help_obs::Observable{String}
    help_shown::Bool
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
_children(obj) = obj isa BMO.AbstractObjectGroup ? collect(BMO.AbstractObject, BMO.shape(obj)) :
                 BMO.AbstractObject[]

"""Returns all objects of the group `obj` that are not groups themselves (recursively), or `[obj]`."""
function _leaves(obj)
    obj isa BMO.AbstractObjectGroup || return BMO.AbstractObject[obj]
    return reduce(vcat, (_leaves(c) for c in BMO.shape(obj)); init = BMO.AbstractObject[])
end

"""Returns the `obj` and all its nested objects and subgroups (recursively)."""
function _descendants(obj)
    out = BMO.AbstractObject[obj]
    for c in _children(obj)
        append!(out, _descendants(c))
    end
    return out
end

"""Returns the chain `[leaf, parent of leaf, …, top-level object]` of the hierarchy of `ctrl.h`."""
function _chain(ctrl, leaf)
    chain = BMO.AbstractObject[leaf]
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
        isnothing(i) || append!(plots, h.handles[i].plots)
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
    R = orientation(obj)
    return (Vector{Float64}(R[:, 2]), Vector{Float64}(R[:, 1]), ctrl.rotation_axis)
end

function _update_help!(ctrl::KinematicController)
    ctrl.help_obs[] = ctrl.help_shown ? _help_text(ctrl.mode[], ctrl.fine_step, ctrl.fine_angle) :
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

function _update_selection_box!(ctrl::KinematicController)
    obj = ctrl.selected[]
    plots = isnothing(obj) ? AbstractPlot[] : _object_plots(ctrl.h, obj)
    if isempty(plots)
        isempty(ctrl.box_obs[]) || (empty!(ctrl.box_obs[]); notify(ctrl.box_obs))
        ctrl.gizmo_visible[] && (ctrl.gizmo_visible[] = false)
        return nothing
    end
    bb = mapreduce(Makie.boundingbox, GeometryBasics.union, plots)
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
    foreach(notify, (ctrl.gizmo_size, ctrl.arrow_pos, ctrl.arrow_dir, ctrl.label_pos, ctrl.ring_pts))
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

function _reset_pose!(ctrl::KinematicController, obj)
    P0, R0 = ctrl.init_poses[obj]
    axis, angle = _axis_angle_from_rotmatrix(R0 * orientation(obj)')
    angle > 1e-12 && rotate3d!(obj, axis, angle)
    translate_to3d!(obj, P0)
    return nothing
end

"""
    _key_step!(ctrl, obj, key, factor)

Moves or rotates the `obj` depending on the mode of the `ctrl`. Returns `false` if the `key` is
not a control key.
"""
function _key_step!(ctrl::KinematicController, obj, key, factor)
    y, x, v = _control_axes(ctrl, obj)
    # Axis and sign of the step for each key, see _help_text
    axis, sign = if ctrl.mode[] == :move
        key == Keyboard.up ? (y, 1) : key == Keyboard.down ? (y, -1) :
        key == Keyboard.right ? (x, 1) : key == Keyboard.left ? (x, -1) :
        key == Keyboard.page_up ? (v, 1) : key == Keyboard.page_down ? (v, -1) : (nothing, 0)
    else
        key == Keyboard.up ? (x, 1) : key == Keyboard.down ? (x, -1) :
        key == Keyboard.left ? (v, 1) : key == Keyboard.right ? (v, -1) :
        key == Keyboard.page_up ? (y, 1) : key == Keyboard.page_down ? (y, -1) : (nothing, 0)
    end
    isnothing(axis) && return false
    if ctrl.mode[] == :move
        translate3d!(obj, (sign * factor * ctrl.fine_step) .* axis)
    else
        rotate3d!(obj, axis, sign * factor * ctrl.fine_angle)
    end
    return true
end

"""Returns the drag plane intersection of the ray through the mouse position."""
function _mouse_plane_hit(scene, ctrl::KinematicController)
    ray = Makie.ray_at_cursor(scene)
    return _ray_plane_intersect(Vector{Float64}(ray.origin), Vector{Float64}(ray.direction),
        ctrl.plane_point, ctrl.plane_normal)
end

_default_pick(ax) = Makie.pick(Makie.get_scene(ax))

"""
    _ray_pick(objects, origin, dir)

Returns the object among `objects` with the nearest `BMO.intersect3d` hit of the ray from `origin`
along `dir`, or `nothing` if the ray misses all `objects`. `ObjectGroup`s and `MultiShape` objects
are picked as a whole through `intersect3d`. Objects for which `intersect3d` errors (e.g. custom
objects without a geometry implementation) or returns `nothing` (e.g. `NonInteractableObject`) are
skipped.
"""
function _ray_pick(objects, origin, dir)
    ray = BMO.Ray(Vector{Float64}(origin), Vector{Float64}(dir))
    best, tmin = nothing, Inf
    for obj in objects
        isect = try
            BMO.intersect3d(obj, ray)
        catch
            nothing
        end
        isnothing(isect) && continue
        t = BMO.length(isect)
        if 0 < t < tmin
            best, tmin = obj, t
        end
    end
    return best
end

"""
    _ray_pick(ctrl::KinematicController, scene)

Returns the object hit by the camera ray at the current cursor position of `scene` among all
objects of the movable objects of `ctrl`, i.e. objects of groups are returned instead of the
groups, see [`_ray_pick(objects, origin, dir)`](@ref).
"""
function _ray_pick(ctrl::KinematicController, scene)
    r = Makie.ray_at_cursor(scene)
    leaves = reduce(vcat, (_leaves(o) for o in ctrl.movable); init = BMO.AbstractObject[])
    return _ray_pick(leaves, r.origin, r.direction)
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

# Mouse controls

- left-drag on an object: moves the object within the plane through its position (`plane_normal`),
  or rotates it around the `rotation_axis` in the rotate mode
- left-click on empty space: deselects the current object

Object groups (e.g. `ObjectGroup`) are selected as a whole by the first click. Each further click on
the selected group selects the next level of the hierarchy towards the object under the cursor,
i.e. a subgroup or a single object, which is then moved on its own. A group is moved and rotated
around its `position`, the group center.

# Keyboard controls

The following keys apply to the selected object, pressing shift multiplies the step size by 10:

| key                   | move mode (`fine_step`) | rotate mode (`fine_angle`) |
|:----------------------|:------------------------|:---------------------------|
| `↑`/`↓`               | along the green arrow   | around the red ring        |
| `→`/`←`               | along the red arrow     | around the blue ring       |
| `page up`/`page down` | along the blue arrow    | around the green ring      |

The first key moves the object in the direction of the arrow, or rotates it in the direction of
the ring. In addition, `r` resets the object to its initial pose and `esc` deselects it. If the
object is part of a group, `esc` selects the enclosing group instead.

The keys `+` and `-` increase or decrease the step size of the current mode (`fine_step` or
`fine_angle`) along the 1-2-5 sequence, e.g. 10 nm → 20 nm → 50 nm → 100 nm. They also work
without a selected object. The current step size is shown in the hint line.

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
  function `ax -> (plot, index)`
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
        pick = nothing
    )
    mode in (:move, :rotate) || throw(ArgumentError("mode must be :move or :rotate, got :$mode"))
    scene = Makie.get_scene(ax)
    movable = BMO.AbstractObject[]
    if isnothing(objects)
        # Top-level objects, since the handles of groups belong to the objects of the group
        for oh in h.handles
            top = _top_level(h, oh.obj)
            any(o -> o === top, movable) || push!(movable, top)
        end
    else
        append!(movable, objects)
    end
    # Initial poses of all levels, such that each object and subgroup can be reset
    init_poses = IdDict{BMO.AbstractObject, Tuple{Point3{Float64}, Matrix{Float64}}}(
        obj => _pose(obj) for top in movable for obj in _descendants(top))

    # Selection box and gizmo, updated via Observables
    box_obs = Observable(Point3f[])
    arrow_pos, arrow_dir, label_pos, ring_pts = _gizmo(mode, zeros(3), ([0, 1, 0], [1, 0, 0], [0, 0, 1]), 1.0)
    arrow_pos, arrow_dir = Observable(arrow_pos), Observable(arrow_dir)
    label_pos, ring_pts = Observable(label_pos), Observable(ring_pts)
    gizmo_size = Observable(1.0)
    gizmo_visible = Observable(false)
    mode_obs = Observable(mode)
    ring_colors = repeat(_AXES_COLORS; inner = 2 * _RING_RES)
    labels = Makie.lift(m -> m == :move ? ["↑", "→", "page up"] : ["page up", "↑", "←"], mode_obs)
    plots = AbstractPlot[
        linesegments!(ax, box_obs; color = :yellow, linewidth = 2),
        # Arrow dimensions relative to the gizmo size, such that ring arrows are mostly tip
        arrows3d!(ax, arrow_pos, arrow_dir; color = _AXES_COLORS, visible = gizmo_visible,
            markerscale = gizmo_size, shaftradius = 0.025, tipradius = 0.1, tiplength = 0.3,
            overdraw = true),
        linesegments!(ax, ring_pts; color = ring_colors, linewidth = 3, overdraw = true,
            visible = Makie.lift((v, m) -> v && m == :rotate, gizmo_visible, mode_obs)),
        text!(ax, label_pos; text = labels, color = _AXES_COLORS, visible = gizmo_visible,
            fontsize = 20, align = (:center, :center), overdraw = true)
    ]
    help_obs = Observable(show_help ? _help_text(mode, fine_step, fine_angle) :
                      _help_hint(mode, fine_step, fine_angle))
    push!(plots, text!(ax, Point2f(0.01, 0.99); text = help_obs, space = :relative,
        align = (:left, :top), fontsize = 14, color = :gray40))

    ctrl = KinematicController(
        ax, h, movable, init_poses, Observable{Union{Nothing, BMO.AbstractObject}}(nothing),
        mode_obs, on_change, normalize(Float64.(plane_normal)), normalize(Float64.(rotation_axis)),
        Float64(rotate_speed), Float64(fine_step), Float64(fine_angle), throttle,
        false, false, zeros(3), zeros(3), (0.0, 0.0), nothing,
        box_obs, arrow_pos, arrow_dir, label_pos, ring_pts, gizmo_size, gizmo_visible,
        help_obs, show_help, plots, Any[], nothing
    )

    # High priority, so that the camera does not receive events while an object is grabbed
    l1 = on(events(scene).mousebutton, priority = 200) do event
        event.button == Mouse.left || return Consume(false)
        if event.action == Mouse.press
            if pick === nothing
                leaf = _ray_pick(ctrl, scene)
                if isnothing(leaf)
                    plot, _ = _default_pick(ax)
                    leaf = isnothing(plot) ? nothing : _pick_leaf(h, plot)
                end
            else
                plot, _ = pick(ax)
                leaf = isnothing(plot) ? nothing : _pick_leaf(h, plot)
            end
            if isnothing(leaf) || !_is_movable(ctrl, leaf)
                ctrl.bg_press_pos = _px(scene)
                return Consume(false)
            end
            # Clicking again on a group selects the next level towards the clicked object
            obj = _drill_select(ctrl, leaf)
            ctrl.selected[] = obj
            ctrl.dragging = true
            ctrl.last_mouse = _px(scene)
            # Keep the offset between object and mouse to avoid a jump at the start
            ctrl.plane_point = Vector{Float64}(position(obj))
            hit = _mouse_plane_hit(scene, ctrl)
            ctrl.grab_offset = isnothing(hit) ? zeros(3) : ctrl.plane_point .- hit
            _update_selection_box!(ctrl)
            return Consume(true)
        elseif event.action == Mouse.release
            if ctrl.dragging
                ctrl.dragging = false
                return Consume(true)
            end
            # Deselect on click, but not after rotating the camera
            if !isnothing(ctrl.bg_press_pos)
                cur = _px(scene)
                moved = hypot((cur .- ctrl.bg_press_pos)...)
                ctrl.bg_press_pos = nothing
                if moved < 3
                    ctrl.selected[] = nothing
                    _update_selection_box!(ctrl)
                end
            end
        end
        return Consume(false)
    end

    l2 = on(events(scene).mouseposition, priority = 200) do _
        obj = ctrl.selected[]
        (isnothing(obj) || !ctrl.dragging) && return Consume(false)
        if ctrl.mode[] == :move
            hit = _mouse_plane_hit(scene, ctrl)
            if !isnothing(hit)
                translate_to3d!(obj, hit .+ ctrl.grab_offset)
                _request_update!(ctrl)
            end
        else
            mp = _px(scene)
            dx = mp[1] - ctrl.last_mouse[1]
            ctrl.last_mouse = mp
            if dx != 0
                rotate3d!(obj, ctrl.rotation_axis, ctrl.rotate_speed * dx)
                _request_update!(ctrl)
            end
        end
        return Consume(true)
    end

    l3 = on(events(scene).keyboardbutton, priority = 200) do event
        event.action in (Keyboard.press, Keyboard.repeat) || return Consume(false)
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
        obj = ctrl.selected[]
        isnothing(obj) && return Consume(false)
        if event.key == Keyboard.escape
            # One level up in the hierarchy of groups, deselect at the top level
            ctrl.selected[] = get(ctrl.h.parent, obj, nothing)
            _update_selection_box!(ctrl)
        elseif event.key == Keyboard.r
            _reset_pose!(ctrl, obj)
            _request_update!(ctrl)
        elseif _key_step!(ctrl, obj, event.key, _shift_pressed(scene) ? 10 : 1)
            _request_update!(ctrl)
        else
            return Consume(false)
        end
        return Consume(true)
    end
    # Typed characters instead of keys, since + and - depend on the keyboard layout
    l4 = on(events(scene).unicode_input, priority = 200) do char
        char in ('+', '-') || return Consume(false)
        _change_step!(ctrl, char == '+' ? 1 : -1)
        return Consume(true)
    end
    push!(ctrl.listeners, l1, l2, l3, l4)

    if throttle
        push!(ctrl.listeners, on(_ -> ctrl.dirty && _apply_update!(ctrl), events(scene).tick))
    end
    return ctrl
end

function Base.close(ctrl::KinematicController)
    foreach(off, ctrl.listeners)
    empty!(ctrl.listeners)
    foreach(p -> delete!(ctrl.ax, p), ctrl.plots)
    empty!(ctrl.plots)
    return nothing
end
