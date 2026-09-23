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

const _HELP_HINT = "h: show controls"

function _help_text(fine_step, fine_angle)
    step = round(fine_step * 1e9, sigdigits = 3)
    angle = round(fine_angle * 1e6, sigdigits = 3)
    return """
    left-drag: move
    shift + left-drag: rotate
    ↑/↓: move along green arrow by $step nm
    ←/→: rotate around blue arrow by $angle µrad
    page up/down: tilt around red arrow by $angle µrad
    shift + key: 10× step
    r: reset
    esc: deselect
    h: hide controls"""
end

"""
    KinematicController

Returned by [`kinematic_controls!`](@ref). The currently selected object is stored in the
`selected` `Observable`. Use `close` to remove the controls.
"""
mutable struct KinematicController{H <: SystemRenderHandle}
    ax::_RenderEnv
    h::H
    movable::Vector{BMO.AbstractObject}
    init_poses::IdDict{BMO.AbstractObject, Tuple{Point3{Float64}, Matrix{Float64}}}
    selected::Observable{Union{Nothing, BMO.AbstractObject}}
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
    rotating::Bool
    plane_point::Vector{Float64}
    grab_offset::Vector{Float64}
    last_mouse::NTuple{2, Float64}
    bg_press_pos::Union{Nothing, NTuple{2, Float64}}
    # selection box and axes of the keyboard controls
    box_obs::Observable{Vector{Point3f}}
    axes_pos::Observable{Vector{Point3f}}
    axes_dir::Observable{Vector{Vec3f}}
    axes_visible::Observable{Bool}
    # controls overlay, toggled via h
    help_obs::Observable{String}
    plots::Vector{AbstractPlot}
    listeners::Vector{Any}
    # last error of on_change, logged only once
    last_error::Union{Nothing, String}
end

function Base.show(io::IO, ctrl::KinematicController)
    obj = ctrl.selected[]
    sel = isnothing(obj) ? "none" : string(nameof(typeof(obj)))
    print(io, "KinematicController(", length(ctrl.movable), " movable, selected = ", sel, ")")
end

_is_movable(ctrl::KinematicController, obj) = any(o -> o === obj, ctrl.movable)

function _object_plots(h::SystemRenderHandle, obj)
    i = findfirst(oh -> oh.obj === obj, h.handles)
    return isnothing(i) ? AbstractPlot[] : h.handles[i].plots
end

function _update_selection_box!(ctrl::KinematicController)
    obj = ctrl.selected[]
    plots = isnothing(obj) ? AbstractPlot[] : _object_plots(ctrl.h, obj)
    if isempty(plots)
        isempty(ctrl.box_obs[]) || (empty!(ctrl.box_obs[]); notify(ctrl.box_obs))
        ctrl.axes_visible[] && (ctrl.axes_visible[] = false)
        return nothing
    end
    bb = mapreduce(Makie.boundingbox, GeometryBasics.union, plots)
    ctrl.box_obs[] = _bbox_wireframe(bb)
    # Move direction, tilt axis and rotation axis of the keyboard controls
    l = 1.5 * maximum(GeometryBasics.widths(bb))
    R = orientation(obj)
    ctrl.axes_pos[] = fill(Point3f(position(obj)), 3)
    ctrl.axes_dir[] = [Vec3f(l * R[:, 2]), Vec3f(l * R[:, 1]), Vec3f(l * ctrl.rotation_axis)]
    ctrl.axes_visible[] || (ctrl.axes_visible[] = true)
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

"""Returns the drag plane intersection of the ray through the mouse position."""
function _mouse_plane_hit(scene, ctrl::KinematicController)
    ray = Makie.ray_at_cursor(scene)
    return _ray_plane_intersect(Vector{Float64}(ray.origin), Vector{Float64}(ray.direction),
        ctrl.plane_point, ctrl.plane_normal)
end

_default_pick(ax) = Makie.pick(Makie.get_scene(ax))

"""
    kinematic_controls!(ax, h::SystemRenderHandle; kwargs...)

Enables mouse and keyboard controls for the objects of the live-rendered system `h`, see
[`live_render!`](@ref). Objects can be grabbed with the mouse and moved or rotated, while the
[`update_render!`](@ref) of `h` and the `on_change` callback are called after each change. The camera
can be used as usual as long as no object is grabbed. Returns a `KinematicController`.

# Mouse controls

- left-drag on an object: moves the object within the plane through its position (`plane_normal`)
- shift + left-drag on an object: rotates the object around the `rotation_axis`
- left-click on empty space: deselects the current object

# Keyboard controls

The following keys apply to the selected object, pressing shift multiplies the step size by 10:

- `↑`/`↓`: moves the object along its local y-axis by `fine_step`
- `←`/`→`: rotates the object around the `rotation_axis` by `fine_angle`
- `page up`/`page down`: tilts the object around its local x-axis by `fine_angle`
- `r`: resets the object to its initial pose
- `esc`: deselects the object

The selected object is marked by a box and three arrows: the direction of `↑` (green), the tilt
axis of `page up` (red) and the `rotation_axis` of `←` (blue). The key `h` shows or hides an overlay
of all controls.

# Keyword args

- `objects = nothing`: the movable objects, all objects of `h` by default
- `on_change = obj -> nothing`: called with the moved object after each change, e.g. to solve the system
- `plane_normal = [0, 0, 1]`: normal of the plane for mouse translation
- `rotation_axis = [0, 0, 1]`: rotation axis for mouse and keyboard rotation
- `rotate_speed = deg2rad(0.5)`: mouse rotation angle per pixel [rad]
- `fine_step = 10e-9`: keyboard translation step [m]
- `fine_angle = 10e-6`: keyboard rotation step [rad]
- `throttle = true`: limits updates to one per frame
- `show_help = false`: shows the controls overlay initially, otherwise only a hint
- `pick = ax -> Makie.pick(Makie.get_scene(ax))`: picking function
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
        throttle::Bool = true,
        show_help::Bool = false,
        pick = _default_pick
    )
    scene = Makie.get_scene(ax)
    movable = isnothing(objects) ? BMO.AbstractObject[oh.obj for oh in h.handles] :
              collect(BMO.AbstractObject, objects)
    init_poses = IdDict{BMO.AbstractObject, Tuple{Point3{Float64}, Matrix{Float64}}}(
        obj => _pose(obj) for obj in movable)
    box_obs = Observable(Point3f[])
    box_plot = linesegments!(ax, box_obs; color = :yellow, linewidth = 2)
    axes_pos = Observable(fill(Point3f(0), 3))
    axes_dir = Observable(fill(Vec3f(0, 0, 1), 3))
    axes_visible = Observable(false)
    axes_colors = [:green, :red, :blue]
    axes_plot = arrows3d!(ax, axes_pos, axes_dir; color = axes_colors, visible = axes_visible,
        shaftradius = 0.015, tipradius = 0.04, tiplength = 0.12, overdraw = true)
    labels_plot = text!(ax, Makie.lift((p, d) -> p .+ 1.15 .* d, axes_pos, axes_dir);
        text = ["↑", "page up", "←"], color = axes_colors, visible = axes_visible,
        fontsize = 16, align = (:center, :center), overdraw = true)
    help = _help_text(fine_step, fine_angle)
    help_obs = Observable(show_help ? help : _HELP_HINT)
    help_plot = text!(ax, Point2f(0.01, 0.99); text = help_obs, space = :relative,
        align = (:left, :top), fontsize = 14, color = :gray40)

    ctrl = KinematicController(
        ax, h, movable, init_poses, Observable{Union{Nothing, BMO.AbstractObject}}(nothing),
        on_change, normalize(Float64.(plane_normal)), normalize(Float64.(rotation_axis)),
        Float64(rotate_speed), Float64(fine_step), Float64(fine_angle), throttle,
        false, false, false, zeros(3), zeros(3), (0.0, 0.0), nothing,
        box_obs, axes_pos, axes_dir, axes_visible, help_obs,
        AbstractPlot[box_plot, axes_plot, labels_plot, help_plot], Any[], nothing
    )

    # High priority, so that the camera does not receive events while an object is grabbed
    l1 = on(events(scene).mousebutton, priority = 200) do event
        event.button == Mouse.left || return Consume(false)
        if event.action == Mouse.press
            plot, _ = pick(ax)
            obj = isnothing(plot) ? nothing : pick_object(h, plot)
            if isnothing(obj) || !_is_movable(ctrl, obj)
                ctrl.bg_press_pos = _px(scene)
                return Consume(false)
            end
            ctrl.selected[] = obj
            if _shift_pressed(scene)
                ctrl.rotating = true
                ctrl.last_mouse = _px(scene)
            else
                # Keep the offset between object and mouse to avoid a jump at the start
                ctrl.dragging = true
                ctrl.plane_point = Vector{Float64}(position(obj))
                hit = _mouse_plane_hit(scene, ctrl)
                ctrl.grab_offset = isnothing(hit) ? zeros(3) : ctrl.plane_point .- hit
            end
            _update_selection_box!(ctrl)
            return Consume(true)
        elseif event.action == Mouse.release
            if ctrl.dragging || ctrl.rotating
                ctrl.dragging = ctrl.rotating = false
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
        isnothing(obj) && return Consume(false)
        if ctrl.dragging
            hit = _mouse_plane_hit(scene, ctrl)
            if !isnothing(hit)
                translate_to3d!(obj, hit .+ ctrl.grab_offset)
                _request_update!(ctrl)
            end
            return Consume(true)
        elseif ctrl.rotating
            mp = _px(scene)
            dx = mp[1] - ctrl.last_mouse[1]
            ctrl.last_mouse = mp
            if dx != 0
                rotate3d!(obj, ctrl.rotation_axis, ctrl.rotate_speed * dx)
                _request_update!(ctrl)
            end
            return Consume(true)
        end
        return Consume(false)
    end

    l3 = on(events(scene).keyboardbutton, priority = 200) do event
        event.action in (Keyboard.press, Keyboard.repeat) || return Consume(false)
        if event.key == Keyboard.h && event.action == Keyboard.press
            help_obs[] = help_obs[] == _HELP_HINT ? help : _HELP_HINT
            return Consume(true)
        end
        obj = ctrl.selected[]
        isnothing(obj) && return Consume(false)
        if event.key == Keyboard.escape
            ctrl.selected[] = nothing
            _update_selection_box!(ctrl)
            return Consume(true)
        end
        step = ctrl.fine_step * (_shift_pressed(scene) ? 10 : 1)
        angle = ctrl.fine_angle * (_shift_pressed(scene) ? 10 : 1)
        if event.key == Keyboard.up
            translate3d!(obj, step .* orientation(obj)[:, 2])
        elseif event.key == Keyboard.down
            translate3d!(obj, -step .* orientation(obj)[:, 2])
        elseif event.key == Keyboard.left
            rotate3d!(obj, ctrl.rotation_axis, angle)
        elseif event.key == Keyboard.right
            rotate3d!(obj, ctrl.rotation_axis, -angle)
        elseif event.key == Keyboard.page_up
            rotate3d!(obj, orientation(obj)[:, 1], angle)
        elseif event.key == Keyboard.page_down
            rotate3d!(obj, orientation(obj)[:, 1], -angle)
        elseif event.key == Keyboard.r
            _reset_pose!(ctrl, obj)
        else
            return Consume(false)
        end
        _request_update!(ctrl)
        return Consume(true)
    end
    push!(ctrl.listeners, l1, l2, l3)

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
