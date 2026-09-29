#=
Camera tools
=#

const _CAMERA_HELP = "g: zoom to selection, click on a beam: inspect it, esc: cancel a long trace"

"""
    _CameraAnimation

State of the animated transition of the camera of a `LiveView` to a new view, see
`_animate_camera!`. Like `_CubeAnimation`, the direction from `lookat` to the eye is rotated by
`angle` about `axis` and the up vector is rolled by `roll`, while `lookat` and the distance `dist`
of the eye are interpolated linearly.
"""
mutable struct _CameraAnimation
    lookat0::Vector{Float64}
    dist0::Float64
    o0::Vector{Float64}
    u0::Vector{Float64}
    # target view, applied as given at the end
    eye1::Vector{Float64}
    lookat1::Vector{Float64}
    up1::Vector{Float64}
    dist1::Float64
    axis::Vector{Float64}
    angle::Float64
    roll::Float64
    duration::Float64
    elapsed::Float64
end

"""Returns the current view `(eye, lookat, up)` of the 3D view of the `gui`."""
function _current_view(gui::LiveView)
    cam = cameracontrols(gui.ax.scene)
    return (Vector{Float64}(cam.eyeposition[]), Vector{Float64}(cam.lookat[]),
        Vector{Float64}(cam.upvector[]))
end

"""Applies the frame at the fraction `t` of the animation `a` to the camera of the `gui`."""
function _apply_camera_frame!(gui::LiveView, a::_CameraAnimation, t)
    if t >= 1
        set_view(gui.ax, a.eye1, a.lookat1, a.up1)
        return nothing
    end
    # Smooth start and stop
    τ = t^2 * (3 - 2t)
    o = normalize(_rotate(a.o0, a.axis, τ * a.angle))
    up = normalize(_rotate(_rotate(a.u0, a.axis, τ * a.angle), o, τ * a.roll))
    lookat = a.lookat0 .+ τ .* (a.lookat1 .- a.lookat0)
    dist = a.dist0 + τ * (a.dist1 - a.dist0)
    set_view(gui.ax, lookat .+ dist .* o, lookat, up)
    return nothing
end

"""Advances the camera animation of the `gui` by `dt` [s]. A new animation of the view cube wins."""
function _step_camera!(gui::LiveView, dt)
    a = gui.camera.animation
    isnothing(a) && return nothing
    if !isnothing(gui.widgets.view_cube) && !isnothing(gui.widgets.view_cube.anim)
        gui.camera.animation = nothing
        return nothing
    end
    a.elapsed += dt
    t = a.duration > 0 ? min(a.elapsed / a.duration, 1.0) : 1.0
    t >= 1 && (gui.camera.animation = nothing)
    _apply_camera_frame!(gui, a, t)
    return nothing
end

"""
    _animate_camera!(gui, eye, lookat, up)

Moves the camera of the `gui` to the view `eye`, `lookat`, `up`, animated over the duration of the
view cube (0.3 s without a view cube), driven by `tick` events.
"""
function _animate_camera!(gui::LiveView, eye, lookat, up)
    eye0, lookat0, up0 = _current_view(gui)
    d0, d1 = eye0 .- lookat0, Vector{Float64}(eye) .- Vector{Float64}(lookat)
    dist0, dist1 = norm(d0), norm(d1)
    o0 = dist0 > 0 ? d0 ./ dist0 : [0.0, 0.0, 1.0]
    o1 = dist1 > 0 ? d1 ./ dist1 : o0
    u0 = something(_orthogonalize(up0, o0), _perp(o0))
    u1 = something(_orthogonalize(Vector{Float64}(up), o1), _perp(o1))
    axis, angle = _rotation_between(o0, o1, cross(o0, u0))
    ur = _rotate(u0, axis, angle)
    roll = atan(dot(cross(ur, u1), o1), dot(ur, u1))
    duration = isnothing(gui.widgets.view_cube) ? 0.3 : gui.widgets.view_cube.duration
    isnothing(gui.widgets.view_cube) || (gui.widgets.view_cube.anim = nothing)
    gui.camera.animation = _CameraAnimation(lookat0, dist0, o0, u0, Vector{Float64}(eye),
        Vector{Float64}(lookat), Vector{Float64}(up), dist1, axis, angle, roll, duration, 0.0)
    duration > 0 || _step_camera!(gui, 0.0)
    return nothing
end

"""
    _zoom_box(gui)

Returns the bounding box of the plots of the selected object of the `gui`, or of all systems and
the visible extras if nothing is selected, or `nothing` if there are no visible plots.
"""
function _zoom_box(gui::LiveView)
    ctrl = gui.controls
    obj = ctrl.selected[]
    isnothing(obj) || return _selection_bbox(ctrl, obj, _object_plots(ctrl.h, obj))
    return _visible_bbox((gui.system_handles..., gui.extras))
end

"""
    _visible_bbox(handles) -> Union{Rect3d, Nothing}

Returns the bounding box of the visible plots of the objects of the `handles` (`SystemRenderHandle`s,
e.g. of the systems and the extras of a live view), or `nothing` if there are none.
"""
function _visible_bbox(handles)
    bbs = [Makie.boundingbox(p) for h in handles for oh in h.handles for p in oh.plots if p.visible[]]
    filter!(_is_finite_box, bbs)
    return isempty(bbs) ? nothing : reduce(GeometryBasics.union, bbs)
end

"""
Size of the scene of the `handles` (see `_visible_bbox`): the largest edge of its bounding box, or
0.125 m for an empty scene. Sets the size of the source markers and of new clip planes.
"""
function _scene_extent(handles)
    bb = _visible_bbox(handles)
    return isnothing(bb) ? 0.125 : Float64(maximum(GeometryBasics.widths(bb)))
end

"""
    _zoom_to_selection!(gui)

Moves the camera of the `gui` such that the bounding sphere of the selected object, or of all
systems and the visible extras if nothing is selected (see `_zoom_box`), fills the view: `lookat`
is set to the center of the bounding box, the eye to the distance `(d/2) / sin(fov/2)`
(orthographic: `d/2`) along the current view direction, where `d` is the diagonal of the box.
"""
function _zoom_to_selection!(gui::LiveView)
    bb = _zoom_box(gui)
    if isnothing(bb)
        gui.status.text[] = "nothing to zoom to"
        return nothing
    end
    center = Vector{Float64}(minimum(bb) .+ GeometryBasics.widths(bb) ./ 2)
    d = max(norm(Vector{Float64}(GeometryBasics.widths(bb))), 1e-6)
    cam = cameracontrols(gui.ax.scene)
    # The orthographic view is scaled by the distance of the eye, which is its half height
    dist = cam.settings.projectiontype[] == Makie.Perspective ? (d / 2) / sind(cam.fov[] / 2) : d / 2
    eye, lookat, up = _current_view(gui)
    o = norm(eye .- lookat) > 0 ? normalize(eye .- lookat) : [0.0, 0.0, 1.0]
    _animate_camera!(gui, center .+ dist .* o, center, up)
    return nothing
end

"""Validates the `views` kwarg of `live_view` and returns a vector of `name => (eye, lookat, up)`."""
function _view_specs(views)
    specs = Pair{String, NTuple{3, Vector{Float64}}}[]
    for v in views
        ok = v isa Pair && v.second isa Tuple && length(v.second) == 3 &&
             all(x -> x isa AbstractVector{<:Real} && length(x) == 3, v.second)
        ok || throw(ArgumentError("invalid view $v, use \"name\" => (eye, lookat, up)"))
        push!(specs, string(v.first) => Tuple(Vector{Float64}.(v.second)))
    end
    return specs
end

"""Returns the options of the views menu, the option `i` sets `views[i]`."""
function _views_options(views)
    isempty(views) && return [("no views", 0)]
    return [(name, i) for (i, (name, _)) in enumerate(views)]
end

"""Sets the camera of the `gui` to the saved view `i`, see `views`."""
function _set_saved_view!(gui::LiveView, i)
    1 <= i <= length(gui.camera.views) || return nothing
    name, (eye, lookat, up) = gui.camera.views[i]
    _animate_camera!(gui, eye, lookat, up)
    gui.status.text[] = "view \"$name\""
    return nothing
end

"""
    _save_view!(gui; io = stdout)

Appends the current view of the `gui` as `"view n"` to the saved views and prints it as the Julia
code of an entry of the `views` kwarg of `live_view` to `io`. Returns the code.
"""
function _save_view!(gui::LiveView; io::IO = stdout)
    n = length(gui.camera.views) + 1
    names = Set(first.(gui.camera.views))
    while "view $n" in names
        n += 1
    end
    name = "view $n"
    view = _current_view(gui)
    push!(gui.camera.views, name => view)
    gui.widgets.views_menu.options[] = _views_options(gui.camera.views)
    code = "$(repr(name)) => ($(join(_vector_code.(view), ", ")))"
    println(io, code)
    gui.status.text[] = "saved $(repr(name)), printed as code"
    return code
end

"""Restores the home view of the `gui`, i.e. the view when the window was shown."""
function _go_home!(gui::LiveView)
    _animate_camera!(gui, gui.camera.home...)
    gui.status.text[] = "home view"
    return nothing
end

"""
Connects the camera tools of the `gui`: the key `g` (zoom to the selection), the home button, the
views menu and the save view button, and the animation and the home view via `tick`.
"""
function _connect_camera!(gui::LiveView)
    listeners = gui.controls.listeners
    scene = gui.ax.scene
    push!(listeners, on(events(scene).keyboardbutton, priority = 200) do event
        (event.action == Keyboard.press && event.key == Keyboard.g) || return Consume(false)
        gui.controls.ignore_keys() && return Consume(false)
        _zoom_to_selection!(gui)
        return Consume(true)
    end)
    push!(listeners, on(events(scene).tick) do tick
        # The home view is the view when the window is shown, e.g. after `set_view`
        if !gui.camera.home_set
            gui.camera.home = _current_view(gui)
            gui.camera.home_set = true
        end
        _step_camera!(gui, tick.delta_time)
        return nothing
    end)
    # The textboxes of the cards are connected with their cards, see `_connect_card!`, those of the
    # layout (see `_layout_boxes`) by the layout
    for obs in (cameracontrols(scene).selected, _open_observable(gui.widgets.menu), gui.widgets.views_menu.is_open)
        _listen!(listeners, _ -> _keep_keyboard!(gui), obs)
    end
    _keep_keyboard!(gui)
    push!(listeners, on(_ -> _go_home!(gui), gui.widgets.home_button.clicks))
    push!(listeners, on(i -> _set_saved_view!(gui, something(i, 0)), gui.widgets.views_menu.i_selected))
    push!(listeners, on(_ -> _save_view!(gui), gui.widgets.save_view_button.clicks))
    gui.controls.help_extra *= "\n" * _CAMERA_HELP
    _update_help!(gui.controls)
    return nothing
end

"""
    _connect_projection!(gui, orthographic)

Connects the orthographic toggle of the `gui` and sets the initial projection.

The orthographic `Camera3D` scales the view with the distance between eye and `lookat` and clips
the scene at `near`/`far` in front of the eye. Zooming in, or moving along the view direction, then
moves the eye into the scene: the near plane cuts through it like a clip plane, while the view
itself stays the same. In orthographic mode the depth range is therefore static and reaches from
behind the eye (negative `near`) past the far side of the scene, updated whenever the eye moves.
The scene bounds are taken when switching to orthographic, with a margin of one radius for later
moves of objects. Makie's `center!`, called e.g. by `reset_limits!` when a plot is added to the
open window or by `save` without `update = false`, replaces the depth range with a positive one;
it is restored, with new scene bounds, before the camera applies it.

Moving along the view direction does not change an orthographic view at all, so the keys for it
(`W`/`S` of `Camera3D`) zoom there, like the zoom keys (`U`/`O`).
"""
function _connect_projection!(gui::LiveView, orthographic::Bool)
    cam = cameracontrols(gui.ax.scene)
    settings, controls = cam.settings, cam.controls
    perspective_depth = (settings.clipping_mode[], cam.near[], cam.far[])
    walk_keys = (controls.forward_key[], controls.backward_key[], controls.zoom_in_key[],
        controls.zoom_out_key[])
    bounds = Ref((zeros(3), 1.0))
    is_ortho() = settings.projectiontype[] != Makie.Perspective
    # Silent updates: every camera move calls `update_cam!` after setting the eye position
    function set_depth!()
        center, radius = bounds[]
        d = norm(Vector{Float64}(cam.eyeposition[]) .- center) + 2 * radius
        cam.near.val = -d
        cam.far.val = d
        return nothing
    end
    # `false` is a key binding that is never pressed
    function set_keys!(ortho::Bool)
        forward, backward, zoom_in, zoom_out = walk_keys
        controls.forward_key[] = ortho ? false : forward
        controls.backward_key[] = ortho ? false : backward
        controls.zoom_in_key[] = ortho ? zoom_in | forward : zoom_in
        controls.zoom_out_key[] = ortho ? zoom_out | backward : zoom_out
        return nothing
    end
    function set_projection!(ortho::Bool)
        if ortho
            bounds[] = _scene_bounds(gui)
            settings.clipping_mode.val = :static
            set_depth!()
        else
            mode, near, far = perspective_depth
            settings.clipping_mode.val = mode
            cam.near.val = near
            cam.far.val = far
        end
        set_keys!(ortho)
        settings.projectiontype[] = ortho ? Makie.Orthographic : Makie.Perspective
        return nothing
    end
    listeners = gui.controls.listeners
    push!(listeners, on(set_projection!, gui.widgets.orthographic_toggle.active))
    push!(listeners, on(cam.eyeposition) do _
        is_ortho() && set_depth!()
        return nothing
    end)
    # Before the listener of the camera, which applies `near` and `far`
    for depth in (cam.near, cam.far)
        push!(listeners, on(depth; priority = 1) do _
            if is_ortho()
                bounds[] = _scene_bounds(gui)
                set_depth!()
            end
            return nothing
        end)
    end
    set_projection!(orthographic)
    return nothing
end

"""Returns the center and radius of the bounding sphere of the plots of the 3D view of the `gui`."""
function _scene_bounds(gui::LiveView)
    bb = Makie.data_limits(gui.ax.scene)
    _is_finite_box(bb) || return (zeros(3), 1.0)
    w = Vector{Float64}(GeometryBasics.widths(bb))
    return Vector{Float64}(minimum(bb)) .+ w ./ 2, max(norm(w) / 2, 1e-6)
end
