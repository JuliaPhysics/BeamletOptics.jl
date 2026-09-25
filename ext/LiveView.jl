using Makie: Figure, Axis, Label, SliderGrid, GridLayout, DataAspect, Relative, colsize!, colgap!,
             heatmap!, autolimits!, limits!, Button, Toggle, Textbox, Menu
import InteractiveUtils

"""
    DetectorPanel

Detector panel of a `LiveView`, which shows the spot diagram or the intensity of a
`Detector` in a 2D `Axis`, see `live_view`.
"""
mutable struct DetectorPanel
    pd::BMO.Detector
    name::String
    # :auto, :spot or :intensity
    mode::Symbol
    kwargs::NamedTuple
    ax::Axis
    xy::Observable{Vector{Point2f}}
    heat_x::Observable{Vector{Float32}}
    heat_y::Observable{Vector{Float32}}
    heat_I::Observable{Matrix{Float32}}
    scatter_plot::AbstractPlot
    heat_plot::AbstractPlot
    # last error of the update, logged only once
    last_error::Union{Nothing, String}
end

Base.show(io::IO, p::DetectorPanel) = print(io, "DetectorPanel(", p.name, ", mode = ", p.mode, ")")

"""
    LiveView

Interactive window returned by [`live_view`](@ref). The `Figure` is stored in `fig`, the `LScene`
of the 3D view in `ax`, the `KinematicController` in `controls` and the view cube in `view_cube`
(or `nothing`). Use `display` to show the window and `close` to remove the controls and the view
cube.

The systems are solved after each change if `auto_trace[]` is `true`, otherwise only via the
`trace_button`, the key `t` or by switching the `auto_trace_toggle` on. `stale` is `true` if the
beams and detector panels do not match the current poses of the objects. If solving takes longer
than `trace_budget` [s], the systems are solved once the movement pauses for `idle_delay` [s].
The clip planes of the 3D view are stored in `clip_planes`, which are applied if `clipping` is
`true`. The `orthographic_toggle` switches the 3D view between perspective and orthographic
projection.

The `export_button` prints the changed poses as Julia code, see [`export_changes`](@ref), and
copies them to the clipboard if `export_clipboard` is `true`. The component `menu` lists the
movable objects `menu_objects`, the objects hidden via the `hide_button` are stored in `hidden`.
The `pose_boxes` of the pose inspector show and set the position `x`, `y`, `z` [mm] of the
selected object and rotate it about the red, green and blue axes of the controls [mrad].
"""
mutable struct LiveView
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
    last_error::Union{Nothing, String}
    # manual tracing, auto_trace is the `active` observable of the toggle
    auto_trace::Observable{Bool}
    stale::Bool
    trace_button::Button
    auto_trace_toggle::Toggle
    # alpha of the beam plots before dimming, restored after tracing
    beam_alphas::IdDict{Any, Any}
    # adaptive tracing, durations of the last solve and panel update [s]
    trace_budget::Float64
    idle_delay::Float64
    solve_time::Float64
    panel_time::Float64
    # a solve is deferred until the movement pauses
    pending::Bool
    pending_obj::Any
    last_change::Float64
    # the panels show a preview on a coarse grid
    coarse::Bool
    labels::IdDict{Any, String}
    step_box::Textbox
    # clip planes, switched on and off via `clipping`
    clip_planes::Vector{LiveClipPlane}
    clipping::Bool
    # edge length of the outline of new clip planes, fixed at construction, since the bounding
    # boxes of clipped plots only cover their visible part
    clip_size::Float64
    # the beams are clipped as well, switched via the toggle
    clip_beams::Bool
    clip_beams_toggle::Toggle
    # view cube in the corner of the 3D view, see `view_cube!`
    view_cube::Union{Nothing, ViewCube}
    # orthographic projection of the 3D view, switched via the toggle
    orthographic_toggle::Toggle
    # export of the changed poses, see `export_changes`
    export_button::Button
    export_clipboard::Bool
    # component menu, the option `i` selects `menu_objects[i]`
    menu::Menu
    menu_objects::Vector{Any}
    # hidden objects, i.e. rendered objects (leaves of groups) whose plots are invisible
    hide_button::Button
    show_all_button::Button
    hidden::Base.IdSet{Any}
    # pose inspector: x, y, z [mm] and the rotations rx, ry, rv [mrad], see `_POSE_FIELDS`
    pose_boxes::Vector{Textbox}
end

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
    close(gui.controls)
    isnothing(gui.view_cube) || close(gui.view_cube)
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

"""Returns all `Detector`s of the `systems`, deduplicated by identity, in the order of discovery."""
function _find_detectors(systems)
    seen = Base.IdSet{Any}()
    dets = BMO.Detector[]
    for sys in systems, obj in BMO.objects(sys)
        obj isa BMO.Detector || continue
        obj in seen && continue
        push!(seen, obj)
        push!(dets, obj)
    end
    return dets
end

const _PANEL_MODES = (:auto, :spot, :intensity)

_panel_spec(pd::BMO.Detector) = (pd, :auto, (;))
_panel_spec(p::Pair{<:BMO.Detector, Symbol}) = (p.first, p.second, (;))
_panel_spec(p::Pair{<:BMO.Detector, <:Tuple{Symbol, NamedTuple}}) = (p.first, p.second...)
function _panel_spec(x)
    throw(ArgumentError("invalid detector panel $x, use `pd`, `pd => mode` or `pd => (mode, kwargs)`"))
end

"""Returns the `(detector, mode, kwargs)` of all panels of the `detectors` kwarg of `live_view`."""
function _panel_specs(detectors, systems)
    if detectors isa Symbol
        detectors == :auto || throw(ArgumentError("detectors must be :auto or a vector, got :$detectors"))
        return [(pd, :auto, (;)) for pd in _find_detectors(systems)]
    end
    specs = [_panel_spec(d) for d in detectors]
    for (_, mode, _) in specs
        mode in _PANEL_MODES || throw(ArgumentError("detector panel mode must be one of $_PANEL_MODES, got :$mode"))
    end
    return specs
end

function DetectorPanel(parent, pd::BMO.Detector, name::String, mode::Symbol, kwargs::NamedTuple)
    ax = Axis(parent; title = "$name: no hits", xlabel = "x [mm]", ylabel = "z [mm]",
        aspect = DataAspect())
    n = get(kwargs, :n, 100)
    xy = Observable(Point2f[])
    heat_x = Observable(zeros(Float32, n))
    heat_y = Observable(zeros(Float32, n))
    heat_I = Observable(zeros(Float32, n, n))
    heat_plot = heatmap!(ax, heat_x, heat_y, heat_I; colorrange = (0.0f0, 1.0f0), visible = false)
    scatter_plot = scatter!(ax, xy; markersize = 3, color = :black, visible = false)
    return DetectorPanel(pd, name, mode, kwargs, ax, xy, heat_x, heat_y, heat_I, scatter_plot,
        heat_plot, nothing)
end

"""Formats `x` with 3 decimal places."""
function _fmt3(x::Real)
    s = string(round(x, digits = 3))
    occursin(r"[eE]", s) && return s
    i = findfirst('.', s)
    isnothing(i) && return s * ".000"
    return s * "0"^max(0, 3 - (length(s) - i))
end

_is_beamlet_hits(h) = h isa AbstractVector{<:BMO.AbstractBeamletHit}

function _resolve_mode(mode::Symbol, h)
    mode == :auto || return mode
    return _is_beamlet_hits(h) ? :intensity : :spot
end

function _clear_panel!(p::DetectorPanel)
    p.scatter_plot.visible[] && (p.scatter_plot.visible[] = false)
    p.heat_plot.visible[] && (p.heat_plot.visible[] = false)
    isempty(p.xy[]) || (p.xy[] = Point2f[])
    return nothing
end

function _hits_title(p::DetectorPanel, h)
    n = length(h)
    return _is_beamlet_hits(h) ? "$(p.name): $n beamlets" : "$(p.name): $n rays"
end

function _update_spot!(p::DetectorPanel, h)
    pts = BMO.spot_diagram(p.pd)
    p.xy[] = [Point2f(1e3 * q[1], 1e3 * q[2]) for q in pts]
    p.heat_plot.visible[] && (p.heat_plot.visible[] = false)
    p.scatter_plot.visible[] || (p.scatter_plot.visible[] = true)
    p.ax.title[] = _hits_title(p, h)
    autolimits!(p.ax)
    _pad_degenerate_limits!(p.ax, p.xy[])
    return nothing
end

"""
    _pad_degenerate_limits!(ax, xy; ε = 1e-6)

Sets the limits of the spot diagram `ax` if the x or z extent of the spots `xy` [mm] is below `ε`,
e.g. for a single ray, since `DataAspect` gives degenerate limits for a zero-width data range. A
degenerate axis gets half the extent of the other axis around its center, or 1 µm if both are
degenerate. Otherwise the limits of `autolimits!` are kept.
"""
function _pad_degenerate_limits!(ax::Axis, xy; ε = 1e-6)
    isempty(xy) && return nothing
    xmin, xmax = extrema(q -> Float64(q[1]), xy)
    zmin, zmax = extrema(q -> Float64(q[2]), xy)
    Δx, Δz = xmax - xmin, zmax - zmin
    (Δx < ε || Δz < ε) || return nothing
    h = max(Δx, Δz) / 2
    h = h < ε ? 1e-3 : h
    hx = Δx < ε ? h : Δx / 2 * 1.05
    hz = Δz < ε ? h : Δz / 2 * 1.05
    cx, cz = (xmin + xmax) / 2, (zmin + zmax) / 2
    limits!(ax, cx - hx, cx + hx, cz - hz, cz + hz)
    return nothing
end

"""Grid size of the intensity of the panel `p`, reduced to a preview while moving objects."""
function _panel_n(p::DetectorPanel, coarse::Bool)
    n = get(p.kwargs, :n, 100)
    return coarse ? max(16, n ÷ 4) : n
end

function _update_intensity!(p::DetectorPanel, h; coarse = false)
    x, z, I = BMO.intensity(p.pd; merge(p.kwargs, (; n = _panel_n(p, coarse)))...)
    # The plot is updated lazily, hence the grid size may change
    p.heat_x[] = Float32.(1e3 .* x)
    p.heat_y[] = Float32.(1e3 .* z)
    p.heat_I[] = Float32.(I)
    Imax = Float32(maximum(I))
    p.heat_plot.colorrange[] = (0.0f0, Imax > 0 ? Imax : 1.0f0)
    p.scatter_plot.visible[] && (p.scatter_plot.visible[] = false)
    p.heat_plot.visible[] || (p.heat_plot.visible[] = true)
    if h isa AbstractVector{<:BMO.GaussianBeamletHit}
        # Optical power from the computed intensity like optical_power, avoids a second field evaluation
        P = BMO.trapz((x, z), I)
        p.ax.title[] = "$(p.name): P = $(_fmt3(1e3 * P)) mW"
    else
        p.ax.title[] = _hits_title(p, h)
    end
    coarse && (p.ax.title[] *= " (preview)")
    autolimits!(p.ax)
    return nothing
end

"""Updates the plots and the title of the panel `p` after the systems have been solved."""
function _update_panel!(p::DetectorPanel; coarse = false)
    try
        h = BMO.hits(p.pd)
        if isnothing(h)
            _clear_panel!(p)
            p.ax.title[] = "$(p.name): no hits"
        else
            mode = _resolve_mode(p.mode, h)
            if mode == :intensity && p.mode == :auto && h isa AbstractVector{<:BMO.AstigmaticGaussianBeamletHit}
                # Fall back to the spot diagram if the field of the hits can not be evaluated
                try
                    _update_intensity!(p, h; coarse)
                catch
                    _update_spot!(p, h)
                end
            elseif mode == :intensity
                _update_intensity!(p, h; coarse)
            else
                _update_spot!(p, h)
            end
        end
        p.last_error = nothing
    catch e
        p.last_error = _log_once(e, p.last_error, "update of the panel \"$(p.name)\"")
        p.ax.title[] = "$(p.name): error"
    end
    return nothing
end

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

_label(gui, obj) = get(gui.labels, obj, string(nameof(typeof(obj))))

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
    _resolve!(gui::LiveView, obj)

Empties all `Detector`s, solves all systems and updates the beams, detector panels and the status
line of the `gui`. The user `on_change` is called with the moved `obj`, or `nothing`.
"""
function _resolve!(gui::LiveView, obj; coarse = false)
    # A detector can be part of several systems, hence empty all before solving
    # Monotonic clock with ns resolution, time() is too coarse on Windows for fast solves
    t0 = time_ns()
    foreach(empty!, _find_detectors(first.(gui.pairs)))
    for (i, (sys, beam)) in enumerate(gui.pairs)
        solve_system!(sys, beam)
        update_render!(gui.beam_handles[i])
    end
    t1 = time_ns()
    foreach(p -> _update_panel!(p; coarse), gui.panels)
    gui.solve_time = 1e-9 * (t1 - t0)
    coarse || (gui.panel_time = 1e-9 * (time_ns() - t1))
    gui.coarse = coarse
    try
        gui.on_change(gui, obj)
        gui.last_error = nothing
    catch e
        gui.last_error = _log_once(e, gui.last_error, "`on_change` callback")
    end
    isnothing(obj) || (gui.status.text[] = _pose_string(gui, obj))
    return nothing
end

const _STALE_ALPHA = 0.3

_beam_plots(h::BeamRenderHandle) = AbstractPlot[h.plot]
_beam_plots(h::GaussianRenderHandle) = AbstractPlot[h.mesh_plot; h.beam_plots]
_beam_plots(h::AstigmaticGroupRenderHandle) = AbstractPlot[h.mesh_plot]
_beam_plots(h) = AbstractPlot[]

"""Dims all beam plots of the `gui` to indicate outdated beams, stores the original `alpha`."""
function _dim_beams!(gui::LiveView)
    for h in gui.beam_handles, plot in _beam_plots(h)
        haskey(plot, :alpha) || continue
        haskey(gui.beam_alphas, plot) || (gui.beam_alphas[plot] = plot.alpha[])
        plot.alpha[] = _STALE_ALPHA
    end
    return nothing
end

"""Restores the `alpha` of all beam plots of the `gui` after dimming."""
function _restore_beams!(gui::LiveView)
    for (plot, alpha) in gui.beam_alphas
        plot.alpha[] = alpha
    end
    empty!(gui.beam_alphas)
    return nothing
end

"""Marks the beams and detector panels of the `gui` as outdated after `obj` (or a slider) changed."""
function _mark_stale!(gui::LiveView, obj; msg = "outdated, press t to trace")
    gui.stale || _dim_beams!(gui)
    gui.stale = true
    gui.status.text[] = isnothing(obj) ? msg : "$(_pose_string(gui, obj)) — $msg"
    return nothing
end

"""
    _solve!(gui::LiveView, obj)

Solves all systems of the `gui` via `_resolve!` and restores the appearance of the beams. If solving
fails, the beams and detector panels are kept marked as outdated. Returns `true` on success.
"""
function _solve!(gui::LiveView, obj; coarse = false)
    gui.pending = false
    try
        _resolve!(gui, obj; coarse)
    catch e
        gui.last_error = _log_once(e, gui.last_error, "solving the systems")
        gui.stale || _dim_beams!(gui)
        gui.stale = true
        gui.status.text[] = "solving the systems failed, see the log"
        return false
    end
    _restore_beams!(gui)
    gui.stale = false
    return true
end

"""
    _on_change!(gui::LiveView, obj)

Called after each change of `obj` (or a slider, then `obj` is `nothing`). Solves the systems, with
a coarse preview of slow detector panels, or marks them as outdated. If solving is slower than the
`trace_budget`, the solve is deferred until the movement pauses, see `_on_idle!`.
"""
function _on_change!(gui::LiveView, obj)
    gui.last_change = time()
    if !gui.auto_trace[]
        _mark_stale!(gui, obj)
    elseif gui.solve_time <= gui.trace_budget
        _solve!(gui, obj; coarse = gui.panel_time > gui.trace_budget)
    else
        _mark_stale!(gui, obj; msg = "tracing when the movement pauses")
        gui.pending = true
        gui.pending_obj = obj
    end
    return nothing
end

"""Solves deferred changes and refines the preview of the detector panels once the movement pauses."""
function _on_idle!(gui::LiveView)
    time() - gui.last_change > gui.idle_delay || return nothing
    if gui.pending && gui.auto_trace[]
        _solve!(gui, gui.pending_obj)
    elseif gui.coarse
        foreach(_update_panel!, gui.panels)
        gui.coarse = false
    end
    return nothing
end

"""Solves all systems of the `gui` on request, with the currently selected object."""
function _trace!(gui::LiveView)
    obj = gui.controls.selected[]
    # A clip plane is not part of the systems
    obj isa LiveClipPlane && (obj = nothing)
    _solve!(gui, obj) && isnothing(obj) && (gui.status.text[] = "traced")
    return nothing
end

"""Connects the trace button, the key `t` and the auto trace toggle of the `gui`."""
function _connect_trace!(gui::LiveView)
    listeners = gui.controls.listeners
    push!(listeners, on(_ -> _trace!(gui), gui.trace_button.clicks))
    push!(listeners, on(events(gui.ax.scene).keyboardbutton, priority = 200) do event
        (event.action == Keyboard.press && event.key == Keyboard.t) || return Consume(false)
        gui.controls.ignore_keys() && return Consume(false)
        _trace!(gui)
        return Consume(true)
    end)
    push!(listeners, on(gui.auto_trace) do active
        active && gui.stale && _trace!(gui)
        return nothing
    end)
    push!(listeners, on(_ -> _on_idle!(gui), events(gui.ax.scene).tick))
    return nothing
end

"""
    _parse_step(s)

Parses a step such as `"250 nm"`, `"0.1 mm"`, `"50 µrad"` or `"1 deg"` into `(:move, step [m])` or
`(:rotate, step [rad])`. Returns `nothing` if `s` is invalid.
"""
function _parse_step(s::AbstractString)
    m = match(r"^\s*([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*(\S+)\s*$", s)
    isnothing(m) && return nothing
    x = parse(Float64, m.captures[1])
    x > 0 || return nothing
    unit = replace(m.captures[2], "u" => "µ", "μ" => "µ")
    units = Dict("nm" => (:move, 1e-9), "µm" => (:move, 1e-6), "mm" => (:move, 1e-3),
        "cm" => (:move, 1e-2), "m" => (:move, 1.0), "µrad" => (:rotate, 1e-6),
        "mrad" => (:rotate, 1e-3), "rad" => (:rotate, 1.0), "deg" => (:rotate, deg2rad(1)),
        "°" => (:rotate, deg2rad(1)))
    haskey(units, unit) || return nothing
    mode, factor = units[unit]
    return mode, x * factor
end

"""Sets the keyboard step and the mode of the controls from the text `s` of the step textbox."""
function _set_step!(gui::LiveView, s)
    step = isnothing(s) ? nothing : _parse_step(s)
    if isnothing(step)
        gui.status.text[] = "invalid step \"$s\", use e.g. 250 nm or 50 µrad"
        return nothing
    end
    ctrl = gui.controls
    mode, x = step
    mode == :move ? (ctrl.fine_step = x) : (ctrl.fine_angle = x)
    if ctrl.mode[] != mode
        ctrl.mode[] = mode
        _update_selection_box!(ctrl)
    end
    _update_help!(ctrl)
    gui.status.text[] = "$mode step: $(_step_string(mode, ctrl.fine_step, ctrl.fine_angle))"
    return nothing
end

"""
    _connect_sliders!(gui, callbacks)

Calls the slider `callbacks` after a value change, then updates all systems and solves again, or
marks them as outdated if auto tracing is off. The updates are throttled to one per frame. The
listeners are added to the controls of the `gui`, such that `close(gui)` removes them.
"""
function _connect_sliders!(gui::LiveView, callbacks)
    pending = Dict{Int, Any}()
    errors = Vector{Union{Nothing, String}}(nothing, length(callbacks))
    listeners = gui.controls.listeners
    for (i, sl) in enumerate(gui.sliders.sliders)
        push!(listeners, on(v -> (pending[i] = v; nothing), sl.value))
    end
    tick = on(events(gui.ax.scene).tick) do _
        isempty(pending) && return nothing
        for i in sort!(collect(keys(pending)))
            v = pending[i]
            try
                callbacks[i](v)
                errors[i] = nothing
            catch e
                errors[i] = _log_once(e, errors[i], "slider callback")
            end
        end
        empty!(pending)
        # The callbacks may have moved objects or sources
        update_render!(gui.controls.h)
        _on_change!(gui, nothing)
        return nothing
    end
    push!(listeners, tick)
    return nothing
end

function _slider_spec(s::Pair)
    label, t = s
    (t isa Tuple && length(t) in (2, 3)) ||
        throw(ArgumentError("invalid slider $s, use \"label\" => (range, callback[, startvalue])"))
    range, callback = t[1], t[2]
    startvalue = length(t) == 3 ? t[3] : first(range)
    return (; label = string(label), range, startvalue), callback
end
_slider_spec(s) = throw(ArgumentError("invalid slider $s, use \"label\" => (range, callback[, startvalue])"))

# Makie ignores all clip planes of a plot beyond the 8th
const _MAX_CLIP_PLANES = 8

const _CLIP_HELP = "p: add clip plane, del: remove, c: clipping on/off, shift+c: flip"

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
    planes = gui.clipping ? Plane3f[_plane3f(p) for p in gui.clip_planes] : Plane3f[]
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
    beam_planes = gui.clip_beams ? planes : Plane3f[]
    for plot in beam_plots
        plot.clip_planes[] == beam_planes || (plot.clip_planes = beam_planes)
    end
    return nothing
end

"""Switches the clipping of the beams of the `gui` on or off, see `clip_beams`."""
function _set_clip_beams!(gui::LiveView, on::Bool)
    gui.clip_beams = on
    gui.clip_beams_toggle.active[] == on || (gui.clip_beams_toggle.active[] = on)
    _apply_clip_planes!(gui)
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
    plane = LiveClipPlane(point, normal, gui.clip_size)
    push!(ctrl.h.handles, _live_render_clip_plane!(gui.ax, plane))
    push!(ctrl.movable, plane)
    ctrl.init_poses[plane] = _pose(plane)
    push!(gui.clip_planes, plane)
    haskey(gui.labels, plane) || (gui.labels[plane] = "Clip plane")
    _apply_clip_planes!(gui)
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
    filter!(p -> p !== plane, gui.clip_planes)
    delete!(gui.labels, plane)
    _apply_clip_planes!(gui)
    gui.status.text[] = "clip plane removed"
    return nothing
end

"""Rotates the clip `plane` by π about its local x-axis, such that the other side is visible."""
function _flip_clip_plane!(gui::LiveView, plane::LiveClipPlane)
    ctrl = gui.controls
    P0, R0 = _pose(plane)
    rotate3d!(plane, plane.dir[:, 1], π)
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
        isempty(gui.clip_planes) && return false
        gui.clipping = !gui.clipping
        _apply_clip_planes!(gui)
        gui.status.text[] = gui.clipping ? "clipping on" : "clipping off"
        return true
    elseif key == Keyboard.p
        # Nothing can be selected in the spectator mode
        ctrl.spectator[] && return false
        if length(gui.clip_planes) >= _MAX_CLIP_PLANES
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
    gui.controls.help_extra = _CLIP_HELP
    _update_help!(gui.controls)
    return nothing
end

#=
Export of the changed poses
=#

"""
    _rotation_axis_angle(R)

Returns the `axis` and `angle` of the rotation matrix `R` such that `rotate3d(axis, angle) ≈ R`.
Unlike `_axis_angle_from_rotmatrix`, which is based on the trace of `R`, it is accurate to the
precision of `Float64` for all angles, including small ones.
"""
function _rotation_axis_angle(R::AbstractMatrix)
    # sin(angle) * axis
    v = [R[3, 2] - R[2, 3], R[1, 3] - R[3, 1], R[2, 1] - R[1, 2]] ./ 2
    s, c = norm(v), (tr(R) - 1) / 2
    angle = atan(s, c)
    iszero(s) && c > 0 && return [0.0, 0.0, 1.0], 0.0
    c > -0.5 && return v ./ s, angle
    # Close to π, where sin(angle) is inaccurate: (R + R') / 2 = c I + (1 - c) axis axis'
    B = (R + R') ./ 2 - c * I
    axis = normalize(B[:, argmax(diag(B))])
    dot(axis, v) < 0 && (axis = -axis)
    return axis, angle
end

"""Sets the pose of `obj` like `_set_pose!`, but with the accurate `_rotation_axis_angle`."""
function _set_pose_exact!(obj, P, R)
    axis, angle = _rotation_axis_angle(R * _pose(obj)[2]')
    angle > 0 && rotate3d!(obj, axis, angle)
    translate_to3d!(obj, P)
    return nothing
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

"""Returns `true` if `s` can be used as the name of a variable."""
_is_variable_name(s) = Base.isidentifier(s) && try
    Meta.parse(s) isa Symbol
catch
    false
end

"""
    _export_names(gui, objects)

Returns the variable names of the `objects` in the code of `export_changes`: the label of an
object if it is a valid, unique variable name, otherwise `obj1`, `obj2`, … by the position `k`
of the object in `objects`, i.e. in the component menu.
"""
function _export_names(gui::LiveView, objects)
    names = IdDict{Any, String}()
    used = Set{String}()
    for (k, obj) in enumerate(objects)
        label = get(gui.labels, obj, "")
        name = _is_variable_name(label) && !(label in used) ? label : "obj$k"
        while name in used
            name *= "_"
        end
        push!(used, name)
        names[obj] = name
    end
    return names
end

"""Formats the vector `v` as Julia code with full precision."""
_vector_code(v) = "[" * join((repr(Float64(x)) for x in v), ", ") * "]"

"""
    _initial_copy(ctrl, group)

Returns a copy of the `group` in which the group and all its objects are reset to their initial
poses of the controls `ctrl`, from the outside in.
"""
function _initial_copy(ctrl::KinematicController, group)
    copy = deepcopy(group)
    for (obj, c) in zip(_descendants(group), _descendants(copy))
        haskey(ctrl.init_poses, obj) && _set_pose_exact!(c, ctrl.init_poses[obj]...)
    end
    return copy
end

"""
    _export_code(gui) -> (code, n)

Returns the Julia code of `export_changes` and the number `n` of changed objects. The objects are
listed in the order of the component menu. Since moving a group moves its objects as well, the
changes of a group and its objects are replayed on an `_initial_copy` of the group, such that the
change of each object of the group is relative to its pose after the preceding lines.
"""
function _export_code(gui::LiveView)
    ctrl = gui.controls
    entries = _menu_entries(ctrl)
    names = _export_names(gui, first.(entries))
    lines = String[
        "# Changed poses of the live view, apply to the objects in their initial poses.",
        "# Each rotation is about the position of the object, groups are moved before their objects."]
    n = 0
    for top in ctrl.movable
        top isa LiveClipPlane && continue
        objs = _descendants(top)
        copies = top isa BMO.AbstractObjectGroup ? _descendants(_initial_copy(ctrl, top)) : nothing
        for (i, obj) in enumerate(objs)
            haskey(ctrl.init_poses, obj) || continue
            P, R = _pose(obj)
            P0, R0 = isnothing(copies) ? ctrl.init_poses[obj] : _pose(copies[i])
            (norm(P - P0) > 1e-15 || norm(R - R0) > 1e-15) || continue
            axis, angle = _rotation_axis_angle(R * R0')
            name = names[obj]
            type = string(nameof(typeof(obj)))
            label = get(gui.labels, obj, nothing)
            push!(lines, "", isnothing(label) ? "# $type" : "# $label ($type)")
            angle > 0 && push!(lines, "rotate3d!($name, $(_vector_code(axis)), $(repr(angle)))")
            push!(lines, "translate_to3d!($name, $(_vector_code(P)))")
            isnothing(copies) || _set_pose_exact!(copies[i], P, R)
            n += 1
        end
    end
    n == 0 && push!(lines, "", "# no changes")
    return join(lines, "\n") * "\n", n
end

"""Copies `code` to the clipboard, returns `false` if there is no clipboard, e.g. headless."""
function _copy_to_clipboard(code::String)
    try
        InteractiveUtils.clipboard(code)
        return true
    catch e
        e isa InterruptException && rethrow()
        return false
    end
end

"""
    export_changes(gui::LiveView; io = stdout, clipboard = false) -> String

Returns the changes of the poses in the `gui` as Julia code, which is printed to `io` and copied to
the clipboard if `clipboard` is `true` and a clipboard is available. The "Export" button of the
`gui` prints the code and copies it to the clipboard.

For each object whose pose differs from its initial pose, i.e. its pose when the window was
opened, the code contains a `rotate3d!` about the position of the object (skipped if the object
was only moved) followed by a `translate_to3d!` to its absolute position [m], with the full
precision of `Float64`. Applied to the objects in their initial poses, e.g. in the script that
created the system, the code reproduces the current poses. The objects of a group are listed after
the group, their changes are relative to the pose after moving the group. Clip planes are not
exported.

The variables are named after the `labels` of [`live_view`](@ref) if they are valid variable names,
otherwise `obj1`, `obj2`, … by the position of the object in the component menu. A comment above
each change names the label and the type of the object.

```julia
gui = live_view(system, beam; labels = Dict(m1 => "m1", lens => "lens"))
# move the objects, then
code = export_changes(gui)
```
"""
function export_changes(gui::LiveView; io::IO = stdout, clipboard::Bool = false)
    code, _ = _export_code(gui)
    print(io, code)
    clipboard && _copy_to_clipboard(code)
    return code
end

"""Prints the changed poses to `stdout` and copies them to the clipboard, see `export_changes`."""
function _export!(gui::LiveView)
    code, n = _export_code(gui)
    print(stdout, code)
    copied = gui.export_clipboard && _copy_to_clipboard(code)
    gui.status.text[] = "exported $n change$(n == 1 ? "" : "s")" * (copied ? " (copied)" : "")
    return nothing
end

#=
Component menu, hide and show
=#

"""Returns the options of the component menu of the `entries`, see `_menu_entries`."""
function _menu_options(labels, entries)
    isempty(entries) && return [("no components", 0)]
    return [("  "^depth * get(labels, obj, string(nameof(typeof(obj)))), i)
            for (i, (obj, depth)) in enumerate(entries)]
end

"""
    _on_menu_select!(gui, i)

Selects the object of the option `i` of the component menu, like a click in the 3D view. The
option `0`, i.e. no selection, is set by `_on_select!` and ignored.
"""
function _on_menu_select!(gui::LiveView, i)
    1 <= i <= length(gui.menu_objects) || return nothing
    ctrl = gui.controls
    obj = gui.menu_objects[i]
    ctrl.selected[] === obj && return nothing
    if ctrl.spectator[]
        gui.status.text[] = "spectator mode, press v to select components"
        gui.menu.i_selected[] = 0
        return nothing
    end
    ctrl.selected[] = obj
    _update_selection_box!(ctrl)
    gui.status.text[] = _pose_string(gui, obj)
    return nothing
end

"""Shows the selected object of the controls in the component menu and in the pose inspector."""
function _on_select!(gui::LiveView)
    obj = gui.controls.selected[]
    i = isnothing(obj) ? nothing : findfirst(o -> o === obj, gui.menu_objects)
    i = something(i, 0)
    gui.menu.i_selected[] == i || (gui.menu.i_selected[] = i)
    _update_inspector!(gui)
    return nothing
end

"""Sets the `visible` attribute of all plots of the rendered objects (leaves) of `obj`."""
function _set_hidden!(gui::LiveView, obj, hide::Bool)
    for leaf in _leaves(obj)
        hide ? push!(gui.hidden, leaf) : delete!(gui.hidden, leaf)
        i = findfirst(oh -> oh.obj === leaf, gui.controls.h.handles)
        isnothing(i) && continue
        for plot in gui.controls.h.handles[i].plots
            plot.visible[] == !hide || (plot.visible[] = !hide)
        end
    end
    return nothing
end

"""
    _toggle_hidden!(gui)

Hides the selected object of the `gui`, i.e. makes its plots invisible and clears the selection,
or shows it again if it is hidden. A hidden object can not be selected in the 3D view, but it
stays in the systems.
"""
function _toggle_hidden!(gui::LiveView)
    ctrl = gui.controls
    obj = ctrl.selected[]
    if isnothing(obj)
        gui.status.text[] = "select a component to hide it"
        return nothing
    elseif obj isa LiveClipPlane
        gui.status.text[] = "clip planes can not be hidden, press c to switch clipping off"
        return nothing
    end
    hide = !all(leaf -> leaf in gui.hidden, _leaves(obj))
    _set_hidden!(gui, obj, hide)
    if hide
        ctrl.selected[] = nothing
        _update_selection_box!(ctrl)
        gui.status.text[] = "$(_label(gui, obj)) hidden, select it in the menu to show it again"
    else
        gui.status.text[] = "$(_label(gui, obj)) shown"
    end
    return nothing
end

"""Shows all hidden objects of the `gui`."""
function _show_all!(gui::LiveView)
    foreach(leaf -> _set_hidden!(gui, leaf, false), collect(gui.hidden))
    gui.status.text[] = "all components shown"
    return nothing
end

#=
Pose inspector
=#

"""Labels of the boxes of the pose inspector: position [mm] and rotations about the axes [mrad]."""
const _POSE_FIELDS = ("x [mm]", "y [mm]", "z [mm]", "rx [mrad]", "ry [mrad]", "rv [mrad]")

"""Gizmo axes (see `_axis_vectors`) of the rotation boxes of the pose inspector."""
const _POSE_AXES = (:x, :y, :v)

"""Label colors of the pose inspector: the colors of the gizmo axes for the rotations."""
const _POSE_COLORS = (:black, :black, :black, :red, :green, :blue)

"""Returns `true` if a textbox or the component menu of the `gui` takes keyboard input."""
function _typing(gui::LiveView)
    gui.step_box.focused[] && return true
    any(tb -> tb.focused[], gui.pose_boxes) && return true
    return gui.menu.is_open[]
end

"""Shows `s` in the textbox `tb` without triggering its listeners, `""` shows the placeholder."""
function _set_box!(tb::Textbox, s::String)
    tb.stored_string.val = isempty(s) ? nothing : s
    tb.displayed_string[] == s || (tb.displayed_string[] = s)
    return nothing
end

"""
    _update_inspector!(gui; force = false)

Shows the position [mm] of the selected object of the `gui` in the position boxes of the pose
inspector and clears the rotation boxes, or clears all boxes if nothing is selected. Focused boxes
are skipped, unless `force`.
"""
function _update_inspector!(gui::LiveView; force::Bool = false)
    obj = gui.controls.selected[]
    p = isnothing(obj) ? nothing : 1e3 .* Vector{Float64}(position(obj))
    for (k, tb) in enumerate(gui.pose_boxes)
        tb.focused[] && !force && continue
        _set_box!(tb, isnothing(p) || k > 3 ? "" : string(round(p[k], digits = 6)))
    end
    return nothing
end

"""Returns `true` if the constraints of `obj` allow a move by `Δ`, see `kinematic_controls!`."""
function _move_allowed(ctrl::KinematicController, obj, Δ)
    haskey(_constraints_of(ctrl, obj), :move) || return true
    allowed = _allowed_axes(ctrl, obj, :move)
    isempty(allowed) && return iszero(Δ)
    A = hcat(_axis_vectors(ctrl, obj, allowed)...)
    return norm(Δ - A * (pinv(A) * Δ)) <= 1e-12 + 1e-9 * norm(Δ)
end

"""
    _apply_pose_input!(gui, k, s)

Applies the input `s` of the box `k` of the pose inspector (see `_POSE_FIELDS`) to the selected
object of the `gui`: for `k ≤ 3`, moves it to the absolute position x, y or z [mm], otherwise
rotates it about the red, green or blue axis of the controls [mrad], like a key step in the rotate
mode. The change is recorded in the undo history and solved like a key step. An invalid input
only shows a message in the status line.
"""
function _apply_pose_input!(gui::LiveView, k::Int, s)
    ctrl = gui.controls
    obj = ctrl.selected[]
    if isnothing(obj)
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
        translate_to3d!(obj, P)
        # `translate_to3d!` moves by `P - position`, which may round
        r = P - Vector{Float64}(position(obj))
        iszero(r) || translate3d!(obj, r)
    else
        sym = _POSE_AXES[k - 3]
        if !(sym in _allowed_axes(ctrl, obj, :rotate))
            gui.status.text[] = "$(_label(gui, obj)) can not rotate about this axis, see the constraints"
            _update_inspector!(gui; force = true)
            return nothing
        end
        iszero(x) || rotate3d!(obj, only(_axis_vectors(ctrl, obj, (sym,))), x / 1e3)
    end
    P1, R1 = _pose(obj)
    ctrl.last_key_step = nothing
    _push_history!(ctrl, obj, P0, R0, P1, R1)
    _request_update!(ctrl)
    _update_inspector!(gui; force = true)
    return nothing
end

"""Connects the export button, the component menu, the hide buttons and the pose inspector."""
function _connect_tools!(gui::LiveView)
    listeners = gui.controls.listeners
    push!(listeners, on(_ -> _export!(gui), gui.export_button.clicks))
    push!(listeners, on(i -> _on_menu_select!(gui, i), gui.menu.i_selected))
    push!(listeners, on(_ -> _on_select!(gui), gui.controls.selected))
    push!(listeners, on(_ -> _toggle_hidden!(gui), gui.hide_button.clicks))
    push!(listeners, on(_ -> _show_all!(gui), gui.show_all_button.clicks))
    for (k, tb) in enumerate(gui.pose_boxes)
        push!(listeners, on(s -> _apply_pose_input!(gui, k, s), tb.stored_string))
    end
    return nothing
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

Additional context, e.g. a static optomechanical assembly, can be added via `render!(gui.ax, ...)`.
The status line shows the pose of the moved object and its change since the window was opened.
The keyboard step can be typed into the textbox below the 3D view, e.g. `250 nm` or `50 µrad`,
where the unit selects the move or rotate mode.

# Component menu, pose inspector and export

The row below the status line holds a menu of all movable objects (the objects of a group
indented after the group, without clip planes), which selects an object like a click in the 3D
view. "hide" makes the plots of the selected object invisible and clears the selection, "hide"
on a hidden object selected in the menu shows it again, "show all" shows all objects. Hidden
objects can not be selected in the 3D view, but are still traced.

The pose inspector shows the position `x`, `y`, `z` [mm] of the selected object, `Enter` in a box
moves the object to the typed absolute coordinate. `rx`, `ry` and `rv` [mrad] rotate the object
about the red, green and blue axis of the controls, like the keys in the rotate mode. Each input is
recorded in the undo history, invalid inputs are reported in the status line.

The "Export" button prints the changed poses as Julia code to `stdout` and copies it to the
clipboard, see [`export_changes`](@ref).

# Adaptive tracing

If solving the systems takes longer than `trace_budget`, the objects still follow the mouse and the
keys immediately, while the beams are dimmed. The systems are solved once the movement pauses for
`idle_delay`. Likewise, detector panels that take longer than `trace_budget` show a preview on a
coarse grid while objects are moved, which is refined once the movement pauses.

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

# Keyword args

- `size = (1400, 800)`: size of the figure
- `auto_trace = true`: solves the systems after each change, otherwise only on request, see
  "Manual tracing"
- `detectors = :auto`: all `Detector`s of all systems. Alternatively a vector of `pd`,
  `pd => mode` or `pd => (mode, kwargs)`, where `mode` is `:auto`, `:spot` or `:intensity` and
  `kwargs` are passed to `intensity`, e.g. `(; n = 200, x_min = -1e-3, x_max = 1e-3, ...)`. An empty
  vector disables the panels.
- `on_change = (gui, obj) -> nothing`: called after each solve with the moved object, or
  `nothing` after a slider change
- `sliders = []`: vector of `"label" => (range, callback)` or `"label" => (range, callback, startvalue)`.
  The `callback` is called with the new value, then the systems are solved again.
- `system_kwargs = (;)`: passed to `live_render!` of each system
- `beam_kwargs = Dict()`: `beam => kwargs` passed to `live_render!` of the beam, by default
  `(; render_every = 5)` for beam groups
- `movable_sources = true`: shows an orange marker at each source, i.e. the beam or beam group of
  each pair, with which the source can be selected and moved like the components
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
- all other kwargs are passed to [`kinematic_controls!`](@ref), e.g. `fine_step`, `plane_normal`
  or `rotation_axis`
"""
function live_view(
        pairs::Pair{<:BMO.AbstractSystem}...;
        size = (1400, 800),
        auto_trace::Bool = true,
        detectors = :auto,
        on_change = (gui, obj) -> nothing,
        sliders = [],
        system_kwargs = (;),
        beam_kwargs = Dict(),
        movable_sources = true,
        labels = Dict(),
        trace_budget = 0.03,
        idle_delay = 0.2,
        clip_planes = [],
        clip_beams::Bool = false,
        view_cube::Bool = true,
        orthographic::Bool = false,
        lighting::Symbol = :studio,
        edges::Union{Nothing, Bool} = nothing,
        kwargs...
    )
    isempty(pairs) && throw(ArgumentError("live_view requires at least one system => beam pair"))
    ps = Pair{BMO.AbstractSystem, Any}[p for p in pairs]
    # several beams may share a system, which is rendered once
    systems = unique(objectid, first.(ps))
    specs = _panel_specs(detectors, systems)
    slider_specs = [_slider_spec(s) for s in sliders]
    clip_specs = _clip_plane_specs(clip_planes)

    fig = Figure(; size)
    ax = LScene(fig[1, 1]; show_axis = false)
    studio_lighting!(ax; preset = lighting)
    # Its click listener runs before the controls, hence clicks on the cube never select objects
    cube = view_cube ? view_cube!(ax) : nothing
    ncols = isempty(specs) ? 1 : 2

    # Detector panels in a near-square grid next to the 3D view
    panels = Any[]
    if !isempty(specs)
        grid = GridLayout(fig[1, 2])
        nc = ceil(Int, sqrt(length(specs)))
        for (i, (pd, mode, kw)) in enumerate(specs)
            parent = grid[(i - 1) ÷ nc + 1, (i - 1) % nc + 1]
            push!(panels, DetectorPanel(parent, pd, get(labels, pd, "Detector $i"), mode, kw))
        end
        colsize!(fig.layout, 1, Relative(0.6))
    end
    slider_grid = if isempty(slider_specs)
        nothing
    else
        SliderGrid(fig[2, 1:ncols], first.(slider_specs)...)
    end
    # Status row: trace button, auto trace toggle and status line
    status_row = GridLayout(fig[isnothing(slider_grid) ? 2 : 3, 1:ncols])
    trace_button = Button(status_row[1, 1]; label = "Trace (t)")
    auto_trace_toggle = Toggle(status_row[1, 2]; active = auto_trace)
    Label(status_row[1, 3], "auto trace")
    clip_beams_toggle = Toggle(status_row[1, 4]; active = clip_beams)
    Label(status_row[1, 5], "clip beams")
    orthographic_toggle = Toggle(status_row[1, 6]; active = orthographic)
    Label(status_row[1, 7], "orthographic")
    step_box = Textbox(status_row[1, 8]; placeholder = "step, e.g. 250 nm", width = 150)
    export_button = Button(status_row[1, 9]; label = "Export")
    status = Label(status_row[1, 10],
        "Click on a component to select it, press h to show the controls"; tellwidth = false)
    # Tool row: component menu (added below, once the movable objects are known), hide buttons and
    # pose inspector
    tool_row = GridLayout(fig[isnothing(slider_grid) ? 3 : 4, 1:ncols])
    hide_button = Button(tool_row[1, 2]; label = "hide")
    show_all_button = Button(tool_row[1, 3]; label = "show all")
    pose_boxes = Textbox[]
    for (k, (field, color)) in enumerate(zip(_POSE_FIELDS, _POSE_COLORS))
        Label(tool_row[1, 2k + 2], field; color)
        push!(pose_boxes, Textbox(tool_row[1, 2k + 3]; placeholder = k <= 3 ? " " : "0", width = 72))
    end
    # Keeps the tool row left-aligned
    # Keeps the tool row left-aligned and compact enough for narrow windows
    Label(tool_row[1, 16], ""; tellwidth = false)
    colgap!(tool_row, 6)

    # `edges` is only passed if given, i.e. custom `render!` methods of user objects do not need to
    # accept it
    sys_kw = isnothing(edges) ? system_kwargs : (; edges, system_kwargs...)
    system_handles = SystemRenderHandle[live_render!(ax, sys; sys_kw...) for sys in systems]
    beam_handles = AbstractRenderHandle[]
    for beam in last.(ps)
        default = beam isa BMO.AbstractBeamGroup ? (; render_every = 5) : (;)
        kw = get(beam_kwargs, beam, default)
        # The planes of the beams are set explicitly by `_apply_clip_planes!`, see `clip_beams`
        push!(beam_handles, live_render!(ax, beam; kw..., clip_planes = Plane3f[]))
    end

    # A single controller for all systems, otherwise several controllers would compete for events
    handles = reduce(vcat, [h.handles for h in system_handles]; init = ObjectRenderHandle[])
    # Size of the systems, before any clip plane shrinks the bounding boxes
    plots = reduce(vcat, (oh.plots for oh in handles); init = AbstractPlot[])
    extent = isempty(plots) ? 0.125 :
             maximum(GeometryBasics.widths(mapreduce(Makie.boundingbox, GeometryBasics.union, plots)))
    if movable_sources
        # Markers of the sources, scaled to the size of the systems
        marker_size = 0.08 * extent
        for src in unique(objectid, last.(ps))
            BMO._is_static(src) || push!(handles, _live_render_source!(ax, src; size = marker_size))
        end
    end
    parent = IdDict{BMO.AbstractObject, BMO.AbstractObject}()
    foreach(h -> merge!(parent, h.parent), system_handles)
    combined = SystemRenderHandle(ax, first(systems), handles, parent)
    gui_ref = Ref{LiveView}()
    # Moving a clip plane only re-applies the planes, the systems are not solved
    change = function (obj)
        gui = gui_ref[]
        obj isa LiveClipPlane ? _on_clip_change!(gui, obj) : _on_change!(gui, obj)
        _update_inspector!(gui)
        return nothing
    end
    # Typing into a textbox or the search of the component menu must not trigger the controls
    controls = kinematic_controls!(ax, combined; on_change = change,
        ignore_keys = () -> isassigned(gui_ref) && _typing(gui_ref[]), kwargs...)

    labels_dict = IdDict{Any, String}(labels)
    entries = _menu_entries(controls)
    menu = Menu(tool_row[1, 1]; options = _menu_options(labels_dict, entries), default = nothing,
        prompt = "select component", width = 170)
    gui = LiveView(fig, ax, ps, system_handles, beam_handles, controls, panels, status, slider_grid,
        on_change, nothing, auto_trace_toggle.active, false, trace_button, auto_trace_toggle,
        IdDict{Any, Any}(), Float64(trace_budget), Float64(idle_delay), 0.0, 0.0, false, nothing,
        0.0, false, labels_dict, step_box, LiveClipPlane[], true, 1.2 * extent,
        clip_beams, clip_beams_toggle, cube, orthographic_toggle, export_button, true, menu,
        Any[first.(entries)...], hide_button, show_all_button, Base.IdSet{Any}(), pose_boxes)
    gui_ref[] = gui
    for (point, normal) in clip_specs
        _add_clip_plane!(gui, point, normal; select = false)
    end
    isnothing(slider_grid) || _connect_sliders!(gui, last.(slider_specs))
    _connect_trace!(gui)
    _connect_clip_planes!(gui)
    _connect_tools!(gui)
    push!(controls.listeners, on(v -> v == gui.clip_beams || _set_clip_beams!(gui, v), clip_beams_toggle.active))
    push!(controls.listeners, on(s -> _set_step!(gui, s), step_box.stored_string))
    push!(controls.listeners, on(v -> _set_orthographic!(gui, v), orthographic_toggle.active))
    _set_orthographic!(gui, orthographic)
    _resolve!(gui, nothing)
    # Initial view from the Front-Right-Top corner, in which the labels of the view cube read
    # correctly. Only set once, later changes of the view, e.g. via `set_view`, are kept.
    cam = cameracontrols(ax.scene)
    lookat = Vector{Float64}(cam.lookat[])
    dist = norm(Vector{Float64}(cam.eyeposition[]) .- lookat)
    o, up = _region_view((1, -1, 1))
    set_view(ax, lookat .+ dist .* o, lookat, up)
    return gui
end

live_view(system::BMO.AbstractSystem, beam; kwargs...) = live_view(system => beam; kwargs...)

"""Switches the 3D view of the `gui` to orthographic (`true`) or perspective (`false`) projection."""
function _set_orthographic!(gui::LiveView, orthographic::Bool)
    settings = cameracontrols(gui.ax.scene).settings
    settings.projectiontype[] = orthographic ? Makie.Orthographic : Makie.Perspective
    return nothing
end
