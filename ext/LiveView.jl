using Makie: Figure, Axis, Label, SliderGrid, GridLayout, DataAspect, Relative, colsize!,
             heatmap!, autolimits!, Button, Toggle, Textbox

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
of the 3D view in `ax` and the `KinematicController` in `controls`. Use `display` to show the
window and `close` to remove the controls.

The systems are solved after each change if `auto_trace[]` is `true`, otherwise only via the
`trace_button`, the key `t` or by switching the `auto_trace_toggle` on. `stale` is `true` if the
beams and detector panels do not match the current poses of the objects. If solving takes longer
than `trace_budget` [s], the systems are solved once the movement pauses for `idle_delay` [s].
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
end

Base.display(gui::LiveView) = display(gui.fig)

function Base.show(io::IO, gui::LiveView)
    print(io, "LiveView(", length(gui.pairs), " systems, ", length(gui.panels), " detector panels)")
end

Base.close(gui::LiveView) = close(gui.controls)

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
    _solve!(gui, obj) && isnothing(obj) && (gui.status.text[] = "traced")
    return nothing
end

"""Connects the trace button, the key `t` and the auto trace toggle of the `gui`."""
function _connect_trace!(gui::LiveView)
    listeners = gui.controls.listeners
    push!(listeners, on(_ -> _trace!(gui), gui.trace_button.clicks))
    push!(listeners, on(events(gui.ax.scene).keyboardbutton, priority = 200) do event
        (event.action == Keyboard.press && event.key == Keyboard.t) || return Consume(false)
        gui.step_box.focused[] && return Consume(false)
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
- `labels = Dict()`: `obj => "name"` for the status line and the titles of the detector panels
- `trace_budget = 0.03`: [s] duration of a solve or panel update, above which tracing is deferred
  or the panels show a preview, see "Adaptive tracing"
- `idle_delay = 0.2`: [s] pause of the movement after which deferred tracing runs
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
        kwargs...
    )
    isempty(pairs) && throw(ArgumentError("live_view requires at least one system => beam pair"))
    ps = Pair{BMO.AbstractSystem, Any}[p for p in pairs]
    systems = first.(ps)
    specs = _panel_specs(detectors, systems)
    slider_specs = [_slider_spec(s) for s in sliders]

    fig = Figure(; size)
    ax = LScene(fig[1, 1]; show_axis = false)
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
    step_box = Textbox(status_row[1, 4]; placeholder = "step, e.g. 250 nm", width = 150)
    status = Label(status_row[1, 5],
        "Click on a component to select it, press h to show the controls"; tellwidth = false)

    system_handles = SystemRenderHandle[live_render!(ax, sys; system_kwargs...) for sys in systems]
    beam_handles = AbstractRenderHandle[]
    for beam in last.(ps)
        default = beam isa BMO.AbstractBeamGroup ? (; render_every = 5) : (;)
        push!(beam_handles, live_render!(ax, beam; get(beam_kwargs, beam, default)...))
    end

    # A single controller for all systems, otherwise several controllers would compete for events
    handles = reduce(vcat, [h.handles for h in system_handles]; init = ObjectRenderHandle[])
    if movable_sources
        # Markers of the sources, scaled to the size of the systems
        plots = reduce(vcat, (oh.plots for oh in handles); init = AbstractPlot[])
        marker_size = isempty(plots) ? 1e-2 :
                      0.08 * maximum(GeometryBasics.widths(mapreduce(Makie.boundingbox, GeometryBasics.union, plots)))
        for src in unique(objectid, last.(ps))
            BMO._is_static(src) || push!(handles, _live_render_source!(ax, src; size = marker_size))
        end
    end
    parent = IdDict{BMO.AbstractObject, BMO.AbstractObject}()
    foreach(h -> merge!(parent, h.parent), system_handles)
    combined = SystemRenderHandle(ax, first(systems), handles, parent)
    gui_ref = Ref{LiveView}()
    # Typing into the step textbox must not trigger the controls
    controls = kinematic_controls!(ax, combined; on_change = obj -> _on_change!(gui_ref[], obj),
        ignore_keys = () -> step_box.focused[], kwargs...)

    gui = LiveView(fig, ax, ps, system_handles, beam_handles, controls, panels, status, slider_grid,
        on_change, nothing, auto_trace_toggle.active, false, trace_button, auto_trace_toggle,
        IdDict{Any, Any}(), Float64(trace_budget), Float64(idle_delay), 0.0, 0.0, false, nothing,
        0.0, false, IdDict{Any, String}(labels), step_box)
    gui_ref[] = gui
    isnothing(slider_grid) || _connect_sliders!(gui, last.(slider_specs))
    _connect_trace!(gui)
    push!(controls.listeners, on(s -> _set_step!(gui, s), step_box.stored_string))
    _resolve!(gui, nothing)
    return gui
end

live_view(system::BMO.AbstractSystem, beam; kwargs...) = live_view(system => beam; kwargs...)
