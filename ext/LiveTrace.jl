#=
Tracing of the live view: preview and full solves, background jobs, stale beams, the keyboard
step and the sliders
=#

_render_every(h::BeamRenderHandle) = h.render_every
_render_every(h::AstigmaticGroupRenderHandle) = h.render_every
_render_every(_) = 1

"""Returns `true` if the `beam` with the render handle `h` is solved as a preview while moving."""
_previewable(beam, h) = beam isa BMO.AbstractBeamGroup && _render_every(h) > 1

"""Returns `true` if the `gui` solves any of its beams as a preview while moving, see `_resolve!`."""
_has_preview(gui::LiveView) = gui.trace.preview_enabled &&
                              any(i -> _previewable(gui.pairs[i].second, gui.beam_handles[i]), eachindex(gui.pairs))

"""
    _solve_from_start!(system, beam)

Solves the `beam` (a beam or a beam group) in the `system` by brute force: from its start ray, with
`retrace = false`, since retracing the path of the last solve is not reliable in the live view.
Each beam is reset to its untraced start state first (`empty!`): without retracing,
`solve_system!` does not trace a beam whose last ray already ends on a surface.
"""
function _solve_from_start!(system, beam::BMO.AbstractBeam)
    empty!(beam)
    solve_system!(system, beam; retrace = false)
    return nothing
end
function _solve_from_start!(system, bg::BMO.AbstractBeamGroup)
    foreach(empty!, BMO.beams(bg))
    solve_system!(system, bg; retrace = false)
    return nothing
end

"""
    _solve_preview!(system, bg, k)

Solves only the rendered beams `beams(bg)[1:k:end]` of the beam group `bg` and resets all other
beams to their untraced start state, such that no outdated paths are rendered or hit a detector.
"""
function _solve_preview!(system, bg::BMO.AbstractBeamGroup, k::Int)
    bms = BMO.beams(bg)
    idx = 1:k:length(bms)
    # Like `solve_system!` of a beam group
    Threads.@threads for i in idx
        _solve_from_start!(system, bms[i])
    end
    for i in eachindex(bms)
        (i - 1) % k == 0 || empty!(bms[i])
    end
    return nothing
end

"""
    _compute(pairs, handles, panels, sinks; coarse = false, preview = false)

Empties all `Detector`s of the systems of the `pairs`, solves the systems and computes the fields
of the detector `panels` (see `_panel_field`), without changing any plot, such that it can run in
a background task, see `_solve!`. `handles` are the render handles of the beams of the `pairs`.
Each source is traced, and each field computed, with its progress output `sinks[k]` (see
`BMO._PROGRESS_SINK`, `nothing` for the terminal): first one per pair, then one per panel.

With `preview`, beam groups rendered with `render_every > 1` are solved only for their rendered
beams, see `_solve_preview!`. Returns `(; previewed, panels, fields, solve_time, field_time)`:
whether a beam group was solved as a preview, the `panels` and their fields and the durations [s].
"""
function _compute(pairs, handles, panels, sinks; coarse = false, preview = false)
    # Monotonic clock with ns resolution, time() is too coarse on Windows for fast solves
    t0 = time_ns()
    # A detector can be part of several systems, hence empty all before solving
    foreach(empty!, _find_detectors(first.(pairs)))
    previewed = preview && any(i -> _previewable(pairs[i].second, handles[i]), eachindex(pairs))
    for (i, (sys, beam)) in enumerate(pairs)
        h = handles[i]
        Base.ScopedValues.with(BMO._PROGRESS_SINK => sinks[i]) do
            if preview && _previewable(beam, h)
                _solve_preview!(sys, beam, _render_every(h))
            else
                _solve_from_start!(sys, beam)
            end
        end
    end
    t1 = time_ns()
    n = length(pairs)
    fields = Any[Base.ScopedValues.with(() -> _panel_field(p, coarse), BMO._PROGRESS_SINK => sinks[n + k])
                 for (k, p) in enumerate(panels)]
    return (; previewed, panels, fields, solve_time = 1e-9 * (t1 - t0),
        field_time = 1e-9 * (time_ns() - t1))
end

"""
    _apply!(gui, r, obj; coarse = false)

Shows the result `r` of `_compute` in the `gui`: updates the beams, the detector panels, the
durations of the adaptive tracing and the status line, and calls the user `on_change` with the
moved `obj` (or `nothing`) after a full solve.

After a preview (see `_compute`), the titles of the detector panels are marked with "(preview)"
and `gui.trace.preview` is set, such that the full solve follows once the movement pauses, see
`_on_idle!`. `on_change` is only called after full solves, the metrics of a preview are not
recorded in the history of the panels.
"""
function _apply!(gui::LiveView, r, obj; coarse = false)
    t0 = time_ns()
    foreach(update_render!, gui.beam_handles)
    t1 = time_ns()
    previewed = r.previewed
    for (p, field) in zip(r.panels, r.fields)
        _apply_panel!(gui, p, field; coarse, preview = previewed)
    end
    solve_time = r.solve_time + 1e-9 * (t1 - t0)
    if previewed
        gui.trace.preview_time = solve_time
    else
        gui.trace.solve_time = solve_time
        coarse || (gui.trace.panel_time = r.field_time + 1e-9 * (time_ns() - t1))
    end
    gui.trace.coarse = coarse
    gui.trace.preview = previewed
    gui.trace.preview_obj = obj
    if !previewed
        try
            gui.on_change(gui, obj)
            gui.last_error = nothing
        catch e
            gui.last_error = _log_once(e, gui.last_error, "`on_change` callback")
        end
        # After `on_change`, which may record what the panels of `add_panel!` show
        _update_user_panels!(gui)
    end
    isnothing(obj) || (gui.status.text[] = _pose_string(gui, obj))
    _on_solved!(gui)
    _on_applied!(gui)
    # e.g. values of the last solve on the cards
    _update_inspector!(gui)
    return nothing
end

"""
    _resolve!(gui::LiveView, obj; coarse = false, preview = false)

Empties all `Detector`s, solves all systems and updates the beams, detector panels and the status
line of the `gui`, see `_compute` and `_apply!`. The user `on_change` is called with the moved
`obj`, or `nothing`. Unlike `_solve!`, it returns only after the solve, which runs on the calling
task; a solve of the `gui` in the background is cancelled first.
"""
function _resolve!(gui::LiveView, obj; coarse = false, preview = false)
    _cancel_solve!(gui)
    panels = _computed_panels(gui, preview)
    sinks = fill(nothing, length(gui.pairs) + length(panels))
    _on_solve_started!(gui)
    r = _compute(gui.pairs, gui.beam_handles, panels, sinks; coarse, preview)
    _apply!(gui, r, obj; coarse)
    return nothing
end

"""
    _start_job(gui, apply, obj, pairs, handles[, panels]; coarse = false, preview = false, timing) -> _SolveJob

Starts `_compute` for the `pairs` (with the beam render `handles`) and the detector `panels` (by
default `_computed_panels`) of the `gui` in a background task, with a progress sink per source and
panel, see `_SolveJob`. `apply` shows the result, `timing` is the duration field that a cancelled
job updates, `:panel_time` for a job that only computes panels, i.e. without `pairs`.
"""
_start_job(gui::LiveView, apply, obj, pairs, handles; preview = false, kwargs...) =
    _start_job(gui, apply, obj, pairs, handles, _computed_panels(gui, preview); preview, kwargs...)

function _start_job(gui::LiveView, apply, obj, pairs, handles, panels; coarse = false,
        preview = false, timing::Symbol)
    isempty(pairs) || _on_solve_started!(gui)
    # The task works on its own copies of the lists, the objects are protected by `_change!`
    pairs, handles, panels = copy(pairs), copy(handles), copy(panels)
    sinks = [BMO._ProgressSink() for _ in 1:(length(pairs) + length(panels))]
    anchors = Point3f[_progress_anchor.(last.(pairs)); _progress_anchor.(getfield.(panels, :pd))]
    done = Base.Event()
    task = Threads.@spawn try
        _compute(pairs, handles, panels, sinks; coarse, preview)
    finally
        notify(done)
    end
    return _SolveJob(task, done, sinks, anchors, apply, obj, timing, time(),
        (; k = 0, t0 = NaN, t = NaN, count = 0))
end

"""
    _progress_anchor(x)

Position of the progress window of the source or detector `x`, see `_SolveJob`. `NaN` for other
types, e.g. custom sources without a position, whose window is shown at the bottom edge of the
view, see `_screen_anchor`.
"""
_progress_anchor(x::Union{BMO.AbstractBeam, BMO.AbstractBeamGroup, BMO.AbstractObject}) =
    Point3f(position(x))
_progress_anchor(_) = Point3f(NaN)

"""
    _run!(gui, job, msg)

Waits for the `job` at most `gui.trace.progress_delay`: if it is done by then, its result is shown at
once, see `_finish!`, and `true` is returned on success. Otherwise the job continues in the
background as `gui.trace.job` (shown once it is done, see `_poll_job!`), the status line shows `msg`
and `false` is returned. Solves up to `progress_delay` thus behave as if they ran on the render
task, only longer ones keep the window responsive and show their progress.
"""
function _run!(gui::LiveView, job::_SolveJob, msg::AbstractString)
    _wait(job.done, gui.trace.progress_delay)
    istaskdone(job.task) && return _finish!(gui, job)
    gui.trace.job = job
    gui.status.text[] = msg
    return false
end

"""Waits for the event `done`, at most `timeout` [s]."""
function _wait(done::Base.Event, timeout::Real)
    isfinite(timeout) || return wait(done)
    timer = Timer(_ -> notify(done), max(timeout, 0.0))
    wait(done)
    close(timer)
    return nothing
end

_running(gui::LiveView) = _running(gui.trace.job)
_running(::Nothing) = false
_running(::_SolveJob) = true

"""
    _cancel_solve!(gui::LiveView)

Cancels the solve of the `gui` that runs in the background, if any: its loops stop after their
current item, see `BMO._ProgressSink`. Waits for the task, discards its result, marks the beams and
detector panels as outdated and counts the elapsed time as the duration of the solve, such that
further changes defer the solve until the movement pauses, see `_on_change!`. A deferred solve,
preview or coarse panel is not completed afterwards, the next change or `t` solves again.
"""
_cancel_solve!(gui::LiveView) = _cancel!(gui, gui.trace.job)
_cancel!(::LiveView, ::Nothing) = nothing

function _cancel!(gui::LiveView, job::_SolveJob)
    foreach(s -> s.cancel[] = true, job.sinks)
    try
        wait(job.task)
    catch
        # cancelled or failed, the result is discarded either way
    end
    gui.trace.job = nothing
    _hide_progress!(gui.trace.progress)
    setproperty!(gui.trace, job.timing, max(getproperty(gui.trace, job.timing), time() - job.t0))
    gui.trace.pending = gui.trace.preview = gui.trace.coarse = false
    _mark_stale!(gui, nothing; msg = _CANCELLED)
    return nothing
end

const _CANCELLED = "trace cancelled, press t to trace"

"""
    _finish!(gui, job)

Shows the result of the finished `job` in the `gui` via `job.apply`, and restores the appearance
of the beams after a solve (a job that only computed detector panels leaves them as they are, see
`_solves`). If solving failed, the error is logged and the beams and detector panels are kept
marked as outdated. Returns `true` on success.
"""
function _finish!(gui::LiveView, job::_SolveJob)
    gui.trace.job === job && (gui.trace.job = nothing)
    _hide_progress!(gui.trace.progress)
    try
        job.apply(fetch(job.task))
    catch e
        _fail!(gui, e)
        return false
    end
    if _solves(job)
        _restore_beams!(gui)
        gui.trace.stale = false
    end
    return true
end

"""Whether the `job` solves the systems, i.e. does not only compute detector panels."""
_solves(job::_SolveJob) = job.timing !== :panel_time

"""Marks the beams and detector panels of the `gui` as outdated after the solve failed with `e`."""
function _fail!(gui::LiveView, e)
    if BMO._is_cancelled(e)
        _mark_stale!(gui, nothing; msg = _CANCELLED)
        return nothing
    end
    gui.last_error = _log_once(_task_error(e), gui.last_error, "solving the systems")
    gui.trace.stale || _dim_beams!(gui)
    gui.trace.stale = true
    gui.status.text[] = "solving the systems failed, see the log"
    return nothing
end

_task_error(e::TaskFailedException) = _task_error(e.task.result)
_task_error(e) = e

"""
    _poll_job!(gui::LiveView)

Called every frame: shows the result of the solve of the `gui` in the background once it is done,
see `_finish!`, and until then the progress window of its running loop, see `_show_loop!`.
"""
_poll_job!(gui::LiveView) = _poll!(gui, gui.trace.job)
_poll!(::LiveView, ::Nothing) = nothing

function _poll!(gui::LiveView, job::_SolveJob)
    if istaskdone(job.task)
        _finish!(gui, job) && isnothing(job.obj) && (gui.status.text[] = "traced")
        return nothing
    end
    shown = any(k -> _show_loop!(gui, job, k, BMO._progress_state(job.sinks[k])),
        eachindex(job.sinks))
    shown || _hide_progress!(gui.trace.progress)
    return nothing
end

"""
    _show_loop!(gui, job, k, state)

Shows the progress window of the loop of the sink `k` of the `job` with the `state` of
`BMO._progress_state`, at the position of its source or detector, once the loop has run for
`gui.trace.progress_delay`, like the terminal bars after `get_progress_threshold()`. Returns `true` if
the window is shown.
"""
_show_loop!(::LiveView, ::_SolveJob, ::Int, ::Nothing) = false

function _show_loop!(gui::LiveView, job::_SolveJob, k::Int, state::NamedTuple)
    t = time()
    t - state.t0 >= gui.trace.progress_delay || return false
    s = job.shown
    (s.k == k && s.t0 == state.t0) || (job.shown = s = (; k, state.t0, t, state.count))
    _show_progress!(gui.trace.progress, job.anchors[k], state.count / max(state.n, 1),
        _progress_label(state, s, t))
    return true
end

"""
    _progress_label(state, shown, t)

Label of the progress window of a loop with the `state` of `BMO._progress_state`, e.g.
"Tracing beams 42 % · 3 s" with the remaining time. It uses the rate since the window appeared
(`shown`, see `_SolveJob`), like the terminal bars, and is left out until that rate is known.
"""
function _progress_label(state, shown, t)
    pct = floor(Int, 100 * state.count / max(state.n, 1))
    label = "$(state.desc) $pct %"
    done = state.count - shown.count
    (done > 0 && t > shown.t) || return label
    left = (state.n - state.count) * (t - shown.t) / done
    return "$label · $(_duration_string(left))"
end

"""Formats the duration `s` [s] as e.g. `3 s` or `2:05` (minutes and seconds)."""
function _duration_string(s::Real)
    s < 59.5 && return "$(round(Int, s)) s"
    m, r = divrem(round(Int, s), 60)
    return "$m:$(lpad(r, 2, '0'))"
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
        haskey(gui.trace.beam_alphas, plot) || (gui.trace.beam_alphas[plot] = plot.alpha[])
        plot.alpha[] = _STALE_ALPHA
    end
    return nothing
end

"""Restores the `alpha` of all beam plots of the `gui` after dimming."""
function _restore_beams!(gui::LiveView)
    for (plot, alpha) in gui.trace.beam_alphas
        plot.alpha[] = alpha
    end
    empty!(gui.trace.beam_alphas)
    return nothing
end

"""Marks the beams and detector panels of the `gui` as outdated after `obj` (or a slider) changed."""
function _mark_stale!(gui::LiveView, obj; msg = "outdated, press t to trace")
    gui.trace.stale || _dim_beams!(gui)
    gui.trace.stale = true
    gui.status.text[] = isnothing(obj) ? msg : "$(_pose_string(gui, obj)) — $msg"
    return nothing
end

"""
    _solve!(gui::LiveView, obj; coarse = false, preview = false)

Solves all systems of the `gui` like `_resolve!`, but in a background task: a solve that takes
longer than `progress_delay` continues in the background, with progress windows at the sources and
detectors, while the window stays responsive, see `_run!`. A running solve is cancelled first.
Afterwards the appearance of the beams is restored; if solving fails, the beams and detector
panels are kept marked as outdated. Returns `true` if the solve succeeded without continuing in
the background.
"""
function _solve!(gui::LiveView, obj; coarse = false, preview = false)
    _cancel_solve!(gui)
    gui.trace.pending = false
    job = _start_job(gui, r -> _apply!(gui, r, obj; coarse), obj, gui.pairs, gui.beam_handles;
        coarse, preview, timing = preview ? :preview_time : :solve_time)
    done = _run!(gui, job, "tracing, Esc cancels")
    if _running(gui)
        # Outdated until the solve in the background is shown, see `_finish!`
        gui.trace.stale || _dim_beams!(gui)
        gui.trace.stale = true
    end
    return done
end

"""Shows the fields `r.fields` of the detector panels `r.panels`, which refine a coarse preview."""
function _refine!(gui::LiveView, r)
    for (p, field) in zip(r.panels, r.fields)
        _update_panel!(p, field; record = false)
    end
    gui.trace.coarse = false
    return nothing
end

"""
    _on_change!(gui::LiveView, obj)

Called after each change of `obj` (or a slider, then `obj` is `nothing`). Solves the systems, with
a preview of beam groups (see `_resolve!`) and a coarse preview of slow detector panels, or marks
them as outdated. If solving (the preview solve, if any) is slower than the `trace_budget`, the
solve is deferred until the movement pauses, see `_on_idle!`.
"""
function _on_change!(gui::LiveView, obj)
    gui.trace.last_change = time()
    preview = _has_preview(gui)
    if !gui.trace.auto[]
        _mark_stale!(gui, obj)
    elseif (preview ? gui.trace.preview_time : gui.trace.solve_time) <= gui.trace.budget
        _solve!(gui, obj; coarse = gui.trace.panel_time > gui.trace.budget, preview)
    else
        _mark_stale!(gui, obj; msg = "tracing when the movement pauses")
        gui.trace.pending = true
        gui.trace.pending_obj = obj
    end
    return nothing
end

"""
Solves deferred changes, completes a preview solve with a full solve and refines the preview of
the detector panels once the movement pauses and no solve runs in the background.
"""
function _on_idle!(gui::LiveView)
    time() - gui.trace.last_change > gui.trace.idle_delay || return nothing
    _running(gui) && return nothing
    if gui.trace.pending && gui.trace.auto[]
        _solve!(gui, gui.trace.pending_obj)
    elseif gui.trace.preview
        # Also if auto tracing was switched off in the meantime, since the preview is incomplete
        _solve!(gui, gui.trace.preview_obj)
    elseif gui.trace.coarse
        job = _start_job(gui, r -> _refine!(gui, r), gui.trace.preview_obj, empty(gui.pairs),
            empty(gui.beam_handles), _shown_panels(gui); timing = :panel_time)
        _run!(gui, job, "computing the detector fields, Esc cancels")
    end
    return nothing
end

"""
Solves all systems of the `gui` on request, with the currently selected object. Ignored while a
solve runs in the background.
"""
function _trace!(gui::LiveView)
    _running(gui) && return nothing
    obj = gui.controls.selected[]
    # A clip plane is not part of the systems
    obj isa LiveClipPlane && (obj = nothing)
    _solve!(gui, obj) && isnothing(obj) && (gui.status.text[] = "traced")
    return nothing
end

"""
Connects the trace button, the key `t`, the auto trace toggle and the key `Esc`, which cancels a
solve in the background, of the `gui`.
"""
function _connect_trace!(gui::LiveView)
    listeners = gui.controls.listeners
    scene = gui.ax.scene
    push!(listeners, on(_ -> _trace!(gui), gui.widgets.trace_button.clicks))
    push!(listeners, on(events(scene).keyboardbutton, priority = 200) do event
        (event.action == Keyboard.press && event.key == Keyboard.t) || return Consume(false)
        gui.controls.ignore_keys() && return Consume(false)
        _trace!(gui)
        return Consume(true)
    end)
    # Before the beam inspection and the controls, which use Esc as well
    push!(listeners, on(events(scene).keyboardbutton, priority = 202) do event
        (event.action == Keyboard.press && event.key == Keyboard.escape) || return Consume(false)
        (_running(gui) && !gui.controls.ignore_keys()) || return Consume(false)
        _cancel_solve!(gui)
        return Consume(true)
    end)
    push!(listeners, on(gui.trace.auto) do active
        active && gui.trace.stale && _trace!(gui)
        return nothing
    end)
    push!(listeners, on(events(scene).tick) do _
        _poll_job!(gui)
        _on_idle!(gui)
        return nothing
    end)
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
                # The callbacks may change any object
                _change!(() -> callbacks[i](v), gui.controls, nothing)
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

