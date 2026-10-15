using Base.ScopedValues: ScopedValue

"""
    _LazyProgress{O}

Thread-safe progress counter for long loops. Where the progress goes depends on the output
type `O`, see `_report!`:

- an `IO`: a `ProgressMeter.Progress` bar is drawn only once the loop has run for `threshold`
  seconds, see [`get_progress_threshold`](@ref),
- a `ProgressSink`: nothing is drawn, the sink exposes the running loop to another task (e.g.
  the live view) and can cancel it.

[`_tick!`](@ref) costs one atomic increment plus the report of the output: one `time()` call
for an `IO`, one atomic load of the cancel flag for a sink. The terminal bar is created lazily
and redrawn at most every `dt` seconds by one thread at a time (`trylock`), so an invisible bar
never contends for a lock. Use it through [`_with_progress`](@ref), which finishes the bar, or
cancels it if the loop throws.

# Fields

- `n`: total number of items
- `desc`: description printed in front of the bar
- `enabled`: if `false`, [`_tick!`](@ref) returns after a single branch
- `output`: stream the bar is drawn to, or the `ProgressSink` the loop reports to
- `count`: number of finished items
- `tnext`: earliest `time()` at which the bar is created or redrawn, `Inf` stops drawing
- `dt`: minimum interval in s between redraws
- `t0`: `time()` at which the counter was created, i.e. the loop started
- `lock`: serializes creating, redrawing and stopping the bar
- `bar`: the ProgressMeter bar, `nothing` until it is first drawn
"""
mutable struct _LazyProgress{O}
    const n::Int
    const desc::String
    const enabled::Bool
    const output::O
    const count::Threads.Atomic{Int}
    const tnext::Threads.Atomic{Float64}
    const dt::Float64
    const t0::Float64
    const lock::Threads.SpinLock
    bar::Nullable{Progress}
end

"""
    _LazyProgress(n, desc; enabled = true, output = stderr, threshold = get_progress_threshold(), dt = 0.2)

Progress counter for `n` items that reports to `output`, an `IO` or a `ProgressSink`. An `IO`
is drawn to once `threshold` seconds have passed; a sink ignores `threshold` and `dt`.
"""
function _LazyProgress(n::Integer, desc::AbstractString; enabled::Bool = true,
        output = stderr, threshold::Real = get_progress_threshold(), dt::Real = 0.2)
    t0 = time()
    return _LazyProgress(n, desc, enabled, output, Threads.Atomic{Int}(0),
        Threads.Atomic{Float64}(t0 + threshold), Float64(dt), t0, Threads.SpinLock(), nothing)
end

"""
    _ProgressCancelled()

Thrown by [`_tick!`](@ref) inside a loop whose `ProgressSink` was cancelled. Check caught
exceptions with `is_cancelled`, which also sees through the wrappers of `Threads.@threads`.
"""
struct _ProgressCancelled <: Exception end

"""
    ProgressSink()

Receiver of the progress of all progress loops of BeamletOptics (e.g. of `solve_system!` with
`progress = true`, see [`_with_progress`](@ref)) that run inside
`with(PROGRESS_SINK => sink) do … end`, including the tasks of `Threads.@threads` started
there. Such loops draw no terminal bar. Part of the developer API, e.g. for a GUI that shows the
progress of a solve in the background and cancels it. Another task reads the running loop with
`progress_state(sink)`, and cancels the loops by setting `sink.cancel[] = true`: the next
[`_tick!`](@ref) then throws `_ProgressCancelled`, i.e. a loop stops after its current item.

# Fields

- `current`: the innermost running loop, `nothing` between loops (atomic)
- `cancel`: cancel request, checked by every tick
"""
mutable struct ProgressSink
    @atomic current::Union{Nothing, _LazyProgress}
    const cancel::Threads.Atomic{Bool}
    ProgressSink() = new(nothing, Threads.Atomic{Bool}(false))
end

"""
    PROGRESS_SINK

Scoped `ProgressSink` of the running code, `nothing` (terminal bars) by default. Set it with
`Base.ScopedValues.with(PROGRESS_SINK => sink) do … end`.
"""
const PROGRESS_SINK = ScopedValue{Union{Nothing, ProgressSink}}(nothing)

"""
    _LazyProgress(progress::Bool, n, desc)

Progress counter for the `progress` keyword of a public function. It reports to the scoped
`PROGRESS_SINK` if one is set, otherwise to `stderr`. It is enabled only if `progress` is
`true` and the output can show it: a sink always, `stderr` only if it is a terminal
(`Base.TTY`), so documentation builds, CI logs and piped output stay clean.
"""
function _LazyProgress(progress::Bool, n::Integer, desc::AbstractString)
    output = _progress_output(PROGRESS_SINK[])
    return _LazyProgress(n, desc; enabled = progress && _can_draw(output), output)
end

_progress_output(::Nothing) = stderr
_progress_output(s::ProgressSink) = s

_can_draw(::Base.TTY) = true
_can_draw(::IO) = false
_can_draw(::ProgressSink) = true

"""
    _tick!(p::_LazyProgress)

Count one finished item and report it to the output of `p`, see `_report!`. Safe to call from
any thread. Throws `_ProgressCancelled` if `p` reports to a cancelled `ProgressSink`.
"""
@inline function _tick!(p::_LazyProgress)
    p.enabled || return nothing
    Threads.atomic_add!(p.count, 1)
    return _report!(p.output, p)
end

"""
    _report!(output::IO, p::_LazyProgress)
    _report!(s::ProgressSink, p::_LazyProgress)

Report a tick of `p`. To an `IO`, the bar is drawn lazily (threshold and redraw interval of
`p`). To a sink nothing is drawn; if the sink is cancelled, `_ProgressCancelled` is thrown.
"""
@inline function _report!(::IO, p::_LazyProgress)
    time() < p.tnext[] && return nothing
    # Drawing prints to an `IO` of unknown type. Behind `invokelatest`, the compiled loops do not
    # depend on `print`, whose methods other packages extend, e.g. Makie: loading them would
    # invalidate every loop with a progress bar.
    Base.invokelatest(_draw!, p)
    return nothing
end

@inline function _report!(s::ProgressSink, ::_LazyProgress)
    s.cancel[] && throw(_ProgressCancelled())
    return nothing
end

@noinline function _draw!(p::_LazyProgress)
    # another thread is drawing, skip this frame
    trylock(p.lock) || return nothing
    try
        t = time()
        t < p.tnext[] && return nothing
        p.tnext[] = t + p.dt
        c = p.count[]
        if isnothing(p.bar)
            # `start` makes the ETA use the rate since the bar appeared
            p.bar = Progress(p.n; desc = p.desc, output = p.output, dt = 0.0, start = c)
        end
        # redraws are throttled by `tnext` already, ProgressMeter must not skip them
        update!(p.bar, c; force = true)
    finally
        unlock(p.lock)
    end
    return nothing
end

"""
    _with_progress(f, p::_LazyProgress)
    _with_progress(f, progress::Bool, n, desc)

Run `f(p)`, which calls [`_tick!`](@ref)`(p)` once per finished item, and return its result.
Afterwards the bar is completed. If `f` throws, e.g. an `InterruptException` from Ctrl-C, the
bar is cancelled before the exception is rethrown, so no half-drawn bar is left behind. In both
cases later ticks draw nothing.

If `p` reports to a `ProgressSink`, `p` is the sink's `current` loop while `f` runs; the
previous loop is restored afterwards, also if `f` throws, so nested loops show the inner loop
while it runs.
"""
function _with_progress(f, p::_LazyProgress)
    previous = _enter!(p.output, p)
    result = try
        f(p)
    catch
        _stop!(p, cancel)
        rethrow()
    finally
        _leave!(p.output, previous)
    end
    _stop!(p, finish!)
    return result
end

function _with_progress(f, progress::Bool, n::Integer, desc::AbstractString)
    return _with_progress(f, _LazyProgress(progress, n, desc))
end

# Register `p` as the running loop of its output, return what `_leave!` restores.
_enter!(::IO, ::_LazyProgress) = nothing
_leave!(::IO, _) = nothing

function _enter!(s::ProgressSink, p::_LazyProgress)
    # a disabled loop (`progress = false`) reports nothing
    p.enabled || return @atomic s.current
    return @atomicswap s.current = p
end

function _leave!(s::ProgressSink, previous)
    @atomic s.current = previous
    return nothing
end

function _stop!(p::_LazyProgress, stop)
    lock(p.lock) do
        p.tnext[] = Inf
        # behind `invokelatest` for the reason given in `_report!`
        isnothing(p.bar) || Base.invokelatest(stop, p.bar)
    end
    return nothing
end

"""
    progress_state(s::ProgressSink)

State of the loop that currently reports to `s`: `nothing` between loops, otherwise
`(; desc, count, n, t0)` with `desc` the description without its trailing `": "` (e.g.
`"Tracing beams"`), `count` finished of `n` items and `t0` the `time()` at which the loop
started. Safe to call from any task while the loop runs.
"""
progress_state(s::ProgressSink) = progress_state(@atomic s.current)
progress_state(::Nothing) = nothing
function progress_state(p::_LazyProgress)
    return (; desc = String(chopsuffix(p.desc, ": ")), count = p.count[], n = p.n, t0 = p.t0)
end

"""
    is_cancelled(e)

`true` if the exception `e` stems from a cancelled `ProgressSink`, also if it is wrapped in a
`TaskFailedException` or a `CompositeException` (as thrown by `Threads.@threads`).
"""
is_cancelled(e) = false
is_cancelled(::_ProgressCancelled) = true
is_cancelled(e::TaskFailedException) = is_cancelled(e.task.result)
is_cancelled(e::CompositeException) = any(is_cancelled, e.exceptions)
