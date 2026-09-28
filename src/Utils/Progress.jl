"""
    _LazyProgress

Thread-safe progress counter for long loops. It draws a `ProgressMeter.Progress` bar only
once the loop has run for `threshold` seconds, see [`get_progress_threshold`](@ref).

[`_tick!`](@ref) costs one atomic increment and one `time()` call per item. The bar is created
lazily and redrawn at most every `dt` seconds by one thread at a time (`trylock`), so an
invisible bar never contends for a lock. Use it through [`_with_progress`](@ref), which
finishes the bar, or cancels it if the loop throws.

# Fields

- `n`: total number of items
- `desc`: description printed in front of the bar
- `enabled`: if `false`, [`_tick!`](@ref) returns after a single branch
- `output`: stream the bar is drawn to
- `count`: number of finished items
- `tnext`: earliest `time()` at which the bar is created or redrawn, `Inf` stops drawing
- `dt`: minimum interval in s between redraws
- `lock`: serializes creating, redrawing and stopping the bar
- `bar`: the ProgressMeter bar, `nothing` until it is first drawn
"""
mutable struct _LazyProgress
    const n::Int
    const desc::String
    const enabled::Bool
    const output::IO
    const count::Threads.Atomic{Int}
    const tnext::Threads.Atomic{Float64}
    const dt::Float64
    const lock::Threads.SpinLock
    bar::Nullable{Progress}
end

"""
    _LazyProgress(n, desc; enabled = true, output = stderr, threshold = get_progress_threshold(), dt = 0.2)

Progress counter for `n` items that draws to `output` once `threshold` seconds have passed.
"""
function _LazyProgress(n::Integer, desc::AbstractString; enabled::Bool = true,
        output::IO = stderr, threshold::Real = get_progress_threshold(), dt::Real = 0.2)
    return _LazyProgress(n, desc, enabled, output, Threads.Atomic{Int}(0),
        Threads.Atomic{Float64}(time() + threshold), dt, Threads.SpinLock(), nothing)
end

"""
    _LazyProgress(progress::Bool, n, desc)

Progress counter for the `progress` keyword of a public function. It is enabled only if
`progress` is `true` and `stderr` is a terminal (`Base.TTY`), so documentation builds, CI logs
and piped output stay clean.
"""
function _LazyProgress(progress::Bool, n::Integer, desc::AbstractString)
    return _LazyProgress(n, desc; enabled = progress && stderr isa Base.TTY)
end

"""
    _tick!(p::_LazyProgress)

Count one finished item. Safe to call from any thread.
"""
@inline function _tick!(p::_LazyProgress)
    p.enabled || return nothing
    Threads.atomic_add!(p.count, 1)
    time() < p.tnext[] && return nothing
    return _draw!(p)
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
"""
function _with_progress(f, p::_LazyProgress)
    result = try
        f(p)
    catch
        _stop!(p, cancel)
        rethrow()
    end
    _stop!(p, finish!)
    return result
end

function _with_progress(f, progress::Bool, n::Integer, desc::AbstractString)
    return _with_progress(f, _LazyProgress(progress, n, desc))
end

function _stop!(p::_LazyProgress, stop)
    lock(p.lock) do
        p.tnext[] = Inf
        isnothing(p.bar) || stop(p.bar)
    end
    return nothing
end
