#=
Hit buffers

The tasks that solve the beams of a group push hits concurrently. Under the lock of a detector each
push waits for the others, e.g. 8 threads were slower than 4 for a detector only. Hence each task
collects its hits and appends them once per detector at its end.
=#

"""
The hits of the current task per [`Detector`](@ref), see [`with_hit_buffers`](@ref). Pushing
`detector => hit` collects the `hit` in the buffer of the `detector`. The buffer of the last
detector is kept at hand, since consecutive hits mostly go to the same detector.
"""
mutable struct HitBuffers
    const buffers::IdDict{Detector, Vector}
    last_detector::Union{Nothing, Detector}
    last_buffer::Vector
    HitBuffers() = new(IdDict{Detector, Vector}(), nothing, [])
end

# Key of the buffers of the current task in the task-local storage
const HIT_BUFFERS_KEY = :BeamletOptics_hit_buffers

"Stores each pushed `detector => hit` directly in its detector, guarded by the `lock` of the detector."
struct LockedHits end

# Where the current task pushes its hits: its buffers within `with_hit_buffers`, otherwise the detectors
current_hit_sink() = get(task_local_storage(), HIT_BUFFERS_KEY, LockedHits())::Union{LockedHits, HitBuffers}

"""
    with_hit_buffers(f)

Calls `f()` such that the hits that the current task pushes to a [`Detector`](@ref) are collected
in a buffer of the task, and appends them to their detectors when `f` returns or throws, one
`append!` per detector. The solve of a beam group calls it in each of its tasks, which then do not
wait for the lock of a detector at every hit. The hits are stored when `with_hit_buffers` returns,
not while `f` runs.
"""
function with_hit_buffers(f)
    buffers = HitBuffers()
    try
        return task_local_storage(f, HIT_BUFFERS_KEY, buffers)
    finally
        for (d, new) in buffers.buffers
            append!(d, new)
        end
    end
end

function push!(::LockedHits, (d, new)::Pair{<:Detector, <:AbstractDetectorHit})
    lock(d.lock) do
        # set data type on first push
        if isnothing(hits(d))
            hits!(d, [new])
            return nothing
        end
        # if the hit type does not match, throws a conversion error
        push!(hits(d), new)
        return nothing
    end
end

function push!(b::HitBuffers, (d, new)::Pair{<:Detector, H}) where {H <: AbstractDetectorHit}
    if b.last_detector !== d
        b.last_buffer = get!(() -> H[], b.buffers, d)
        b.last_detector = d
    end
    push!(b.last_buffer::Vector{H}, new)
    return nothing
end
