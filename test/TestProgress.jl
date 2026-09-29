module TestProgress

using BeamletOptics
using Test
using LinearAlgebra
using Base.ScopedValues: with

const BMO = BeamletOptics

const mm = 1e-3

function lens_system()
    lens = SphericalLens(100mm, Inf, 5mm, 25.4mm, x -> 1.5)
    pd = Detector(20mm)
    translate3d!(pd, [0, 100mm, 0])
    return System([lens, pd]), pd
end

ray_source() = UniformDiscSource([0, -10mm, 0], [0, 1, 0], 10mm, 1e-6; num_rays = 500)

ray_path(bg) = [(position(r), BMO.direction(r)) for b in BMO.beams(bg) for r in BMO.rays(b)]

function ticked_loop(p, n)
    return BMO._with_progress(p) do prog
        Threads.@threads for _ in 1:n
            BMO._tick!(prog)
        end
    end
end

"""
    ProbeSystem(inner, sink)

Delegates the per-beam solves of a beam group to `inner` and records the loop that reports to
`sink` (and the time) at each beam.
"""
struct ProbeSystem{S <: BMO.AbstractSystem} <: BMO.AbstractSystem
    inner::S
    sink::BMO.ProgressSink
    seen::Vector{Any}
    times::Vector{Float64}
    lock::ReentrantLock
end
ProbeSystem(inner, sink) = ProbeSystem(inner, sink, Any[], Float64[], ReentrantLock())

function BMO.solve_system!(probe::ProbeSystem, beam::BMO.AbstractBeam; kwargs...)
    lock(probe.lock) do
        push!(probe.seen, @atomic probe.sink.current)
        push!(probe.times, time())
    end
    return solve_system!(probe.inner, beam; kwargs...)
end

@testset "Progress bars" begin
    @testset "_LazyProgress" begin
        @testset "invisible below threshold" begin
            buf = IOBuffer()
            p = BMO._LazyProgress(100, "Test: "; output = buf, threshold = Inf)
            ticked_loop(p, 100)
            @test p.count[] == 100
            @test isnothing(p.bar)
            @test isempty(take!(buf))
        end

        @testset "drawn above threshold" begin
            buf = IOBuffer()
            p = BMO._LazyProgress(100, "Test: "; output = buf, threshold = 0.0)
            ticked_loop(p, 100)
            out = String(take!(buf))
            @test p.count[] == 100
            @test occursin("Test: ", out)
            @test occursin("100%", out)
        end

        @testset "disabled" begin
            buf = IOBuffer()
            p = BMO._LazyProgress(100, "Test: "; enabled = false, output = buf, threshold = 0.0)
            ticked_loop(p, 100)
            @test p.count[] == 0
            @test isempty(take!(buf))
        end

        @testset "cancelled on exception" begin
            buf = IOBuffer()
            p = BMO._LazyProgress(10, "Test: "; output = buf, threshold = 0.0)
            @test_throws ErrorException BMO._with_progress(p) do prog
                BMO._tick!(prog)
                error("loop failed")
            end
            out = String(take!(buf))
            @test occursin("Aborted", out)
            @test endswith(out, "\n")
            # a stopped bar is never redrawn
            BMO._tick!(p)
            @test isempty(take!(buf))
        end

        @testset "result and terminal gate" begin
            @test BMO._with_progress(_ -> 42, false, 1, "Test: ") == 42
            @test !BMO._LazyProgress(false, 10, "Test: ").enabled
            redirect_stderr(devnull) do
                @test !BMO._LazyProgress(true, 10, "Test: ").enabled
            end
        end
    end

    @testset "solve_system! on beam groups" begin
        sys_off, _ = lens_system()
        sys_on, pd_on = lens_system()
        bg_off, bg_on = ray_source(), ray_source()
        solve_system!(sys_off, bg_off; progress = false)
        solve_system!(sys_on, bg_on; progress = true)
        @test ray_path(bg_on) == ray_path(bg_off)
        # threads push hits in varying order, so fields are only compared on the same detector
        @test BMO.hits(pd_on) isa Vector{<:BMO.RayHit}
        @test electric_field(pd_on; progress = true) == electric_field(pd_on; progress = false)
        # `progress` is consumed by the beam group method, other kwargs still reach each beam
        sys, _ = lens_system()
        @test isnothing(solve_system!(sys, ray_source(); progress = true, r_max = 5))
    end

    @testset "electric_field for beamlet and polarized hits" begin
        field_matches(pd) = electric_field(pd; n = 50, progress = true) ==
                            electric_field(pd; n = 50, progress = false)

        @testset "GaussianBeamletHit" begin
            pd = Detector(20mm)
            translate3d!(pd, [0, 50mm, 0])
            solve_system!(System([pd]), GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 1e-6, 1mm))
            @test BMO.hits(pd) isa Vector{<:BMO.GaussianBeamletHit}
            @test field_matches(pd)
        end

        @testset "AstigmaticGaussianBeamletHit" begin
            sys, pd = lens_system()
            bg = SphericalGaussianBeamletSource([0, -100mm, 0], [0, 1, 0], deg2rad(4), 1e-6;
                num_rings = 3, num_rays = 100)
            solve_system!(sys, bg; progress = true)
            @test BMO.hits(pd) isa Vector{<:BMO.AstigmaticGaussianBeamletHit}
            @test field_matches(pd)
        end

        @testset "PolarizedRayHit" begin
            λ = 1e-6
            pd = Detector(1.0)
            for dir in (normalize([0.1, 1, 0]), normalize([-0.1, 1, 0]))
                ray = PolarizedRay(-λ .* dir, dir, λ, [0, 0, 1.0])
                BMO.intersection!(ray, BMO.Intersection(λ, BMO.Point3(-dir)))
                push!(pd, BMO.PolarizedRayHit(ray, BMO.optical_path_length(ray)))
            end
            @test BMO.hits(pd) isa Vector{<:BMO.PolarizedRayHit}
            @test field_matches(pd)
        end
    end

    @testset "progress sink" begin
        @testset "loop state and nesting" begin
            s = BMO.ProgressSink()
            @test isnothing(BMO.progress_state(s))
            outer = BMO._LazyProgress(3, "Outer: "; output = s)
            inner = BMO._LazyProgress(2, "Inner: "; output = s)
            BMO._with_progress(outer) do po
                BMO._tick!(po)
                @test BMO.progress_state(s).desc == "Outer"
                BMO._with_progress(inner) do pi
                    BMO._tick!(pi)
                    st = BMO.progress_state(s)
                    @test st.desc == "Inner"
                    @test (st.count, st.n) == (1, 2)
                    @test st.t0 == inner.t0
                end
                # the outer loop is shown again after the inner one
                @test BMO.progress_state(s).desc == "Outer"
            end
            @test isnothing(BMO.progress_state(s))
            # restored also if the loop throws
            @test_throws ErrorException BMO._with_progress(_ -> error("loop failed"), outer)
            @test isnothing(BMO.progress_state(s))
        end

        @testset "no drawing, ignores the threshold" begin
            s = BMO.ProgressSink()
            p = BMO._LazyProgress(100, "Test: "; output = s, threshold = 0.0)
            ticked_loop(p, 100)
            @test p.count[] == 100
            @test isnothing(p.bar)
        end

        @testset "scoped output and gate" begin
            s = BMO.ProgressSink()
            # the sink is used even if stderr is no terminal
            redirect_stderr(devnull) do
                with(BMO.PROGRESS_SINK => s) do
                    p = BMO._LazyProgress(true, 10, "Test: ")
                    @test p.output === s
                    @test p.enabled
                    @test !BMO._LazyProgress(false, 10, "Test: ").enabled
                end
                @test BMO._LazyProgress(true, 10, "Test: ").output isa IO
            end
        end

        @testset "@inferred _tick!" begin
            p_io = BMO._LazyProgress(10, "Test: "; output = IOBuffer(), threshold = Inf)
            p_sink = BMO._LazyProgress(10, "Test: "; output = BMO.ProgressSink())
            @test isnothing(@inferred BMO._tick!(p_io))
            @test isnothing(@inferred BMO._tick!(p_sink))
        end

        @testset "is_cancelled" begin
            failed(e) = (t = Threads.@spawn throw(e); try wait(t) catch end; t)
            c = BMO._ProgressCancelled()
            @test BMO.is_cancelled(c)
            @test !BMO.is_cancelled(ErrorException("x"))
            @test BMO.is_cancelled(TaskFailedException(failed(c)))
            @test !BMO.is_cancelled(TaskFailedException(failed(ErrorException("x"))))
            @test BMO.is_cancelled(CompositeException([TaskFailedException(failed(c))]))
            @test !BMO.is_cancelled(CompositeException([ErrorException("x")]))
        end
    end

    @testset "solve_system! with a progress sink" begin
        sources() = CollimatedSource([0, -10mm, 0], [0, 1, 0], 10mm, 1e-6;
            num_rings = 2, num_rays = 40)

        @testset "reports all items" begin
            sys, _ = lens_system()
            bg = sources()
            n = length(bg)
            s = BMO.ProgressSink()
            probe = ProbeSystem(sys, s)
            t_start = time()
            stderr_out = mktemp() do path, io
                redirect_stderr(io) do
                    with(BMO.PROGRESS_SINK => s) do
                        solve_system!(probe, bg)
                    end
                end
                close(io)
                read(path, String)
            end
            # nothing is written to the terminal
            @test isempty(stderr_out)
            @test length(probe.seen) == n
            p = only(unique(probe.seen))
            @test p isa BMO._LazyProgress{BMO.ProgressSink}
            @test p.count[] == n
            st = BMO.progress_state(p)
            @test st.desc == "Tracing beams"
            @test st.n == n
            @test t_start <= st.t0 <= minimum(probe.times)
            @test isnothing(BMO.progress_state(s))
            @test all(b -> !isempty(BMO.rays(b)), BMO.beams(bg))
        end

        @testset "progress = false reports nothing" begin
            sys, _ = lens_system()
            bg = sources()
            s = BMO.ProgressSink()
            probe = ProbeSystem(sys, s)
            with(BMO.PROGRESS_SINK => s) do
                solve_system!(probe, bg; progress = false)
            end
            @test length(probe.seen) == length(bg)
            @test all(isnothing, probe.seen)
        end

        @testset "cancelled" begin
            sys, _ = lens_system()
            s = BMO.ProgressSink()
            s.cancel[] = true
            err = try
                with(BMO.PROGRESS_SINK => s) do
                    solve_system!(sys, sources())
                end
                nothing
            catch e
                e
            end
            # thrown through the tasks of `Threads.@threads`
            @test err isa CompositeException
            @test BMO.is_cancelled(err)
            @test isnothing(BMO.progress_state(s))
        end
    end
end

end
