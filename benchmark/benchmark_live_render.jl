# Benchmark: legacy `render!` redraw vs. `live_render!` + `update_render!`
#
# Each frame moves a component, re-solves the system, resyncs the plots and renders an image
# via `colorbuffer`. Run with a GLMakie-enabled environment:
#
#   julia --project=<env with BeamletOptics + GLMakie> benchmark/benchmark_live_render.jl

using BeamletOptics, GLMakie, Statistics, Printf

GLMakie.activate!(visible = false)

const N_WARMUP = 2
const N_RUNS = 10

median_ms(f; runs = N_RUNS) = (for _ in 1:N_WARMUP; f(); end; 1e3 * median([@elapsed(f()) for _ in 1:runs]))

"""Legacy frame: delete all plots, render system and beam from scratch."""
function legacy_frame!(fig, ax, system, beam; kwargs...)
    for p in collect(ax.scene.plots)
        delete!(ax.scene, p)
    end
    render!(ax, system)
    render!(ax, beam; kwargs...)
    colorbuffer(fig)
end

"""Legacy frame (lower bound): components drawn once and left in place, only the beam is redrawn."""
function legacy_rays_frame!(fig, ax, n0, beam; kwargs...)
    for p in ax.scene.plots[(n0 + 1):end]
        delete!(ax.scene, p)
    end
    render!(ax, beam; kwargs...)
    colorbuffer(fig)
end

function live_frame!(fig, hsys, hbeam)
    update_render!(hsys)
    update_render!(hbeam)
    colorbuffer(fig)
end

function run_scenario(name, setup, step!; beam_kwargs = (;))
    system, beam = setup()
    step!(system, beam)

    t_solve = median_ms(() -> step!(system, beam))
    println("\n", name)

    # legacy, full redraw
    fig = Figure(); ax = LScene(fig[1, 1])
    render!(ax, system); render!(ax, beam; beam_kwargs...)
    n_legacy = length(ax.scene.plots)
    t_legacy = median_ms(runs = 3) do
        step!(system, beam)
        legacy_frame!(fig, ax, system, beam; beam_kwargs...)
    end
    @printf "  %-44s %9.2f ms   (%d plots)\n" "legacy: full redraw + frame" t_legacy n_legacy
    flush(stdout)

    # legacy, beam only (visually wrong: components don't follow, lower bound)
    fig = Figure(); ax = LScene(fig[1, 1])
    render!(ax, system); n0 = length(ax.scene.plots)
    render!(ax, beam; beam_kwargs...)
    t_legacy_rays = median_ms(runs = 3) do
        step!(system, beam)
        legacy_rays_frame!(fig, ax, n0, beam; beam_kwargs...)
    end

    # live
    fig = Figure(); ax = LScene(fig[1, 1])
    hsys = live_render!(ax, system)
    hbeam = live_render!(ax, beam; beam_kwargs...)
    n_live = length(ax.scene.plots)
    t_live = median_ms() do
        step!(system, beam)
        live_frame!(fig, hsys, hbeam)
    end
    t_update = median_ms() do
        step!(system, beam)
        update_render!(hsys); update_render!(hbeam)
    end
    t_frame = median_ms(() -> colorbuffer(fig))

    @printf "  %-44s %9.2f ms\n" "solve_system! (incl. kinematics)" t_solve
    @printf "  %-44s %9.2f ms\n" "legacy: beam-only redraw + frame (wrong)" t_legacy_rays
    @printf "  %-44s %9.2f ms   (%d plots)\n" "live: update_render! + frame" t_live n_live
    @printf "  %-44s %9.2f ms\n" "live: solve + update_render!, no frame" t_update
    @printf "  %-44s %9.2f ms\n" "bare frame (colorbuffer, no changes)" t_frame
    @printf "  %-44s %9.1fx\n" "speedup full redraw -> live" t_legacy / t_live
end

# Scenario 1: lens + mirror, 200-ray collimated source
function lens_setup()
    lens = SphericalLens(50e-3, -50e-3, 10e-3, 25e-3, 1.5168)
    mir = RoundPlanoMirror(25.4e-3, 5e-3)
    translate3d!(mir, [0, 60e-3, 0]); zrotate3d!(mir, deg2rad(45))
    system = System([lens, mir])
    src = CollimatedSource([0.0, -20e-3, 0.0], [0.0, 1.0, 0.0], 12e-3; num_rings = 5, num_rays = 200)
    return system, src
end
lens_step!(system, src) = (zrotate3d!(system.objects[2], 1e-4); solve_system!(system, src))

# Scenario 2: Michelson interferometer with a Gaussian beamlet
function michelson_setup()
    l_0 = 0.1
    m1 = SquarePlanoMirror2D(BeamletOptics.inch)
    m2 = SquarePlanoMirror2D(BeamletOptics.inch)
    bs = ThinBeamsplitter(BeamletOptics.inch, reflectance = 0.5)
    pd = Detector(BeamletOptics.inch / 5)
    translate3d!(m1, [l_0, 0, 0]); translate3d!(m2, [0, l_0, 0]); translate3d!(pd, [-l_0, 0, 0])
    zrotate3d!(bs, deg2rad(45)); zrotate3d!(m1, deg2rad(90)); zrotate3d!(pd, deg2rad(90))
    system = System([m1, m2, bs, pd])
    beam = GaussianBeamlet([0, -l_0, 0], [0, 1.0, 0], 632.8e-9, 1e-3)
    return system, beam
end
function michelson_step!(system, beam)
    zrotate3d!(system.objects[1], 1e-6)
    empty!(system.objects[4])
    solve_system!(system, beam)
end

run_scenario("Lens + mirror, CollimatedSource (200 rays)", lens_setup, lens_step!; beam_kwargs = (; render_every = 1))
run_scenario("Michelson, GaussianBeamlet", michelson_setup, michelson_step!)
