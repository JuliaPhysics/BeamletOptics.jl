# Benchmark: responsiveness of `live_view` with a heavy scene
#
# Synthetic scene (no external files, fully deterministic): a `NonInteractableObject` with a
# ~1 M triangle mesh, 10 spherical lenses, 4 detectors that are all hit and a 1000-ray
# `UniformDiscSource`. The live view is shown in a GLMakie window of 1600×950 and measured for
#
# - frame time: 60 camera rotation steps (6° about the up axis each), see `frame!` for what one
#   frame includes;
# - idle redraws: frames actually rendered during 2 s without input, see `count_redraws`;
# - solve + panel time: full re-solve incl. the detector panels after moving one lens, see
#   `solve_times`.
#
# Run with the docs environment (it has GLMakie):
#
#   julia --project=docs benchmark/benchmark_live_view.jl [--layout compact|app] [--out file.toml]
#
# `--layout` is passed to `live_view` only if it is not `compact`, i.e. the script also works on
# versions without the `layout` kwarg.

using BeamletOptics, GLMakie, Statistics, Printf, LinearAlgebra

const BMO = BeamletOptics
const BME = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)
const mm = 1e-3

const WINDOW_SIZE = (1600, 950)
const N_WARMUP_FRAMES = 10
const N_FRAMES = 60
const IDLE_SECONDS = 2.0
const N_WARMUP_SOLVES = 2
const N_SOLVES = 10

#=
Arguments
=#

function parse_args(args)
    layout, out = :compact, nothing
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--layout" && i < length(args)
            layout = Symbol(args[i + 1]); i += 2
        elseif startswith(a, "--layout=")
            layout = Symbol(split(a, '='; limit = 2)[2]); i += 1
        elseif a == "--out" && i < length(args)
            out = args[i + 1]; i += 2
        elseif startswith(a, "--out=")
            out = split(a, '='; limit = 2)[2]; i += 1
        else
            error("unknown argument $a, use [--layout compact|app] [--out file.toml]")
        end
    end
    layout in (:compact, :app) || error("--layout must be compact or app, got $layout")
    return layout, out
end

#=
Scene
=#

"""
UV sphere of radius `r` [m] around `center` with `nθ` bands and `nφ` segments as a
`BeamletOptics.Mesh`, `2 nθ nφ - 2 nφ` triangles (one triangle per segment at the poles).
"""
function sphere_mesh(center, r; nθ = 500, nφ = 1000)
    nv = (nθ + 1) * nφ
    V = Matrix{Float64}(undef, nv, 3)
    idx(i, j) = i * nφ + mod(j, nφ) + 1
    for i in 0:nθ, j in 0:(nφ - 1)
        θ, φ = π * i / nθ, 2π * j / nφ
        V[idx(i, j), :] = center .+ r .* (sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ))
    end
    F = Matrix{Int}(undef, 2nθ * nφ - 2nφ, 3)
    k = 0
    for i in 0:(nθ - 1), j in 0:(nφ - 1)
        a, b, c, d = idx(i, j), idx(i, j + 1), idx(i + 1, j), idx(i + 1, j + 1)
        # the first and last band touch the poles, where one of the two triangles is degenerate
        i > 0 && (F[k += 1, :] = [a, c, b])
        i < nθ - 1 && (F[k += 1, :] = [b, c, d])
    end
    @assert k == size(F, 1)
    return BMO.Mesh{Float64}(V, F, Matrix{Float64}(I, 3, 3), BMO.Point3{Float64}(center...), 1.0)
end

"""
Heavy test scene along the optical axis +y: 5 pairs of a biconvex and a biconcave lens
(f ≈ ±100 mm), 3 transparent detectors between them and a stopping detector at the end, a
tessellated sphere (radius 40 mm) next to the beam path and a 1000-ray disc source.
"""
function heavy_scene()
    lenses = BMO.AbstractObject[]
    for k in 1:10
        lens = isodd(k) ? SphericalLens(100mm, -100mm, 4mm, 25.4mm, 1.5) :
               SphericalLens(-100mm, 100mm, 2mm, 25.4mm, 1.5)
        translate3d!(lens, [0, 25mm * (k - 1), 0])
        push!(lenses, lens)
    end
    detectors = [Detector(25mm, false), Detector(25mm, false), Detector(25mm, false), Detector(25mm)]
    for (pd, y) in zip(detectors, (64mm, 139mm, 214mm, 260mm))
        translate3d!(pd, [0, y, 0])
    end
    housing = NonInteractableObject(sphere_mesh([80mm, 130mm, 0.0], 40mm))
    system = System([lenses; detectors; housing])
    source = UniformDiscSource([0, -30mm, 0], [0, 1, 0], 8mm, 1064e-9; num_rays = 1000)
    labels = Dict{Any, String}()
    foreach(((k, l),) -> labels[l] = "L$k", enumerate(lenses))
    foreach(((k, d),) -> labels[d] = "PD$k", enumerate(detectors))
    labels[housing] = "housing"
    return (; system, source, lenses, detectors, housing, labels)
end

#=
Measurements
=#

"""
One full frame, the loop body of `GLMakie.render_asap` plus a GPU sync:

- `pollevents(screen, RegularRenderTick)`: GLFW events and the `events(fig).tick` listeners, i.e.
  the per-frame CPU work of the live view (which `colorbuffer` would skip, it sends a
  `BackendTick` that does not reach `events.tick`);
- `poll_updates`: resolves the changed render objects (compute graph → GPU buffers);
- `render_frame` + `SwapBuffers`: draws the figure;
- `glFinish`: waits until the GPU is done, so the time is not just the CPU submission.

The render loop of the screen is stopped meanwhile, such that it does not render in between.
"""
function frame!(screen)
    GLMakie.pollevents(screen, Makie.RegularRenderTick)
    GLMakie.poll_updates(screen)
    GLMakie.render_frame(screen)
    GLMakie.GLFW.SwapBuffers(GLMakie.to_native(screen))
    GLMakie.ModernGL.glFinish()
    return nothing
end

"""Rotates the camera of the 3D view by `angle` [rad] about its up vector through `lookat`."""
function orbit!(gui, angle)
    cam = Makie.cameracontrols(gui.ax.scene)
    eye, lookat, up = Vector{Float64}(cam.eyeposition[]), Vector{Float64}(cam.lookat[]),
                      normalize(Vector{Float64}(cam.upvector[]))
    # Rodrigues rotation of the eye offset about `up`
    v = eye .- lookat
    v = v .* cos(angle) .+ cross(up, v) .* sin(angle) .+ up .* dot(up, v) .* (1 - cos(angle))
    set_view(gui.ax, lookat .+ v, lookat, up)
    return nothing
end

"""Frame times [ms] of `n` camera rotation steps of 6°, i.e. a full turn for `n = 60`."""
function frame_times(gui, screen; n = N_FRAMES, warmup = N_WARMUP_FRAMES)
    GLMakie.stop_renderloop!(screen; close_after_renderloop = false)
    try
        for _ in 1:warmup
            orbit!(gui, 2π / n)
            frame!(screen)
        end
        return map(1:n) do _
            t0 = time_ns()
            orbit!(gui, 2π / n)
            frame!(screen)
            1e-6 * (time_ns() - t0)
        end
    finally
        GLMakie.start_renderloop!(screen)
    end
end

"""
Counts the ticks of the render loop during `seconds` and how many of them follow a rendered frame.
GLMakie renders on demand: each loop iteration sends a tick with the state of the previous
iteration, `RegularRenderTick` if it rendered a frame, `SkippedRenderTick` if nothing changed. With
`control = true`, the status line is changed every 100 ms meanwhile, which must cause about
10 redraws per second, i.e. shows that the counter sees rendered frames.
"""
function count_redraws(gui; seconds = IDLE_SECONDS, control = false)
    rendered, ticks = Ref(0), Ref(0)
    listener = on(events(gui.fig).tick) do tick
        ticks[] += 1
        tick.state == Makie.RegularRenderTick && (rendered[] += 1)
        return nothing
    end
    k = Ref(0)
    timer = control ? Timer(_ -> (gui.status.text[] = "control $(k[] += 1)"), 0.1; interval = 0.1) :
            nothing
    sleep(seconds)
    isnothing(timer) || close(timer)
    off(listener)
    return (; rendered = rendered[], ticks = ticks[], changes = k[])
end

"""
Durations [ms] of full re-solves after moving the first lens by ±10 µm along x: the lens is moved
like the keyboard controls do (`_change!` + `update_render!` of the controls), then `_resolve!`
solves the systems and updates the beams and all detector panels synchronously, as the deferred
full solve does once a movement pauses. `solve` and `panel` are the `solve_time` and `panel_time`
the live view measures itself (see `_apply!`), `total` the wall time of `_resolve!`.
"""
function solve_times(gui, lens; n = N_SOLVES, warmup = N_WARMUP_SOLVES)
    ctrl = gui.controls
    rows = map(1:(warmup + n)) do k
        dx = (isodd(k) ? 10 : -10) * 1e-6
        BME._change!(() -> translate3d!(lens, [dx, 0, 0]), ctrl, lens)
        update_render!(ctrl.h)
        t0 = time_ns()
        BME._resolve!(gui, lens)
        total = 1e-6 * (time_ns() - t0)
        (; solve = 1e3 * gui.solve_time, panel = 1e3 * gui.panel_time, total)
    end
    return rows[(warmup + 1):end]
end

#=
Main
=#

function main(args)
    layout, out = parse_args(args)
    GLMakie.activate!()
    scene = heavy_scene()
    ntri = size(BMO.faces(scene.housing.shape), 1)

    kw = (; size = WINDOW_SIZE, labels = scene.labels)
    t_build = @elapsed gui = layout === :compact ?
                             live_view(scene.system => scene.source; kw...) :
                             live_view(scene.system => scene.source; layout, kw...)
    screen = display(gui)
    # Let the window appear and the first frames (and the home view) settle
    sleep(2.0)
    hits = [isnothing(BMO.hits(pd)) ? 0 : length(BMO.hits(pd)) for pd in scene.detectors]
    any(iszero, hits) && @warn "a detector panel is empty" hits

    ft = frame_times(gui, screen)
    sleep(1.0)
    idle = count_redraws(gui)
    control = count_redraws(gui; control = true)
    st = solve_times(gui, first(scene.lenses))
    sleep(0.5)
    idle_after = count_redraws(gui)

    r = (; layout = String(layout), julia = string(VERSION), threads = Threads.nthreads(),
        window = collect(WINDOW_SIZE), framebuffer = collect(size(screen)),
        mesh_triangles = ntri, detector_hits = hits, build_s = t_build,
        frame_mean_ms = mean(ft), frame_median_ms = median(ft), frame_min_ms = minimum(ft),
        frame_max_ms = maximum(ft),
        idle_redraws = idle.rendered, idle_ticks = idle.ticks,
        idle_redraws_after_solve = idle_after.rendered,
        control_redraws = control.rendered, control_changes = control.changes,
        solve_median_ms = median(getfield.(st, :solve)),
        panel_median_ms = median(getfield.(st, :panel)),
        solve_panel_median_ms = median(getfield.(st, :solve) .+ getfield.(st, :panel)),
        resolve_wall_median_ms = median(getfield.(st, :total)))

    println()
    @printf "live_view benchmark, layout = %s, Julia %s, %d threads\n" r.layout r.julia r.threads
    @printf "  window %d×%d (framebuffer %d×%d), mesh %d triangles, detector hits %s\n" r.window... r.framebuffer... r.mesh_triangles string(hits)
    @printf "  %-38s %10.2f s\n" "build (live_view)" r.build_s
    @printf "  %-38s %10.2f ms\n" "frame time, mean ($N_FRAMES rotation steps)" r.frame_mean_ms
    @printf "  %-38s %10.2f ms\n" "frame time, median" r.frame_median_ms
    @printf "  %-38s %6.2f - %.2f ms\n" "frame time, min - max" r.frame_min_ms r.frame_max_ms
    @printf "  %-38s %10d   (of %d ticks)\n" "idle redraws in $(IDLE_SECONDS) s" r.idle_redraws r.idle_ticks
    @printf "  %-38s %10d\n" "idle redraws in $(IDLE_SECONDS) s after solves" r.idle_redraws_after_solve
    @printf "  %-38s %10d   (%d changes, counter check)\n" "control: redraws, label every 100 ms" r.control_redraws r.control_changes
    @printf "  %-38s %10.2f ms\n" "solve time, median ($N_SOLVES)" r.solve_median_ms
    @printf "  %-38s %10.2f ms\n" "panel time, median" r.panel_median_ms
    @printf "  %-38s %10.2f ms\n" "solve + panel time, median" r.solve_panel_median_ms
    @printf "  %-38s %10.2f ms\n" "_resolve! wall time, median" r.resolve_wall_median_ms
    r.control_redraws > 0 || @warn "the redraw counter saw no frames in the control run"

    isnothing(out) || write_toml(out, r)
    close(gui)
    close(screen)
    return r
end

"""Writes the flat NamedTuple `r` of numbers, strings and vectors as TOML."""
function write_toml(path, r)
    val(x::AbstractString) = repr(String(x))
    val(x::Integer) = string(x)
    val(x::Real) = @sprintf("%.4f", x)
    val(x::AbstractVector) = "[" * join(val.(x), ", ") * "]"
    open(path, "w") do io
        for (k, v) in pairs(r)
            println(io, k, " = ", val(v))
        end
    end
    println("  written to ", path)
    return nothing
end

main(ARGS)
