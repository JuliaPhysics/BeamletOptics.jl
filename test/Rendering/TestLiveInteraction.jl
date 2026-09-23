module TestLiveInteraction

using BeamletOptics
using Makie
using GeometryBasics
using Test
using LinearAlgebra

const BMO = BeamletOptics

@testset "Kinematic controls" begin
    Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)
    @test !isnothing(Ext)

    @testset "_ray_plane_intersect" begin
        # Straight down the -z axis onto the z=0 plane
        p = Ext._ray_plane_intersect([0.0, 0, 5], [0.0, 0, -1], [0.0, 0, 0], [0.0, 0, 1])
        @test p ≈ [0.0, 0, 0]

        # Off-axis ray/plane
        p = Ext._ray_plane_intersect([1.0, 2, 5], [0.0, 0, -1], [0.0, 0, 0], [0.0, 0, 1])
        @test p ≈ [1.0, 2, 0]

        # Parallel to the plane: no intersection
        @test isnothing(Ext._ray_plane_intersect([0.0, 0, 5], [1.0, 0, 0], [0.0, 0, 0], [0.0, 0, 1]))

        # Intersection behind the ray origin
        @test isnothing(Ext._ray_plane_intersect([0.0, 0, -5], [0.0, 0, -1], [0.0, 0, 0], [0.0, 0, 1]))
    end

    @testset "_axis_angle_from_rotmatrix round-trips through rotate3d" begin
        cases = [
            ([0.0, 0, 1], deg2rad(30)),
            ([1.0, 0, 0], deg2rad(90)),
            (normalize([1.0, 2, 3]), deg2rad(179)), # near pi
            ([0.0, 1, 0], 1e-9),                     # near 0
            (normalize([1.0, 1, 1]), deg2rad(120)),
        ]
        for (axis, θ) in cases
            R = BMO.rotate3d(axis, θ)
            got_axis, got_angle = Ext._axis_angle_from_rotmatrix(R)
            R2 = BMO.rotate3d(got_axis, got_angle)
            @test isapprox(R2, R; atol = 1e-6)
        end
    end

    @testset "reset restores the original pose" begin
        mir = RoundPlanoMirror(0.025, 0.005)
        P0, R0 = BMO.position(mir), Matrix{Float64}(BMO.orientation(mir))

        translate3d!(mir, [0.03, -0.01, 0.02])
        zrotate3d!(mir, deg2rad(37))
        xrotate3d!(mir, deg2rad(-12))

        R = Matrix{Float64}(BMO.orientation(mir))
        Rd = R0 * R'
        axis, angle = Ext._axis_angle_from_rotmatrix(Rd)
        angle > 1e-12 && rotate3d!(mir, axis, angle)
        translate_to3d!(mir, P0)

        @test isapprox(collect(BMO.position(mir)), collect(P0); atol = 1e-9)
        @test isapprox(Matrix{Float64}(BMO.orientation(mir)), R0; atol = 1e-6)
    end

    @testset "_bbox_wireframe" begin
        bb = GeometryBasics.Rect3d(GeometryBasics.Point3d(0, 0, 0), GeometryBasics.Vec3d(1, 2, 3))
        pts = Ext._bbox_wireframe(bb)
        @test length(pts) == 24 # 12 segments
        xs = [p[1] for p in pts]
        @test isapprox(minimum(xs), 0.0; atol = 1e-6)
        @test isapprox(maximum(xs), 1.0; atol = 1e-6)
    end

    @testset "_ray_pick pure helper" begin
        origin = [0.0, -1.0, 0.0]
        dir = [0.0, 1.0, 0.0]

        @testset "occluder in front is skipped, mirror behind is picked" begin
            cube = BMO.CubeMesh(1)
            translate3d!(cube, -[0.5, 0.5, 0.5])
            translate3d!(cube, [0.0, -0.5, 0.0]) # in front of the mirror, on the ray
            occluder = NonInteractableObject(cube) # not part of any system
            mir = RoundPlanoMirror(0.025, 0.005) # at the origin, behind the occluder
            @test Ext._ray_pick([occluder, mir], origin, dir) === mir
        end

        @testset "nearest of two movables on the ray is picked" begin
            near = RoundPlanoMirror(0.025, 0.005)
            far = RoundPlanoMirror(0.025, 0.005)
            translate3d!(far, [0.0, 2.0, 0.0])
            @test Ext._ray_pick([far, near], origin, dir) === near
            @test Ext._ray_pick([near, far], origin, dir) === near
        end

        @testset "ray misses all movables" begin
            m1 = RoundPlanoMirror(0.025, 0.005)
            m2 = RoundPlanoMirror(0.025, 0.005)
            translate3d!(m2, [0.0, 2.0, 0.0])
            @test isnothing(Ext._ray_pick([m1, m2], [0.1, -1.0, 0.0], dir))
        end

        @testset "objects whose intersect3d errors are skipped, not rethrown" begin
            struct _BrokenRayPickObject <: BMO.AbstractObject{Float64} end
            mir = RoundPlanoMirror(0.025, 0.005)
            @test Ext._ray_pick([_BrokenRayPickObject(), mir], origin, dir) === mir
        end
    end

    # Event handling, picking requires a screen and is replaced
    function _fixture()
        m1 = RoundPlanoMirror(0.025, 0.005)
        m2 = RoundPlanoMirror(0.025, 0.005)
        translate3d!(m2, [0.2, 0, 0])
        sys = System([m1, m2])
        fig = Figure()
        ax = LScene(fig[1, 1])
        h = live_render!(ax, sys)
        return fig, ax, h, m1, m2
    end

    @testset "default pick (nothing) selects via ray casting" begin
        fig, ax, h, m1, m2 = _fixture()
        scene = ax.scene
        # Default camera looks at [0,0,0] (m1's position) from the center of the viewport
        vp = scene.viewport[]
        cx, cy = vp.origin[1] + vp.widths[1] / 2, vp.origin[2] + vp.widths[2] / 2
        ctrl = Ext.kinematic_controls!(ax, h; throttle = false) # no pick kwarg: ray picking
        events(scene).mouseposition[] = (cx, cy)
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        @test ctrl.selected[] === m1
        close(ctrl)
    end

    @testset "ray miss falls back to Makie.pick" begin
        fig, ax, h, m1, m2 = _fixture()
        scene = ax.scene
        ctrl = Ext.kinematic_controls!(ax, h; throttle = false)
        # Mouse position at the corner: the ray misses both mirrors, the fallback finds no plot
        # without a backend
        events(scene).mouseposition[] = (1.0, 1.0)
        @test isnothing(Ext._ray_pick(ctrl, scene))
        @test_logs events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        @test ctrl.selected[] === nothing
        @test !ctrl.dragging
        close(ctrl)
    end

    @testset "grab, drag-translate, release" begin
        fig, ax, h, m1, m2 = _fixture()
        scene = ax.scene
        P0 = collect(Float64.(BMO.position(m1)))

        pick_m1 = ax2 -> (h.handles[1].plots[1], 0)
        ctrl = Ext.kinematic_controls!(ax, h; throttle = false, pick = pick_m1)
        @test ctrl isa Ext.KinematicController

        events(scene).mouseposition[] = (100.0, 100.0)
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        @test ctrl.selected[] === m1
        @test ctrl.dragging

        events(scene).mouseposition[] = (140.0, 160.0)
        @test collect(Float64.(BMO.position(m1))) != P0 # object followed the drag
        @test length(h.handles[1].plots) > 0 # no plot churn

        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release)
        @test !ctrl.dragging
        @test ctrl.selected[] === m1 # stays selected after releasing

        close(ctrl)
    end

    @testset "grab consumes the press (camera does not see it)" begin
        fig, ax, h, m1, m2 = _fixture()
        scene = ax.scene
        pick_m1 = ax2 -> (h.handles[1].plots[1], 0)
        ctrl = Ext.kinematic_controls!(ax, h; throttle = false, pick = pick_m1)

        probe = Ref(0)
        on(events(scene).mousebutton, priority = -1000) do event
            probe[] += 1
            return Consume(false)
        end

        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        @test probe[] == 0 # low-priority listener never reached: event was consumed

        close(ctrl)
    end

    @testset "background click without drag deselects" begin
        fig, ax, h, m1, m2 = _fixture()
        scene = ax.scene
        pick_none = ax2 -> (nothing, 0)
        ctrl = Ext.kinematic_controls!(ax, h; throttle = false, pick = pick_none)
        ctrl.selected[] = m1 # pretend something was already selected

        events(scene).mouseposition[] = (50.0, 50.0)
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        @test ctrl.selected[] === m1 # press alone doesn't deselect yet
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release)
        @test ctrl.selected[] === nothing

        close(ctrl)
    end

    @testset "keyboard fine controls" begin
        fig, ax, h, m1, m2 = _fixture()
        scene = ax.scene
        pick_m1 = ax2 -> (h.handles[1].plots[1], 0)
        ctrl = Ext.kinematic_controls!(ax, h; throttle = false, pick = pick_m1, fine_step = 1e-3, fine_angle = 1e-3)
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release)
        @test ctrl.selected[] === m1

        P0 = collect(Float64.(BMO.position(m1)))
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.up, Keyboard.press)
        P1 = collect(Float64.(BMO.position(m1)))
        @test isapprox(norm(P1 .- P0), 1e-3; atol = 1e-9)

        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.down, Keyboard.press)
        P2 = collect(Float64.(BMO.position(m1)))
        @test isapprox(P2, P0; atol = 1e-9)

        # move mode: ← moves along -x, page up along the rotation axis
        R0 = Matrix{Float64}(BMO.orientation(m1))
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.left, Keyboard.press)
        @test collect(Float64.(BMO.position(m1))) ≈ P0 .- 1e-3 .* R0[:, 1]
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.page_up, Keyboard.press)
        @test collect(Float64.(BMO.position(m1))) ≈ P0 .- 1e-3 .* R0[:, 1] .+ [0, 0, 1e-3]
        @test Matrix{Float64}(BMO.orientation(m1)) ≈ R0

        # rotate mode: ← rotates positively around the rotation axis, position is kept
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.m, Keyboard.press)
        @test ctrl.mode[] == :rotate
        P1 = collect(Float64.(BMO.position(m1)))
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.left, Keyboard.press)
        @test Matrix{Float64}(BMO.orientation(m1)) ≈ BMO.rotate3d([0, 0, 1], 1e-3) * R0
        @test collect(Float64.(BMO.position(m1))) ≈ P1
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.m, Keyboard.press)
        @test ctrl.mode[] == :move

        # shift multiplies the step by 10
        push!(events(scene).keyboardstate, Keyboard.left_shift)
        P3 = collect(Float64.(BMO.position(m1)))
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.up, Keyboard.press)
        P4 = collect(Float64.(BMO.position(m1)))
        @test isapprox(norm(P4 .- P3), 1e-2; atol = 1e-8)
        delete!(events(scene).keyboardstate, Keyboard.left_shift)

        # reset
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.r, Keyboard.press)
        @test isapprox(collect(Float64.(BMO.position(m1))), P0; atol = 1e-9)
        @test isapprox(Matrix{Float64}(BMO.orientation(m1)), R0; atol = 1e-6)

        # escape deselects
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.escape, Keyboard.press)
        @test ctrl.selected[] === nothing

        close(ctrl)
    end

    @testset "on_change and update_render! are called" begin
        fig, ax, h, m1, m2 = _fixture()
        scene = ax.scene
        pick_m1 = ax2 -> (h.handles[1].plots[1], 0)
        changed = Ref{Any}(nothing)
        ctrl = Ext.kinematic_controls!(
            ax, h; throttle = false, pick = pick_m1, on_change = o -> (changed[] = o)
        )
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.up, Keyboard.press)
        @test changed[] === m1
        close(ctrl)
    end

    @testset "errors in on_change are logged once, not rethrown" begin
        fig, ax, h, m1, m2 = _fixture()
        scene = ax.scene
        pick_m1 = ax2 -> (h.handles[1].plots[1], 0)
        ctrl = Ext.kinematic_controls!(
            ax, h; throttle = false, pick = pick_m1, on_change = o -> error("no hits")
        )
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        @test_logs (:error, r"on_change") events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.up, Keyboard.press)
        # the same error again is not logged a second time
        @test_logs events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.up, Keyboard.press)
        close(ctrl)
    end

    @testset "throttle coalesces updates to one per tick" begin
        fig, ax, h, m1, m2 = _fixture()
        scene = ax.scene
        pick_m1 = ax2 -> (h.handles[1].plots[1], 0)
        n_updates = Ref(0)
        ctrl = Ext.kinematic_controls!(
            ax, h; throttle = true, pick = pick_m1, on_change = o -> (n_updates[] += 1)
        )
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        for _ in 1:5
            events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.up, Keyboard.press)
        end
        @test n_updates[] == 0 # nothing applied yet, only marked dirty
        @test ctrl.dirty
        notify(events(scene).tick)
        @test n_updates[] == 1
        @test !ctrl.dirty
        close(ctrl)
    end

    @testset "left-drag rotates in the rotate mode" begin
        fig, ax, h, m1, m2 = _fixture()
        scene = ax.scene
        pick_m1 = ax2 -> (h.handles[1].plots[1], 0)
        ctrl = Ext.kinematic_controls!(ax, h; throttle = false, pick = pick_m1, mode = :rotate,
            rotate_speed = 1e-2)
        P0 = collect(Float64.(BMO.position(m1)))
        R0 = Matrix{Float64}(BMO.orientation(m1))
        events(scene).mouseposition[] = (100.0, 100.0)
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        events(scene).mouseposition[] = (110.0, 100.0)
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release)
        @test Matrix{Float64}(BMO.orientation(m1)) ≈ BMO.rotate3d([0, 0, 1], 0.1) * R0
        @test collect(Float64.(BMO.position(m1))) ≈ P0
        # rings are shown instead of arrows
        @test ctrl.plots[3].visible[]
        close(ctrl)
        @test_throws ArgumentError Ext.kinematic_controls!(ax, h; mode = :fly)
    end

    @testset "h toggles the controls overlay" begin
        fig, ax, h, m1, m2 = _fixture()
        scene = ax.scene
        ctrl = Ext.kinematic_controls!(ax, h; throttle = false, fine_step = 20e-9)
        hint = Ext._help_hint(:move, 20e-9, 10e-6)
        @test ctrl.help_obs[] == hint
        # works without a selected object
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.h, Keyboard.press)
        @test occursin("step: 20 nm", ctrl.help_obs[])
        @test occursin("+/-: change step", ctrl.help_obs[])
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.h, Keyboard.press)
        @test ctrl.help_obs[] == hint
        close(ctrl)
        ctrl = Ext.kinematic_controls!(ax, h; show_help = true)
        @test ctrl.help_obs[] != Ext._help_hint(:move, 10e-9, 10e-6)
        close(ctrl)
    end

    @testset "step size keys" begin
        @test Ext._next_step(3e-8, 1) ≈ 5e-8
        @test Ext._next_step(3e-8, -1) ≈ 2e-8
        @test Ext._next_step(1e-8, 1) ≈ 2e-8
        @test Ext._next_step(5e-8, 1) ≈ 1e-7
        @test Ext._next_step(1e-7, -1) ≈ 5e-8
        @test Ext._next_step(1.0, -1) ≈ 0.5

        fig, ax, h, m1, m2 = _fixture()
        scene = ax.scene
        ctrl = Ext.kinematic_controls!(ax, h; throttle = false, fine_step = 10e-9, fine_angle = 10e-6)
        @test ctrl.help_obs[] == "move mode, step 10 nm, +/-: step, m: switch mode, h: show controls"

        # move mode, no object selected: 1-2-5 sequence on fine_step only
        @test isnothing(ctrl.selected[])
        steps = Float64[]
        for c in ('+', '+', '+', '-')
            events(scene).unicode_input[] = c
            push!(steps, ctrl.fine_step)
        end
        @test steps ≈ [20e-9, 50e-9, 100e-9, 50e-9]
        @test ctrl.fine_angle == 10e-6
        @test ctrl.help_obs[] == "move mode, step 50 nm, +/-: step, m: switch mode, h: show controls"

        # help overlay shows the new step as well
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.h, Keyboard.press)
        @test occursin("step: 50 nm", ctrl.help_obs[])
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.h, Keyboard.press)

        # rotate mode: fine_angle only
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.m, Keyboard.press)
        steps = Float64[]
        for c in ('+', '+', '+', '-')
            events(scene).unicode_input[] = c
            push!(steps, ctrl.fine_angle)
        end
        @test steps ≈ [20e-6, 50e-6, 100e-6, 50e-6]
        @test ctrl.fine_step ≈ 50e-9
        @test occursin("rotate mode, step 50 µrad", ctrl.help_obs[])

        # clamping at the bounds
        ctrl.fine_angle = π / 4
        events(scene).unicode_input[] = '+'
        @test ctrl.fine_angle ≈ π / 4
        events(scene).unicode_input[] = '-'
        @test ctrl.fine_angle ≈ 0.5
        ctrl.fine_angle = 1e-9
        events(scene).unicode_input[] = '-'
        @test ctrl.fine_angle ≈ 1e-9
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.m, Keyboard.press)
        ctrl.fine_step = 1.0
        events(scene).unicode_input[] = '+'
        @test ctrl.fine_step ≈ 1.0
        ctrl.fine_step = 1e-12
        events(scene).unicode_input[] = '-'
        @test ctrl.fine_step ≈ 1e-12

        # + and - are consumed, other characters are passed on
        probe = Char[]
        on(events(scene).unicode_input, priority = -1000) do c
            push!(probe, c)
            return Consume(false)
        end
        events(scene).unicode_input[] = '+'
        events(scene).unicode_input[] = 'x'
        @test probe == ['x']

        # also with a selected object
        ctrl.selected[] = m1
        ctrl.fine_step = 10e-9
        events(scene).unicode_input[] = '+'
        @test ctrl.fine_step ≈ 20e-9

        # the listener is removed by close
        close(ctrl)
        events(scene).unicode_input[] = '+'
        @test ctrl.fine_step ≈ 20e-9
    end

    @testset "close disconnects listeners and removes the selection box" begin
        fig, ax, h, m1, m2 = _fixture()
        scene = ax.scene
        pick_m1 = ax2 -> (h.handles[1].plots[1], 0)
        n0 = length(ax.scene.plots)
        ctrl = Ext.kinematic_controls!(ax, h; throttle = false, pick = pick_m1)
        @test length(ax.scene.plots) == n0 + 5 # selection box, gizmo and controls overlay

        close(ctrl)
        @test isempty(ctrl.listeners)
        @test length(ax.scene.plots) == n0

        # listeners are gone: this must not error and must not select anything
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        @test ctrl.selected[] === nothing
    end
end

end # module
