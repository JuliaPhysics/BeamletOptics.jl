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

        R0 = Matrix{Float64}(BMO.orientation(m1))
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.left, Keyboard.press)
        R1 = Matrix{Float64}(BMO.orientation(m1))
        @test !isapprox(R1, R0; atol = 1e-9)

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

    @testset "h toggles the controls overlay" begin
        fig, ax, h, m1, m2 = _fixture()
        scene = ax.scene
        ctrl = Ext.kinematic_controls!(ax, h; throttle = false, fine_step = 20e-9)
        @test ctrl.help_obs[] == Ext._HELP_HINT
        # works without a selected object
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.h, Keyboard.press)
        @test occursin("20.0 nm", ctrl.help_obs[])
        events(scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.h, Keyboard.press)
        @test ctrl.help_obs[] == Ext._HELP_HINT
        close(ctrl)
        ctrl = Ext.kinematic_controls!(ax, h; show_help = true)
        @test ctrl.help_obs[] != Ext._HELP_HINT
        close(ctrl)
    end

    @testset "close disconnects listeners and removes the selection box" begin
        fig, ax, h, m1, m2 = _fixture()
        scene = ax.scene
        pick_m1 = ax2 -> (h.handles[1].plots[1], 0)
        n0 = length(ax.scene.plots)
        ctrl = Ext.kinematic_controls!(ax, h; throttle = false, pick = pick_m1)
        @test length(ax.scene.plots) == n0 + 4 # selection box, axes, labels and controls overlay

        close(ctrl)
        @test isempty(ctrl.listeners)
        @test length(ax.scene.plots) == n0

        # listeners are gone: this must not error and must not select anything
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        @test ctrl.selected[] === nothing
    end
end

end # module
