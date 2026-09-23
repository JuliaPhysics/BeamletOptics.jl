module TestLiveView

using BeamletOptics
using Makie
using Test

const BMO = BeamletOptics

@testset "Live view" begin
    Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)
    @test !isnothing(Ext)

    # Beam along +y, mirror at 45° reflects it along +x onto the detector
    function _fixture()
        m = RoundPlanoMirror(25e-3, 5e-3)
        zrotate3d!(m, deg2rad(45))
        translate3d!(m, [0, 0.1, 0])
        pd = Detector(5e-3)
        zrotate3d!(pd, -π / 2)
        translate3d!(pd, [0.1, 0.1, 0])
        return m, pd
    end

    # Picking requires a screen and is replaced by a custom `pick`, see the tests below
    function _select!(gui)
        scene = gui.ax.scene
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release)
        return nothing
    end

    _key!(gui, key) = (events(gui.ax.scene).keyboardbutton[] = Makie.KeyEvent(key, Keyboard.press))

    _mean_x(pts) = sum(p -> p[1], pts) / length(pts)

    _gauss() = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 1e-6, 0.5e-3)

    @testset "construction and panels" begin
        m, pd = _fixture()
        sys = System([m, pd])
        gauss = _gauss()
        gui = live_view(sys, gauss)
        @test gui isa Ext.LiveView
        @test length(gui.panels) == 1
        @test occursin("Detector 1: P =", gui.panels[1].ax.title[])
        @test occursin("mW", gui.panels[1].ax.title[])
        @test gui.panels[1].heat_plot.visible[]
        @test !gui.panels[1].scatter_plot.visible[]
        @test sprint(show, gui) == "LiveView(1 systems, 1 detector panels)"
        @test isnothing(gui.sliders)
        close(gui)

        m, pd = _fixture()
        cs = CollimatedSource([0.0, 0, 0], [0.0, 1, 0], 2e-3, 1e-6; num_rings = 2, num_rays = 40)
        gui = live_view(System([m, pd]) => cs)
        @test occursin("40 rays", gui.panels[1].ax.title[])
        @test gui.panels[1].scatter_plot.visible[]
        @test length(gui.panels[1].xy[]) == 40
        close(gui)

        # explicit modes and kwargs, a solved beam can not be reused for new objects
        m, pd = _fixture()
        gui = live_view(System([m, pd]), _gauss();
            detectors = [pd => (:intensity, (; n = 20, x_min = -2.5e-3, x_max = 2.5e-3,
                z_min = -2.5e-3, z_max = 2.5e-3))])
        @test size(gui.panels[1].heat_I[]) == (20, 20)
        @test gui.panels[1].heat_x[] ≈ collect(LinRange(-2.5f0, 2.5f0, 20))
        close(gui)
        gui = live_view(System([m, pd]), _gauss(); detectors = [pd => :spot])
        @test gui.panels[1].scatter_plot.visible[]
        close(gui)

        # no panels
        gui = live_view(System([m, pd]), _gauss(); detectors = [])
        @test isempty(gui.panels)
        close(gui)

        @test_throws ArgumentError live_view()
        @test_throws ArgumentError live_view(System([m, pd]), gauss; detectors = [pd => :fancy])
        @test_throws ArgumentError live_view(System([m, pd]), gauss; detectors = :none)
    end

    @testset "shared detector is emptied once per solve" begin
        m, pd = _fixture()
        sys1 = System([m, pd])
        sys2 = System([pd])
        b1 = Beam([0.0, 0, 0], [0.0, 1, 0])
        b2 = Beam([0.05, 0.1, 0.001], [1.0, 0, 0])
        gui = live_view(sys1 => b1, sys2 => b2; throttle = false, mode = :rotate, fine_angle = 1e-3)
        @test length(gui.panels) == 1 # deduplicated
        @test sprint(show, gui) == "LiveView(2 systems, 1 detector panels)"
        # both systems are handled by one controller
        @test length(gui.controls.h.handles) == 3
        @test length(BMO.hits(pd)) == 2
        gui.controls.selected[] = m
        _key!(gui, Keyboard.left)
        _key!(gui, Keyboard.left)
        n_gui = length(BMO.hits(pd))
        # reference: fresh solves of both pairs
        empty!(pd)
        solve_system!(sys1, b1)
        n1 = length(BMO.hits(pd))
        empty!(pd)
        solve_system!(sys2, b2)
        n2 = length(BMO.hits(pd))
        @test n1 == 1 && n2 == 1
        @test n_gui == n1 + n2
        close(gui)
    end

    @testset "key step moves the spot, no hits without error" begin
        m, pd = _fixture()
        sys = System([m, pd])
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        n_calls = Ref(0)
        gui_ref = Ref{Any}(nothing)
        gui = live_view(sys, beam; throttle = false, mode = :rotate, fine_angle = 1e-2,
            on_change = (g, obj) -> (n_calls[] += 1), pick = ax -> (gui_ref[].controls.h.handles[1].plots[1], 0))
        gui_ref[] = gui
        @test gui.controls.h.handles[1].obj === m
        n0 = n_calls[]
        _select!(gui)
        @test gui.controls.selected[] === m
        x0 = _mean_x(gui.panels[1].xy[])
        _key!(gui, Keyboard.left)
        @test n_calls[] == n0 + 1
        @test length(BMO.hits(pd)) == 1
        x1 = _mean_x(gui.panels[1].xy[])
        @test abs(x1 - x0) > 1 # [mm], 20 mrad deflection over 100 mm
        @test startswith(gui.status.text[], "$(nameof(typeof(m))) at (")

        # move the mirror out of the beam
        _key!(gui, Keyboard.m)
        gui.controls.fine_step = 0.1
        @test_logs _key!(gui, Keyboard.page_up)
        @test isnothing(BMO.hits(pd))
        @test gui.panels[1].ax.title[] == "Detector 1: no hits"
        @test !gui.panels[1].scatter_plot.visible[]
        @test !gui.panels[1].heat_plot.visible[]
        @test isempty(gui.panels[1].xy[])

        # close removes all listeners
        n1 = n_calls[]
        close(gui)
        @test isempty(gui.controls.listeners)
        _key!(gui, Keyboard.page_down)
        @test n_calls[] == n1
    end

    @testset "sliders" begin
        m, pd = _fixture()
        sys = System([m, pd])
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        called = Float64[]
        n_calls = Ref(0)
        callback = v -> (push!(called, v); translate_to3d!(pd, [0.1, 0.1, v * 1e-3]))
        gui = live_view(sys, beam; sliders = ["detector z [mm]" => (0:0.1:2, callback)],
            on_change = (g, obj) -> (obj === nothing && (n_calls[] += 1)))
        @test isempty(called) # not called at construction
        @test length(gui.sliders.sliders) == 1
        n0 = n_calls[]
        z0 = gui.panels[1].xy[][1][2]
        gui.sliders.sliders[1].value[] = 1.0
        @test isempty(called) # throttled until the next tick
        notify(events(gui.ax.scene).tick)
        @test called == [1.0]
        @test n_calls[] == n0 + 1
        z1 = gui.panels[1].xy[][1][2]
        @test abs(abs(z1 - z0) - 1) < 1e-3 # [mm]
        # no update without changes
        notify(events(gui.ax.scene).tick)
        @test n_calls[] == n0 + 1

        # errors of the callback are logged once
        close(gui)
        m, pd = _fixture()
        gui = live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); sliders = ["x" => (0:0.1:1, v -> error("slider"), 0.5)])
        @test gui.sliders.sliders[1].value[] == 0.5
        gui.sliders.sliders[1].value[] = 0.1
        @test_logs (:error, r"slider callback") notify(events(gui.ax.scene).tick)
        gui.sliders.sliders[1].value[] = 0.2
        @test_logs notify(events(gui.ax.scene).tick)
        close(gui)
        gui.sliders.sliders[1].value[] = 0.3
        @test_logs notify(events(gui.ax.scene).tick)
    end

    @testset "occluding mesh does not block ray picking" begin
        # Mirror at the origin, which the default camera looks at
        m = RoundPlanoMirror(0.025, 0.005)
        sys = System([m])
        beam = Beam([1.0, 1.0, 1.0], [1.0, 0, 0]) # unrelated to the mirror, live_view needs a beam
        gui = live_view(sys, beam; detectors = [])

        # A MeshDummy-like housing in front of the mirror along the camera ray, not part of the
        # system: added directly to the axis, as described for `render!(gui.ax, ...)`
        cube = BMO.CubeMesh(1)
        translate3d!(cube, -[0.5, 0.5, 0.5])
        translate3d!(cube, [1.5, 1.5, 1.5]) # between the default camera eye [3,3,3] and the mirror
        occluder = NonInteractableObject(cube)
        render!(gui.ax, occluder)

        scene = gui.ax.scene
        vp = scene.viewport[]
        cx, cy = vp.origin[1] + vp.widths[1] / 2, vp.origin[2] + vp.widths[2] / 2
        events(scene).mouseposition[] = (cx, cy)
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release)
        @test gui.controls.selected[] === m # the occluding housing does not block the selection
        close(gui)
    end

    @testset "panel power matches optical_power" begin
        m, pd = _fixture()
        # small area and coarse grid, such that the edges contribute to the integral
        area = (; n = 5, x_min = -0.3e-3, x_max = 0.3e-3, z_min = -0.3e-3, z_max = 0.3e-3)
        gui = live_view(System([m, pd]), _gauss(); detectors = [pd => (:intensity, area)])
        P = optical_power(pd; area...)
        @test gui.panels[1].ax.title[] == "Detector 1: P = $(Ext._fmt3(1e3 * P)) mW"
        close(gui)
    end

    @testset "failed solve keeps the beams outdated" begin
        m, pd = _fixture()
        gui = live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); auto_trace = false,
            throttle = false)
        @test !gui.stale
        # solve_system! fails for this beam
        gui.pairs[1] = gui.pairs[1].first => nothing
        @test_logs (:error, r"solving the systems") Ext._trace!(gui)
        @test gui.stale
        @test occursin("failed", gui.status.text[])
        close(gui)
    end

    @testset "manual trace" begin
        _spot(pd) = sort(collect(BMO.spot_diagram(pd)); by = p -> (p[1], p[2]))

        # beam points and spot of a fresh solve of the system
        function _fresh(sys, pd)
            b = Beam([0.0, 0, 0], [0.0, 1, 0])
            empty!(pd)
            solve_system!(sys, b)
            pts = Point3f[]
            Ext._collect_segments!(pts, b; flen = 1.0, render_every = 5)
            return pts, length(BMO.hits(pd)), _spot(pd)
        end

        @testset "key steps, t, trace button" begin
            m, pd = _fixture()
            sys = System([m, pd])
            beam = Beam([0.0, 0, 0], [0.0, 1, 0])
            n_calls = Ref(0)
            gui = live_view(sys, beam; auto_trace = false, throttle = false, mode = :rotate,
                fine_angle = 1e-2, on_change = (g, obj) -> (n_calls[] += 1))
            @test !gui.auto_trace[]
            @test !gui.auto_trace_toggle.active[]
            @test gui.trace_button isa Makie.Button
            # the initial solve runs anyway
            @test n_calls[] == 1
            @test length(BMO.hits(pd)) == 1
            @test !gui.stale
            bh = gui.beam_handles[1]

            gui.controls.selected[] = m
            R0 = Matrix{Float64}(BMO.orientation(m))
            n_hits = length(BMO.hits(pd))
            pts0 = copy(bh.points[])
            xy0 = copy(gui.panels[1].xy[])
            _key!(gui, Keyboard.left)
            # the object moves, but nothing is solved
            @test Matrix{Float64}(BMO.orientation(m)) ≈ BMO.rotate3d([0, 0, 1], 1e-2) * R0
            @test gui.system_handles[1].handles[1].R ≈ Matrix{Float64}(BMO.orientation(m))
            @test n_calls[] == 1
            @test length(BMO.hits(pd)) == n_hits
            @test bh.points[] == pts0
            @test gui.panels[1].xy[] == xy0
            @test gui.stale
            @test occursin("outdated, press t to trace", gui.status.text[])
            @test startswith(gui.status.text[], "$(nameof(typeof(m))) at (")
            @test bh.plot.alpha[] ≈ 0.3

            # t solves the system
            _key!(gui, Keyboard.t)
            @test n_calls[] == 2
            @test !gui.stale
            @test bh.plot.alpha[] ≈ 1.0
            pts_gui, n_gui, xy_gui = copy(bh.points[]), length(BMO.hits(pd)), copy(gui.panels[1].xy[])
            spot_gui = _spot(pd)
            @test xy_gui != xy0
            pts, n, spot = _fresh(sys, pd)
            @test n_gui == n
            @test pts_gui ≈ pts
            @test spot_gui ≈ spot

            # the trace button solves the system as well, rotate back to keep the spot on the detector
            _key!(gui, Keyboard.right)
            @test gui.stale
            @test n_calls[] == 2
            notify(gui.trace_button.clicks)
            @test n_calls[] == 3
            @test !gui.stale
            pts_gui, xy_gui2 = copy(bh.points[]), copy(gui.panels[1].xy[])
            @test xy_gui2 != xy_gui
            pts, n, spot = _fresh(sys, pd)
            @test n == 1
            @test pts_gui ≈ pts

            # switching auto trace on without changes does not solve
            gui.auto_trace_toggle.active[] = true
            @test n_calls[] == 3
            gui.auto_trace_toggle.active[] = false

            close(gui)
            _key!(gui, Keyboard.t)
            @test n_calls[] == 3
        end

        @testset "sliders and auto trace toggle" begin
            m, pd = _fixture()
            sys = System([m, pd])
            beam = Beam([0.0, 0, 0], [0.0, 1, 0])
            n_calls = Ref(0)
            called = Float64[]
            callback = v -> (push!(called, v); translate_to3d!(pd, [0.1, 0.1, v * 1e-3]))
            gui = live_view(sys, beam; auto_trace = false, throttle = false,
                sliders = ["detector z [mm]" => (0:0.1:2, callback)],
                on_change = (g, obj) -> (n_calls[] += 1))
            @test n_calls[] == 1
            xy0 = copy(gui.panels[1].xy[])
            gui.sliders.sliders[1].value[] = 1.0
            notify(events(gui.ax.scene).tick)
            # the callback and update_render! run, but no solve
            @test called == [1.0]
            @test gui.system_handles[1].handles[2].P ≈ BMO.position(pd)
            @test n_calls[] == 1
            @test gui.panels[1].xy[] == xy0
            @test gui.stale
            @test gui.status.text[] == "outdated, press t to trace"

            # switching auto trace on solves once, since the state is outdated
            gui.auto_trace_toggle.active[] = true
            @test gui.auto_trace[]
            @test n_calls[] == 2
            @test !gui.stale
            @test gui.panels[1].xy[] != xy0

            # auto tracing resumes
            gui.controls.selected[] = m
            _key!(gui, Keyboard.up)
            @test n_calls[] == 3
            @test !gui.stale
            gui.sliders.sliders[1].value[] = 0.5
            notify(events(gui.ax.scene).tick)
            @test n_calls[] == 4
            @test !gui.stale
            close(gui)
        end

        @testset "auto trace by default, t traces as well" begin
            m, pd = _fixture()
            n_calls = Ref(0)
            gui = live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); throttle = false,
                on_change = (g, obj) -> (n_calls[] += 1))
            @test gui.auto_trace[]
            @test gui.auto_trace_toggle.active[]
            gui.controls.selected[] = m
            _key!(gui, Keyboard.up)
            @test n_calls[] == 2
            @test !gui.stale
            _key!(gui, Keyboard.t)
            @test n_calls[] == 3
            close(gui)
        end
    end

    @testset "failing on_change is logged once" begin
        m, pd = _fixture()
        sys = System([m, pd])
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        gui_ref = Ref{Any}(nothing)
        @test_logs (:error, r"on_change") begin
            gui_ref[] = live_view(sys, beam; throttle = false, on_change = (g, obj) -> error("boom"),
                pick = ax -> (gui_ref[].controls.h.handles[1].plots[1], 0))
            gui = gui_ref[]
            _select!(gui)
            _key!(gui, Keyboard.up)
        end
        @test gui_ref[].last_error !== nothing
        close(gui_ref[])
    end
end

end # module
