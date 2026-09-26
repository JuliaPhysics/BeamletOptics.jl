module TestLiveView

using BeamletOptics
using Makie
using LinearAlgebra: normalize, dot, norm
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

    # The tests expect each change to be solved immediately. Adaptive tracing depends on the
    # measured solve time, which exceeds the default `trace_budget` on slow runners (e.g. CI with
    # coverage), hence it is disabled unless a test sets `trace_budget` itself.
    _live_view(args...; kwargs...) = live_view(args...; merge((; trace_budget = Inf), kwargs)...)

    @testset "construction and panels" begin
        m, pd = _fixture()
        sys = System([m, pd])
        gauss = _gauss()
        gui = _live_view(sys, gauss)
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
        gui = _live_view(System([m, pd]) => cs)
        @test occursin("40 rays", gui.panels[1].ax.title[])
        @test gui.panels[1].scatter_plot.visible[]
        @test length(gui.panels[1].xy[]) == 40
        close(gui)

        # explicit modes and kwargs, a solved beam can not be reused for new objects
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), _gauss();
            detectors = [pd => (:intensity, (; n = 20, x_min = -2.5e-3, x_max = 2.5e-3,
                z_min = -2.5e-3, z_max = 2.5e-3))])
        @test size(gui.panels[1].heat_I[]) == (20, 20)
        @test gui.panels[1].heat_x[] ≈ collect(LinRange(-2.5f0, 2.5f0, 20))
        close(gui)
        gui = _live_view(System([m, pd]), _gauss(); detectors = [pd => :spot])
        @test gui.panels[1].scatter_plot.visible[]
        close(gui)

        # no panels
        gui = _live_view(System([m, pd]), _gauss(); detectors = [])
        @test isempty(gui.panels)
        close(gui)

        @test_throws ArgumentError _live_view()
        @test_throws ArgumentError _live_view(System([m, pd]), gauss; detectors = [pd => :fancy])
        @test_throws ArgumentError _live_view(System([m, pd]), gauss; detectors = :none)
    end

    @testset "shared detector is emptied once per solve" begin
        m, pd = _fixture()
        sys1 = System([m, pd])
        sys2 = System([pd])
        b1 = Beam([0.0, 0, 0], [0.0, 1, 0])
        b2 = Beam([0.05, 0.1, 0.001], [1.0, 0, 0])
        gui = _live_view(sys1 => b1, sys2 => b2; throttle = false, mode = :rotate, fine_angle = 1e-3)
        @test length(gui.panels) == 1 # deduplicated
        @test sprint(show, gui) == "LiveView(2 systems, 1 detector panels)"
        # both systems and the markers of both sources are handled by one controller
        @test length(gui.controls.h.handles) == 5
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
        gui = _live_view(sys, beam; throttle = false, mode = :rotate, fine_angle = 1e-2,
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
        gui = _live_view(sys, beam; sliders = ["detector z [mm]" => (0:0.1:2, callback)],
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
        gui = _live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); sliders = ["x" => (0:0.1:1, v -> error("slider"), 0.5)])
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
        # unrelated to the mirror, live_view needs a beam; its marker is not on the camera ray
        beam = Beam([1.0, -1.0, 0.0], [1.0, 0, 0])
        gui = _live_view(sys, beam; detectors = [])
        # Makie's default camera instead of the initial Front-Right-Top view of `live_view`
        set_view(gui.ax, [3.0, 3, 3], [0.0, 0, 0], [0.0, 0, 1])

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

    @testset "ray picking with the orthographic projection" begin
        # Mirror away from the origin, the camera looks at it, as after zooming to a component
        m = RoundPlanoMirror(0.025, 0.005)
        translate3d!(m, [0.2, 0.3, 0.1])
        gui = _live_view(System([m]), Beam([1.0, -1.0, 0.0], [1.0, 0, 0]); detectors = [],
            orthographic = true)
        c = collect(BMO.position(m))
        set_view(gui.ax, c .+ 0.5 .* [1.0, -1, 1], c, [0.0, 0, 1])
        scene = gui.ax.scene
        cam = Makie.cameracontrols(scene)
        @test cam.settings.projectiontype[] == Makie.Orthographic
        vp = scene.viewport[]
        events(scene).mouseposition[] = (vp.origin[1] + vp.widths[1] / 2, vp.origin[2] + vp.widths[2] / 2)
        # the ray through the center starts at the eye and points to lookat
        origin, dir = Ext._cursor_ray(scene)
        @test isapprox(origin, collect(cam.eyeposition[]); atol = 1e-6)
        @test isapprox(dir, normalize(collect(cam.lookat[] - cam.eyeposition[])); atol = 1e-6)
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release)
        @test gui.controls.selected[] === m
        close(gui)
    end

    @testset "movable sources" begin
        _marker(gui, src) = gui.controls.h.handles[findfirst(oh -> oh.obj === src, gui.controls.h.handles)]

        # Beam: select the marker, move the source along its direction
        m, pd = _fixture()
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        gui_ref = Ref{Any}(nothing)
        gui = _live_view(System([m, pd]), beam; throttle = false, fine_step = 1e-3,
            pick = ax -> (_marker(gui_ref[], beam).plots[1], 0))
        gui_ref[] = gui
        marker = _marker(gui, beam)
        @test beam in gui.controls.movable
        _select!(gui)
        @test gui.controls.selected[] === beam
        _key!(gui, Keyboard.up)
        @test collect(BMO.position(beam)) ≈ [0, 1e-3, 0]
        @test marker.P ≈ [0, 1e-3, 0]
        # solved again from the new start point
        @test length(BMO.hits(pd)) == 1
        @test gui.beam_handles[1].points[][1] ≈ Point3f(0, 1e-3, 0)
        @test startswith(gui.status.text[], "Beam at (")
        # rotate the source, the spot moves on the detector
        x0 = _mean_x(gui.panels[1].xy[])
        _key!(gui, Keyboard.m)
        _key!(gui, Keyboard.page_up)
        _key!(gui, Keyboard.m)
        _key!(gui, Keyboard.backspace)
        @test collect(BMO.position(beam)) ≈ [0, 0, 0]
        @test collect(BMO.direction(beam)) ≈ [0, 1, 0]
        close(gui)

        # beam group: the whole source is moved, all beams are solved without preview tracing
        m, pd = _fixture()
        cs = CollimatedSource([0.0, 0, 0], [0.0, 1, 0], 2e-3, 1e-6; num_rings = 2, num_rays = 40)
        gui_ref = Ref{Any}(nothing)
        gui = _live_view(System([m, pd]) => cs; throttle = false, fine_step = 1e-3, preview = false,
            pick = ax -> (_marker(gui_ref[], cs).plots[1], 0))
        gui_ref[] = gui
        _select!(gui)
        @test gui.controls.selected[] === cs
        z0 = sum(p -> p[2], gui.panels[1].xy[]) / 40
        _key!(gui, Keyboard.page_up) # along the vertical axis
        @test collect(BMO.position(cs)) ≈ [0, 0, 1e-3]
        @test length(BMO.hits(pd)) == 40
        @test sum(p -> p[2], gui.panels[1].xy[]) / 40 ≈ z0 + 1 atol = 1e-3 # [mm]
        close(gui)

        # hide and show the markers via the toggle and the key `s`, a selected source is deselected
        m, pd = _fixture()
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        gui_ref = Ref{Any}(nothing)
        gui = _live_view(System([m, pd]), beam; throttle = false,
            pick = ax -> (_marker(gui_ref[], beam).plots[1], 0))
        gui_ref[] = gui
        marker = _marker(gui, beam)
        @test gui.sources_toggle.active[]
        @test all(p -> p.visible[], marker.plots)
        _select!(gui)
        @test gui.controls.selected[] === beam
        _key!(gui, Keyboard.s)
        @test !gui.sources_toggle.active[]
        @test !any(p -> p.visible[], marker.plots)
        @test isnothing(gui.controls.selected[])
        @test startswith(gui.status.text[], "sources hidden")
        # "show all" keeps the markers hidden, the mirror stays visible
        Ext._show_all!(gui)
        @test !any(p -> p.visible[], marker.plots)
        gui.sources_toggle.active[] = true
        @test all(p -> p.visible[], marker.plots)
        close(gui)
        # initially hidden
        m, pd = _fixture()
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        gui = _live_view(System([m, pd]), beam; show_sources = false)
        @test !any(p -> p.visible[], _marker(gui, beam).plots)
        close(gui)

        # no markers
        m, pd = _fixture()
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        gui = _live_view(System([m, pd]), beam; movable_sources = false)
        @test !any(oh -> oh.obj === beam, gui.controls.h.handles)
        @test !(beam in gui.controls.movable)
        close(gui)
    end

    @testset "labels, pose and step textbox" begin
        m, pd = _fixture()
        gui_ref = Ref{Any}(nothing)
        gui = _live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); throttle = false,
            fine_step = 1e-3, labels = Dict(m => "Mirror 1", pd => "PD"),
            pick = ax -> (gui_ref[].controls.h.handles[1].plots[1], 0))
        gui_ref[] = gui
        @test startswith(gui.panels[1].ax.title[], "PD: ")
        _select!(gui)
        _key!(gui, Keyboard.up)
        @test startswith(gui.status.text[], "Mirror 1 at (")
        @test occursin("moved by 1 mm, rotated by 0 µrad", gui.status.text[])

        gui.step_box.stored_string[] = "250 nm"
        @test gui.controls.fine_step ≈ 250e-9
        @test gui.controls.mode[] == :move
        gui.step_box.stored_string[] = "50 µrad"
        @test gui.controls.fine_angle ≈ 50e-6
        @test gui.controls.mode[] == :rotate
        gui.step_box.stored_string[] = "fast"
        @test occursin("invalid step", gui.status.text[])
        @test gui.controls.fine_angle ≈ 50e-6

        # typing into the textbox does not trigger the controls
        R0 = Matrix{Float64}(BMO.orientation(m))
        gui.step_box.focused[] = true
        _key!(gui, Keyboard.left)
        _key!(gui, Keyboard.m)
        @test Matrix{Float64}(BMO.orientation(m)) == R0
        @test gui.controls.mode[] == :rotate
        gui.step_box.focused[] = false
        _key!(gui, Keyboard.left)
        @test Matrix{Float64}(BMO.orientation(m)) ≈ BMO.rotate3d([0, 0, 1], 50e-6) * R0
        close(gui)

        @test Ext._length_string(0.9999999e-3) == "1 mm"
        @test Ext._length_string(2.5e-7) == "250 nm"
        @test Ext._length_string(12.0) == "12000 mm"
        @test Ext._angle_string(5e-5) == "50 µrad"
        @test Ext._angle_string(0.1) == "100 mrad"
        @test Ext._angle_string(deg2rad(90)) == "90 °"
        @test Ext._parse_step("250 nm")[1] == :move
        @test Ext._parse_step("250 nm")[2] ≈ 250e-9
        @test Ext._parse_step("0.5um")[2] ≈ 0.5e-6
        @test Ext._parse_step("1e-3 m")[2] ≈ 1e-3
        @test Ext._parse_step("2 deg")[2] ≈ deg2rad(2)
        @test Ext._parse_step("3 mrad")[1] == :rotate
        @test isnothing(Ext._parse_step("10"))
        @test isnothing(Ext._parse_step("-1 nm"))
        @test isnothing(Ext._parse_step("1 inch"))
    end

    @testset "adaptive tracing" begin
        # slow systems: the solve is deferred until the movement pauses
        m, pd = _fixture()
        gui_ref = Ref{Any}(nothing)
        gui = _live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); throttle = false,
            mode = :rotate, fine_angle = 1e-2, trace_budget = 0.0, idle_delay = 0.1,
            pick = ax -> (gui_ref[].controls.h.handles[1].plots[1], 0))
        gui_ref[] = gui
        gui.solve_time = 1.0 # pretend that the solve is slow, independent of the machine
        _select!(gui)
        pts0 = copy(gui.beam_handles[1].points[])
        _key!(gui, Keyboard.left)
        @test gui.pending
        @test gui.stale
        @test gui.beam_handles[1].points[] == pts0
        @test occursin("tracing when the movement pauses", gui.status.text[])
        # still moving
        notify(events(gui.ax.scene).tick)
        @test gui.pending
        gui.last_change -= 1
        notify(events(gui.ax.scene).tick)
        @test !gui.pending
        @test !gui.stale
        @test gui.beam_handles[1].points[] != pts0
        close(gui)

        # slow panels: a coarse preview while moving, refined once the movement pauses
        m, pd = _fixture()
        gui_ref = Ref{Any}(nothing)
        gui = _live_view(System([m, pd]), _gauss(); throttle = false, mode = :rotate,
            fine_angle = 1e-4, idle_delay = 0.1, detectors = [pd => (:intensity, (; n = 40))],
            trace_budget = 0.5, pick = ax -> (gui_ref[].controls.h.handles[1].plots[1], 0))
        gui_ref[] = gui
        @test size(gui.panels[1].heat_I[]) == (40, 40)
        # pretend that the solve is fast and the panels are slow, independent of the machine
        gui.solve_time = 0.0
        gui.panel_time = 1.0
        _select!(gui)
        _key!(gui, Keyboard.left)
        @test gui.coarse
        @test size(gui.panels[1].heat_I[]) == (16, 16)
        @test endswith(gui.panels[1].ax.title[], "(preview)")
        gui.last_change -= 1
        notify(events(gui.ax.scene).tick)
        @test !gui.coarse
        @test size(gui.panels[1].heat_I[]) == (40, 40)
        @test !endswith(gui.panels[1].ax.title[], "(preview)")
        close(gui)
    end

    @testset "panel power matches optical_power" begin
        m, pd = _fixture()
        # small area and coarse grid, such that the edges contribute to the integral
        area = (; n = 5, x_min = -0.3e-3, x_max = 0.3e-3, z_min = -0.3e-3, z_max = 0.3e-3)
        gui = _live_view(System([m, pd]), _gauss(); detectors = [pd => (:intensity, area)])
        P = optical_power(pd; area...)
        @test gui.panels[1].ax.title[] == "Detector 1: P = $(Ext._fmt3(1e3 * P)) mW"
        close(gui)
    end

    @testset "failed solve keeps the beams outdated" begin
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); auto_trace = false,
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
            gui = _live_view(sys, beam; auto_trace = false, throttle = false, mode = :rotate,
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
            gui = _live_view(sys, beam; auto_trace = false, throttle = false,
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
            gui = _live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); throttle = false,
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

    @testset "clip planes" begin
        _handle(gui, obj) = gui.controls.h.handles[findfirst(oh -> oh.obj === obj, gui.controls.h.handles)]
        # all plots of the system objects, including the nested plots of recipes
        _nested(p) = AbstractPlot[p; reduce(vcat, _nested.(p.plots); init = AbstractPlot[])]
        _optics_plots(gui) = reduce(vcat, (_nested(p) for h in gui.system_handles for oh in h.handles
                                           for p in oh.plots); init = AbstractPlot[])
        # markers of the sources and clip planes
        _marker_plots(gui) = reduce(vcat, (_nested(p) for oh in gui.controls.h.handles
                                           if !(oh.obj isa BMO.AbstractObject) for p in oh.plots); init = AbstractPlot[])
        _beam_plots(gui) = reduce(vcat, (_nested(p) for h in gui.beam_handles for p in Ext._beam_plots(h)))
        _control_plots(gui) = reduce(vcat, (_nested(p) for p in gui.controls.plots[1:4]))
        _planes(gui) = only(unique(p.clip_planes[] for p in _optics_plots(gui)))
        _shift_key!(gui, key) = (push!(events(gui.ax.scene).keyboardstate, Keyboard.left_shift);
                                 _key!(gui, key);
                                 delete!(events(gui.ax.scene).keyboardstate, Keyboard.left_shift))
        _beam() = Beam([0.0, 0, 0], [0.0, 1, 0])
        P1 = [Plane3f(Point3f(0, 0.1, 0), Vec3f(0, 1, 0))]

        @testset "kwarg, markers and beams are not clipped" begin
            m, pd = _fixture()
            n_calls = Ref(0)
            pick_plot = Ref{Any}(nothing)
            gui = _live_view(System([m, pd]), _beam(); throttle = false, fine_step = 1e-3,
                clip_planes = [[0, 0.1, 0] => [0, 1, 0]], on_change = (g, obj) -> (n_calls[] += 1),
                pick = ax -> (pick_plot[], 0))
            @test length(gui.clip_planes) == 1
            plane = gui.clip_planes[1]
            @test isnothing(gui.controls.selected[])
            @test gui.clipping
            @test plane in gui.controls.movable
            @test gui.ax.scene.theme.clip_planes[] == P1
            @test !isempty(_optics_plots(gui))
            @test all(p -> p.clip_planes[] == P1, _optics_plots(gui))
            for plots in (_marker_plots(gui), _beam_plots(gui), _control_plots(gui))
                @test !isempty(plots)
                @test all(p -> p.clip_planes[] == Plane3f[], plots)
            end

            # a plot added later gets the planes, also after a move
            added = Ext._capture_new_plots(gui.ax) do
                cube = BMO.CubeMesh(0.01)
                render!(gui.ax, NonInteractableObject(cube))
            end
            @test !isempty(added)
            @test all(p -> p.clip_planes[] == P1, reduce(vcat, _nested.(added)))

            # the outline does not select the plane, the handle does
            marker = _handle(gui, plane)
            pick_plot[] = marker.plots[1]
            @test marker.plots[1] isa Makie.Lines
            _select!(gui)
            @test isnothing(gui.controls.selected[])
            pick_plot[] = marker.plots[2]
            _select!(gui)
            @test gui.controls.selected[] === plane
            @test all(isfinite, reduce(vcat, collect.(gui.controls.box_obs[])))

            # a key step along the normal moves the plane, without solving the systems
            n0, spot0 = n_calls[], copy(BMO.spot_diagram(pd))
            stale0 = gui.stale
            _key!(gui, Keyboard.up)
            @test collect(position(plane)) ≈ [0, 0.101, 0]
            @test abs(_planes(gui)[1].distance - P1[1].distance - 1e-3) < 1e-6
            @test _planes(gui) == [Ext._plane3f(plane)]
            @test all(p -> p.clip_planes[] == _planes(gui), reduce(vcat, _nested.(added)))
            @test marker.P ≈ position(plane)
            @test n_calls[] == n0
            @test BMO.spot_diagram(pd) == spot0
            @test gui.stale == stale0
            @test startswith(gui.status.text[], "Clip plane at (")

            # shift+c flips the selected plane, c switches clipping off and on
            _shift_key!(gui, Keyboard.c)
            @test Ext._normal(plane) ≈ [0, -1, 0]
            @test _planes(gui)[1].normal ≈ Vec3f(0, -1, 0)
            @test collect(position(plane)) ≈ [0, 0.101, 0]
            @test gui.clipping
            flipped = _planes(gui)
            _key!(gui, Keyboard.c)
            @test !gui.clipping
            @test _planes(gui) == Plane3f[]
            @test all(p -> p.clip_planes[] == Plane3f[], reduce(vcat, _nested.(added)))
            _key!(gui, Keyboard.c)
            @test gui.clipping
            @test _planes(gui) == flipped
            @test n_calls[] == n0

            # reset of the plane
            _key!(gui, Keyboard.backspace)
            @test collect(position(plane)) ≈ [0, 0.1, 0]
            @test _planes(gui)[1].normal ≈ P1[1].normal
            @test _planes(gui)[1].distance ≈ P1[1].distance
            @test n_calls[] == n0

            # `v` still toggles the spectator mode with a clip plane selected
            @test gui.controls.selected[] === plane
            _key!(gui, Keyboard.v)
            @test gui.controls.spectator[]
            @test isnothing(gui.controls.selected[])
            _key!(gui, Keyboard.v)
            @test !gui.controls.spectator[]
            close(gui)

            # clip_beams
            m, pd = _fixture()
            gui = _live_view(System([m, pd]), _gauss(); clip_planes = [[0, 0.1, 0] => [0, 1, 0]],
                clip_beams = true, beam_kwargs = Dict())
            @test !isempty(_beam_plots(gui))
            @test all(p -> p.clip_planes[] == P1, _beam_plots(gui))
            @test all(p -> p.clip_planes[] == Plane3f[], _marker_plots(gui))
            @test gui.clip_beams_toggle.active[]
            # switched off and on again at runtime via the toggle
            gui.clip_beams_toggle.active[] = false
            @test !gui.clip_beams
            @test all(p -> p.clip_planes[] == Plane3f[], _beam_plots(gui))
            gui.clip_beams_toggle.active[] = true
            @test all(p -> p.clip_planes[] == P1, _beam_plots(gui))
            # clipping off: the beams follow as well
            _key!(gui, Keyboard.c)
            @test all(p -> p.clip_planes[] == Plane3f[], _beam_plots(gui))
            close(gui)

            # default: the beams are not clipped until the toggle is switched on
            m, pd = _fixture()
            gui = _live_view(System([m, pd]), _gauss(); clip_planes = [[0, 0.1, 0] => [0, 1, 0]],
                beam_kwargs = Dict())
            @test !gui.clip_beams_toggle.active[]
            @test all(p -> p.clip_planes[] == Plane3f[], _beam_plots(gui))
            gui.clip_beams_toggle.active[] = true
            @test all(p -> p.clip_planes[] == P1, _beam_plots(gui))
            close(gui)

            # invalid planes
            m, pd = _fixture()
            @test_throws ArgumentError _live_view(System([m, pd]), _beam();
                clip_planes = [[0, 0, i] => [0, 0, 1] for i in 1:9])
            @test_throws ArgumentError _live_view(System([m, pd]), _beam(); clip_planes = [[0, 0, 0] => [0, 0, 0]])
            @test_throws ArgumentError Ext.LiveClipPlane([0, 0, 0], [0, 0, 0], 1.0)
        end

        @testset "drag of a plane with its normal along the rotation axis" begin
            # the normal (local y-axis) is parallel to the rotation axis z, e.g. for a plane added
            # in the top view, hence the allowed axes of the drag are linearly dependent
            m, pd = _fixture()
            pick_plot = Ref{Any}(nothing)
            gui = _live_view(System([m, pd]), _beam(); throttle = false,
                clip_planes = [[0, 0.1, 0] => [0, 0, 1]], pick = ax -> (pick_plot[], 0))
            plane = gui.clip_planes[1]
            pick_plot[] = _handle(gui, plane).plots[2]
            scene = gui.ax.scene
            events(scene).mouseposition[] = (100.0, 100.0)
            _select!(gui)
            @test gui.controls.selected[] === plane
            P0 = collect(Float64.(BMO.position(plane)))
            events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
            events(scene).mouseposition[] = (140.0, 160.0)
            @test gui.controls.dragging
            P1 = collect(Float64.(BMO.position(plane)))
            @test all(isfinite, P1)
            @test P1 != P0
            events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release)
            close(gui)
        end

        @testset "fully clipped selection" begin
            m, pd = _fixture()
            # everything below y = 0.5 is clipped, i.e. all objects
            gui = _live_view(System([m, pd]), _beam(); throttle = false,
                clip_planes = [[0, 0.5, 0] => [0, 1, 0]])
            @test Makie.boundingbox(gui.system_handles[1].handles[1].plots[1]) == Makie.Rect3d()
            gui.controls.selected[] = m
            Ext._update_selection_box!(gui.controls)
            pts = gui.controls.box_obs[]
            @test !isempty(pts)
            @test all(p -> all(isfinite, p), pts)
            @test all(p -> all(isfinite, p), gui.controls.arrow_pos[])
            close(gui)
        end

        @testset "add and remove at runtime" begin
            m, pd = _fixture()
            n_calls = Ref(0)
            gui = _live_view(System([m, pd]), _beam(); throttle = false, fine_step = 1e-3,
                fine_angle = 1e-2, on_change = (g, obj) -> (n_calls[] += 1))
            @test isempty(gui.clip_planes)
            # no clip planes: nothing is written
            @test all(p -> p.clip_planes[] == Plane3f[], _optics_plots(gui))
            @test occursin("p: add clip plane", Ext._help_text(:move, 1e-9, 1e-6) * "\n" * gui.controls.help_extra)
            gui.controls.help_shown = true
            Ext._update_help!(gui.controls)
            @test occursin("shift+c: flip", gui.controls.help_obs[])
            n0 = n_calls[]

            # nothing selected: through the lookat point of the camera, along the view direction
            cam = Makie.cameracontrols(gui.ax.scene)
            cam.eyeposition[] = Vec3f(0.3, -0.2, 0.25)
            cam.lookat[] = Vec3f(0.05, 0.1, 0.0)
            lookat, eye = Vector{Float64}(cam.lookat[]), Vector{Float64}(cam.eyeposition[])
            view_dir = (lookat - eye) / sqrt(sum(abs2, lookat - eye))
            _key!(gui, Keyboard.p)
            @test length(gui.clip_planes) == 1
            plane = gui.clip_planes[1]
            @test gui.controls.selected[] === plane
            @test collect(position(plane)) ≈ lookat
            @test maximum(abs.(Ext._normal(plane) - view_dir)) < 1e-6
            @test all(p -> p.clip_planes[] == [Ext._plane3f(plane)], _optics_plots(gui))
            @test all(p -> p.clip_planes[] == Plane3f[], _marker_plots(gui))

            # moved and rotated with keys like a component
            _key!(gui, Keyboard.up)
            @test collect(position(plane)) ≈ lookat + 1e-3 * view_dir
            x_axis = plane.dir[:, 1]
            _key!(gui, Keyboard.m)
            _key!(gui, Keyboard.up)
            @test Ext._normal(plane) ≈ BMO.rotate3d(x_axis, 1e-2) * view_dir
            @test all(p -> p.clip_planes[] == [Ext._plane3f(plane)], _optics_plots(gui))
            _key!(gui, Keyboard.m)
            @test n_calls[] == n0

            # a selected component: through its position
            gui.controls.selected[] = m
            _key!(gui, Keyboard.p)
            @test length(gui.clip_planes) == 2
            plane2 = gui.clip_planes[2]
            @test gui.controls.selected[] === plane2
            @test collect(position(plane2)) ≈ collect(position(m))
            @test maximum(abs.(Ext._normal(plane2) - view_dir)) < 1e-6
            @test _planes(gui) == Ext._plane3f.([plane, plane2])

            # Delete with a component selected does nothing
            gui.controls.selected[] = m
            _key!(gui, Keyboard.delete)
            @test length(gui.clip_planes) == 2
            @test gui.controls.selected[] === m

            # Delete removes the marker and the plane
            marker_plots = copy(_handle(gui, plane2).plots)
            gui.controls.selected[] = plane2
            _key!(gui, Keyboard.delete)
            @test gui.clip_planes == [plane]
            @test isnothing(gui.controls.selected[])
            @test !(plane2 in gui.controls.movable)
            @test !haskey(gui.controls.init_poses, plane2)
            @test !any(oh -> oh.obj === plane2, gui.controls.h.handles)
            @test !any(p -> any(q -> q === p, gui.ax.scene.plots), marker_plots)
            @test _planes(gui) == [Ext._plane3f(plane)]
            # the undo history of the removed plane is dropped, undo still works
            @test all(e -> e.obj !== plane2, gui.controls.undo_stack)
            gui.controls.selected[] = plane
            _key!(gui, Keyboard.delete)
            @test isempty(gui.clip_planes)
            @test all(p -> p.clip_planes[] == Plane3f[], _optics_plots(gui))
            @test gui.ax.scene.theme.clip_planes[] == Plane3f[]
            @test n_calls[] == n0

            # at most 8 planes
            for _ in 1:8
                _key!(gui, Keyboard.p)
            end
            @test length(gui.clip_planes) == 8
            @test_logs _key!(gui, Keyboard.p)
            @test length(gui.clip_planes) == 8
            @test gui.status.text[] == "at most 8 clip planes"
            @test length(_planes(gui)) == 8
            @test n_calls[] == n0
            close(gui)
        end
    end

    @testset "orthographic toggle" begin
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]))
        settings = cameracontrols(gui.ax.scene).settings
        @test !gui.orthographic_toggle.active[]
        @test settings.projectiontype[] == Makie.Perspective
        gui.orthographic_toggle.active[] = true
        @test settings.projectiontype[] == Makie.Orthographic
        gui.orthographic_toggle.active[] = false
        @test settings.projectiontype[] == Makie.Perspective
        close(gui)

        m, pd = _fixture()
        gui = _live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); orthographic = true)
        @test gui.orthographic_toggle.active[]
        @test cameracontrols(gui.ax.scene).settings.projectiontype[] == Makie.Orthographic
        close(gui)
    end

    @testset "spot panel limits" begin
        # A single ray gives a zero-width spot diagram, which must not collapse the limits
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]))
        p = gui.panels[1]
        @test length(p.xy[]) == 1
        for lims in (p.ax.targetlimits[], p.ax.finallimits[])
            @test all(isfinite, lims.origin) && all(isfinite, lims.widths)
            @test all(lims.widths .>= 2e-3 * (1 - 1e-6))
        end
        close(gui)

        # Many spots keep the limits of `autolimits!`
        m, pd = _fixture()
        cs = CollimatedSource([0.0, 0, 0], [0.0, 1, 0], 2e-3, 1e-6; num_rings = 2, num_rays = 40)
        gui = _live_view(System([m, pd]) => cs)
        p = gui.panels[1]
        lims = p.ax.targetlimits[]
        autolimits!(p.ax)
        @test p.ax.targetlimits[] == lims
        close(gui)

        # Only the degenerate axis is padded, by half the extent of the other axis
        ax = Axis(Figure()[1, 1])
        Ext._pad_degenerate_limits!(ax, [Point2f(1, 0), Point2f(1, 2)])
        lims = ax.targetlimits[]
        @test lims.origin ≈ [0.0, 1 - 1.05] && lims.widths ≈ [2.0, 2 * 1.05]
    end

    @testset "initial view from the Front-Right-Top corner" begin
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); throttle = false)
        cam = cameracontrols(gui.ax.scene)
        dir = normalize(Vector{Float64}(cam.eyeposition[] .- cam.lookat[]))
        @test isapprox(dir, normalize([1.0, -1, 1]); atol = 1e-6)
        @test abs(dot(Vector{Float64}(cam.upvector[]), dir)) < 1e-6
        # A later view of the user is kept, also after tracing
        set_view(gui.ax, [0.0, 0, 1], [0.0, 0, 0], [0.0, 1, 0])
        _key!(gui, Keyboard.t)
        @test cam.eyeposition[] ≈ Vec3f(0, 0, 1)
        @test cam.upvector[] ≈ Vec3f(0, 1, 0)
        close(gui)
    end

    @testset "failing on_change is logged once" begin
        m, pd = _fixture()
        sys = System([m, pd])
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        gui_ref = Ref{Any}(nothing)
        @test_logs (:error, r"on_change") begin
            gui_ref[] = _live_view(sys, beam; throttle = false, on_change = (g, obj) -> error("boom"),
                pick = ax -> (gui_ref[].controls.h.handles[1].plots[1], 0))
            gui = gui_ref[]
            _select!(gui)
            _key!(gui, Keyboard.up)
        end
        @test gui_ref[].last_error !== nothing
        close(gui_ref[])
    end

    _ctrl_z!(gui) = (push!(events(gui.ax.scene).keyboardstate, Keyboard.left_control);
                     _key!(gui, Keyboard.z);
                     delete!(events(gui.ax.scene).keyboardstate, Keyboard.left_control))

    # Angle between the rotation matrices R1 and R2
    _angle(R1, R2) = Ext._rotation_axis_angle(R1 * R2')[2]

    @testset "export changes" begin
        m, pd = _fixture()
        m2 = RoundPlanoMirror(0.02, 0.004)
        translate3d!(m2, [0.2, 0, 0])
        g = ObjectGroup([RoundPlanoMirror(0.02, 0.004), RoundPlanoMirror(0.02, 0.004)])
        translate3d!(g.objects[2], [0, 0.05, 0])
        translate3d!(g, [0.3, 0, 0])
        lens = SphericalLens(0.05, -0.05, 5e-3, 25.4e-3)
        translate3d!(lens, [0, 0.05, 0])
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        cs = CollimatedSource([0.0, 0, 0], [0.0, 1, 0], 2e-3, 1e-6; num_rings = 2, num_rays = 40)
        # all objects of the menu and their fresh copies in the initial poses
        objs = Any[m, pd, m2, g, g.objects[1], g.objects[2], lens, beam, cs]
        fresh = deepcopy(objs)
        sys = System([m, pd, m2, g, lens])
        gui = _live_view(sys => beam, sys => cs; throttle = false, detectors = [],
            labels = Dict(m => "m1", lens => "the lens", g => "end", pd => "PD"))
        gui.export_clipboard = false
        entries = Ext._menu_entries(gui.controls)
        @test all(first.(entries) .=== objs)
        @test last.(entries) == [0, 0, 0, 0, 1, 1, 0, 0, 0]

        code = export_changes(gui; io = devnull)
        @test occursin("# no changes", code)
        @test !occursin("translate_to3d!", code)

        # mirror: selected via the menu, moved and rotated via keys
        gui.menu.i_selected[] = 1
        @test gui.controls.selected[] === m
        _key!(gui, Keyboard.up)
        _key!(gui, Keyboard.m)
        _key!(gui, Keyboard.left)
        _key!(gui, Keyboard.page_up)
        # group and one of its objects, lens, sources
        translate3d!(g, [0.01, 0.002, -0.001])
        zrotate3d!(g, 0.1)
        xrotate3d!(g, 0.02)
        xrotate3d!(g.objects[1], 2e-9)
        translate3d!(g.objects[1], [0, 0, 1e-3])
        translate3d!(lens, [1e-3, 0, 0])
        translate3d!(beam, [1e-3, 0, 0])
        rotate3d!(beam, [0.0, 0, 1], 0.01)
        translate3d!(cs, [0, 0, 2e-3])
        rotate3d!(cs, [1.0, 0, 0], 0.02)
        # clip planes are not exported
        plane = Ext._add_clip_plane!(gui, [0, 0.05, 0], [0, 1, 0])
        translate3d!(plane, [0, 0.01, 0])

        buf = IOBuffer()
        code = export_changes(gui; io = buf)
        @test String(take!(buf)) == code
        names = Ext._export_names(gui, objs)
        @test names[m] == "m1"
        @test names[lens] == "obj7" && names[g] == "obj4" && names[pd] == "PD" && names[m2] == "obj3"
        @test occursin("# m1 (Mirror)\nrotate3d!(m1, [", code)
        @test occursin("# the lens (Lens)\ntranslate_to3d!(obj7, [", code)
        @test occursin("# end (ObjectGroup)", code)
        @test occursin("\n# Beam\n", code)
        # unmoved objects, including the object of the moved group, and the clip plane
        for obj in (pd, m2, g.objects[2])
            @test !occursin("($(names[obj]),", code)
        end
        @test !occursin("Clip", code)
        @test count("translate_to3d!(", code) == 6

        # the code reproduces the poses on the fresh copies
        mod = Module()
        Core.eval(mod, :(using BeamletOptics))
        for (obj, f) in zip(objs, fresh)
            Core.eval(mod, :($(Symbol(names[obj])) = $f))
        end
        include_string(mod, code)
        for (obj, f) in zip(objs, fresh)
            P, R = Ext._pose(obj)
            Pf, Rf = Ext._pose(f)
            @test maximum(abs, P - Pf) < 1e-12
            @test _angle(R, Rf) < 1e-12
        end

        # the export button prints the same code
        out = mktemp() do path, io
            redirect_stdout(io) do
                notify(gui.export_button.clicks)
            end
            flush(io)
            read(path, String)
        end
        @test out == code
        @test gui.status.text[] == "exported 6 changes"
        close(gui)

        # accurate axis and angle, also for small angles and close to π
        for (axis, angle) in (([1.0, 2, 3], 1e-9), ([0.0, 0, 1], 0.3), ([1.0, -1, 0.5], π - 1e-9),
                ([0.3, 0.2, -1], 2.5), ([0.0, 1, 0], π))
            R = BMO.rotate3d(axis, angle)
            a, θ = Ext._rotation_axis_angle(R)
            @test θ ≈ angle rtol = 1e-12
            @test maximum(abs, BMO.rotate3d(a, θ) - R) < 1e-14
        end
        @test Ext._rotation_axis_angle([1.0 0 0; 0 1 0; 0 0 1])[2] == 0
    end

    @testset "pose inspector" begin
        m, pd = _fixture()
        lens = SphericalLens(0.05, -0.05, 5e-3, 25.4e-3)
        translate3d!(lens, [0.01, 0.05, 0.002])
        gui_ref = Ref{Any}(nothing)
        gui = _live_view(System([lens, m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); throttle = false,
            fine_angle = 1e-3, fine_step = 1e-3, pick = ax -> (gui_ref[].controls.h.handles[2].plots[1], 0))
        gui_ref[] = gui
        boxes = gui.pose_boxes
        texts() = [tb.displayed_string[] for tb in boxes]
        @test texts() == fill("", 6)

        # nothing selected: input is ignored
        boxes[1].stored_string[] = "5"
        @test collect(BMO.position(lens)) == [0.01, 0.05, 0.002]
        @test texts() == fill("", 6)

        # the boxes follow the selection
        gui.menu.i_selected[] = 1
        @test gui.controls.selected[] === lens
        @test texts() == ["10.0", "50.0", "2.0", "", "", ""]
        _select!(gui)
        @test gui.controls.selected[] === m
        @test texts() == ["0.0", "100.0", "0.0", "", "", ""]
        gui.menu.i_selected[] = 1

        # absolute position [mm]
        P0 = collect(BMO.position(lens))
        R0 = Matrix(BMO.orientation(lens))
        boxes[1].stored_string[] = "12.5"
        P = collect(BMO.position(lens))
        @test P[1] == 0.0125
        @test P[2:3] == P0[2:3]
        @test BMO.orientation(lens) == R0
        @test texts() == ["12.5", "50.0", "2.0", "", "", ""]
        @test occursin("at (12.5, 50.0, 2.0) mm", gui.status.text[])

        # a rotation of 1 mrad about the blue axis equals one key step
        ref = deepcopy(lens)
        _key!(gui, Keyboard.m)
        @test gui.controls.mode[] == :rotate
        @test Ext._key_step!(gui.controls, ref, Keyboard.left, 1)
        boxes[6].stored_string[] = "1"
        @test BMO.orientation(lens) == BMO.orientation(ref)
        @test collect(BMO.position(lens)) == P
        @test texts() == ["12.5", "50.0", "2.0", "", "", ""]
        # red and green axes, like the keys up and page up
        for (k, key) in ((4, Keyboard.up), (5, Keyboard.page_up))
            Ext._key_step!(gui.controls, ref, key, 1)
            boxes[k].stored_string[] = "1"
            @test BMO.orientation(lens) ≈ BMO.orientation(ref) atol = 1e-15
        end

        # undo
        for _ in 1:3
            _ctrl_z!(gui)
        end
        @test _angle(Matrix(BMO.orientation(lens)), R0) < 1e-12
        _ctrl_z!(gui)
        @test collect(BMO.position(lens)) ≈ P0 atol = 1e-15
        @test texts()[1] == "10.0"

        # invalid input changes nothing
        for s in ("abc", "1 mm", "NaN")
            boxes[2].stored_string[] = s
            @test startswith(gui.status.text[], "invalid input \"$s\"")
            @test collect(BMO.position(lens)) ≈ P0 atol = 1e-15
        end

        # keys are ignored while a box is focused, the focused box is not overwritten
        _key!(gui, Keyboard.m)
        @test gui.controls.mode[] == :move
        boxes[1].focused[] = true
        boxes[1].displayed_string[] = "3"
        _key!(gui, Keyboard.up)
        _key!(gui, Keyboard.m)
        @test gui.controls.mode[] == :move
        @test collect(BMO.position(lens)) ≈ P0 atol = 1e-15
        translate3d!(lens, [0, 0, 1e-3])
        Ext._request_update!(gui.controls)
        @test texts()[1:3] == ["3", "50.0", "3.0"]
        boxes[1].focused[] = false
        # the boxes follow key steps
        _key!(gui, Keyboard.up)
        @test texts()[1:3] == ["10.0", "51.0", "3.0"]

        # deselect
        _key!(gui, Keyboard.escape)
        @test isnothing(gui.controls.selected[])
        @test texts() == fill("", 6)
        close(gui)

        # constraints are respected
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); throttle = false,
            constraints = Dict(m => (; move = (:v,), rotate = ())))
        gui.menu.i_selected[] = 1
        P0, R0 = collect(BMO.position(m)), Matrix(BMO.orientation(m))
        gui.pose_boxes[1].stored_string[] = "5"
        gui.pose_boxes[6].stored_string[] = "5"
        @test collect(BMO.position(m)) == P0
        @test BMO.orientation(m) == R0
        gui.pose_boxes[3].stored_string[] = "5"
        @test collect(BMO.position(m)) ≈ [0, 0.1, 0.005]
        close(gui)
    end

    @testset "component menu, hide and show" begin
        m, pd = _fixture()
        g = ObjectGroup([RoundPlanoMirror(0.02, 0.004), RoundPlanoMirror(0.02, 0.004)])
        translate3d!(g.objects[2], [0, 0.05, 0])
        translate3d!(g, [0.3, 0, 0])
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        gui_ref = Ref{Any}(nothing)
        gui = _live_view(System([m, pd, g]), beam; throttle = false,
            labels = Dict(m => "M1", g.objects[1] => "G1"),
            pick = ax -> (gui_ref[].controls.h.handles[1].plots[1], 0))
        gui_ref[] = gui
        # every movable object once, the objects of the group after the group, indented
        @test first.(gui.menu.options[]) == ["M1", "Detector", "ObjectGroup", "  G1", "  Mirror", "Beam"]
        @test all(gui.menu_objects .=== Any[m, pd, g, g.objects[1], g.objects[2], beam])
        @test gui.menu.i_selected[] == 0

        # the menu selects, the selection in the 3D view updates the menu
        gui.menu.i_selected[] = 4
        @test gui.controls.selected[] === g.objects[1]
        @test !isempty(gui.controls.box_obs[])
        @test startswith(gui.status.text[], "G1 at (")
        _select!(gui)
        @test gui.controls.selected[] === m
        @test gui.menu.i_selected[] == 1
        _key!(gui, Keyboard.escape)
        @test gui.menu.i_selected[] == 0

        # keys are ignored while the menu is open, e.g. for its search
        gui.menu.is_open[] = true
        _key!(gui, Keyboard.m)
        @test gui.controls.mode[] == :move
        gui.menu.is_open[] = false

        # ray picking of the mirror in the center of the view
        set_view(gui.ax, [0.3, -0.2, 0.3], [0.0, 0.1, 0.0], [0.0, 0, 1])
        scene = gui.ax.scene
        vp = scene.viewport[]
        events(scene).mouseposition[] = (vp.origin[1] + vp.widths[1] / 2, vp.origin[2] + vp.widths[2] / 2)
        @test first(Ext._ray_pick(gui.controls, scene)) === m

        # hide
        notify(gui.hide_button.clicks)
        @test startswith(gui.status.text[], "select a component")
        gui.menu.i_selected[] = 1
        notify(gui.hide_button.clicks)
        @test isnothing(gui.controls.selected[])
        @test gui.menu.i_selected[] == 0
        @test m in gui.hidden
        @test all(p -> !p.visible[], gui.controls.h.handles[1].plots)
        @test isnothing(first(Ext._ray_pick(gui.controls, scene)))
        # also not via the `pick` function
        _select!(gui)
        @test isnothing(gui.controls.selected[])
        # still traced
        Ext._trace!(gui)
        @test length(BMO.hits(pd)) == 1

        # a hidden object can be selected in the menu and shown again
        gui.menu.i_selected[] = 1
        @test gui.controls.selected[] === m
        notify(gui.hide_button.clicks)
        @test !(m in gui.hidden)
        @test all(p -> p.visible[], gui.controls.h.handles[1].plots)
        @test gui.controls.selected[] === m

        # groups: all objects, show all
        gui.menu.i_selected[] = 3
        notify(gui.hide_button.clicks)
        gui.menu.i_selected[] = 1
        notify(gui.hide_button.clicks)
        @test length(gui.hidden) == 3
        leaf_plots = [p for oh in gui.controls.h.handles if oh.obj in (m, g.objects...) for p in oh.plots]
        @test all(p -> !p.visible[], leaf_plots)
        notify(gui.show_all_button.clicks)
        @test isempty(gui.hidden)
        @test all(p -> p.visible[], leaf_plots)
        @test first(Ext._ray_pick(gui.controls, scene)) === m
        close(gui)
    end

    _tick!(gui, dt = 1.0) = (events(gui.ax.scene).tick[] = Makie.Tick(Makie.RegularRenderTick, 0, 0.0, dt))
    _source() = CollimatedSource([0.0, 0, 0], [0.0, 1, 0], 2e-3, 1e-6; num_rings = 2, num_rays = 40)
    _traced(b) = !isnothing(BMO.intersection(first(BMO.rays(b))))

    @testset "preview tracing" begin
        m, pd = _fixture()
        cs = _source()
        n_calls = Ref(0)
        gui = _live_view(System([m, pd]) => cs; throttle = false, fine_step = 1e-4, idle_delay = 10.0,
            on_change = (g, obj) -> (n_calls[] += 1))
        @test gui.preview_enabled
        @test length(BMO.hits(pd)) == 40
        @test !endswith(gui.panels[1].ax.title[], "(preview)")
        n0 = n_calls[]
        gui.controls.selected[] = m
        _key!(gui, Keyboard.up)
        # only the rendered beams are traced, the others are reset
        @test gui.preview
        @test findall(_traced, BMO.beams(cs)) == collect(1:5:40)
        @test length(BMO.hits(pd)) == 8
        @test endswith(gui.panels[1].ax.title[], "8 rays (preview)")
        @test n_calls[] == n0
        pts = copy(gui.beam_handles[1].points[])
        # still moving
        _tick!(gui)
        @test gui.preview
        # the full solve once the movement pauses
        gui.last_change -= 100
        _tick!(gui)
        @test !gui.preview
        @test all(_traced, BMO.beams(cs))
        @test length(BMO.hits(pd)) == 40
        @test gui.panels[1].ax.title[] == "Detector 1: 40 rays"
        @test n_calls[] == n0 + 1
        # the beam plot shows the same rendered subset
        @test gui.beam_handles[1].points[] == pts
        # t solves fully as well
        _key!(gui, Keyboard.up)
        @test gui.preview
        _key!(gui, Keyboard.t)
        @test !gui.preview
        @test length(BMO.hits(pd)) == 40
        close(gui)

        # the budget applies to the preview solve while moving
        m, pd = _fixture()
        gui = _live_view(System([m, pd]) => _source(); throttle = false, fine_step = 1e-4,
            trace_budget = 0.5)
        gui.solve_time = 1.0 # pretend that the full solve is slow
        gui.controls.selected[] = m
        _key!(gui, Keyboard.up)
        @test !gui.pending
        @test gui.preview
        @test gui.preview_time < 0.5
        close(gui)

        # without preview tracing, all beams are traced immediately
        m, pd = _fixture()
        cs = _source()
        gui = _live_view(System([m, pd]) => cs; throttle = false, fine_step = 1e-4, preview = false)
        gui.controls.selected[] = m
        _key!(gui, Keyboard.up)
        @test !gui.preview
        @test all(_traced, BMO.beams(cs))
        @test length(BMO.hits(pd)) == 40
        close(gui)

        # a group rendered completely, single beams and Gaussian beamlets are not affected
        for beam in (_source(), Beam([0.0, 0, 0], [0.0, 1, 0]), _gauss())
            m, pd = _fixture()
            gui = _live_view(System([m, pd]) => beam; throttle = false, fine_step = 1e-4,
                beam_kwargs = beam isa BMO.AbstractBeamGroup ? Dict(beam => (; render_every = 1)) : Dict())
            gui.controls.selected[] = m
            _key!(gui, Keyboard.up)
            @test !gui.preview
            @test !endswith(gui.panels[1].ax.title[], "(preview)")
            close(gui)
        end
    end

    @testset "detector metrics" begin
        # synthetic spot
        c = (0.3e-3, -0.2e-3)
        pts = [Point2(c[1] + dx, c[2] + dz) for (dx, dz) in ((1e-3, 0), (-1e-3, 0), (0, 2e-3), (0, -2e-3))]
        mt = Ext._spot_metrics(pts)
        @test mt.n == 4
        @test abs(mt.cx - c[1]) < 1e-9 && abs(mt.cz - c[2]) < 1e-9
        @test abs(mt.rms - sqrt(2.5) * 1e-3) < 1e-9
        @test abs(mt.rmax - 2e-3) < 1e-9

        # ray hits of a panel
        m, pd = _fixture()
        gui = _live_view(System([m, pd]) => _source())
        p = gui.panels[1]
        spots = BMO.spot_diagram(pd)
        cx, cz = sum(q -> q[1], spots) / 40, sum(q -> q[2], spots) / 40
        @test abs(p.metrics.cx - cx) < 1e-9 && abs(p.metrics.cz - cz) < 1e-9
        @test abs(p.metrics.rms - sqrt(sum(q -> (q[1] - cx)^2 + (q[2] - cz)^2, spots) / 40)) < 1e-9
        @test p.metrics.n == 40
        @test p.centroid[] ≈ [Point2f(1e3 * cx, 1e3 * cz)]
        @test startswith(p.ax.subtitle[], "N = 40, c = (")
        @test occursin("rms ", p.ax.subtitle[])
        @test isnothing(p.profiles_ax) && isempty(p.history_axes)
        close(gui)

        # Gaussian beamlet: power and 1/e² radius
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), _gauss())
        p = gui.panels[1]
        @test isapprox(p.metrics.P, optical_power(pd); rtol = 1e-3)
        # waist 0.5 mm at the source, 200 mm to the detector
        zR = π * 0.5e-3^2 / 1e-6
        w = 0.5e-3 * sqrt(1 + (0.2 / zR)^2)
        @test isapprox(p.metrics.wx, w; rtol = 0.02)
        @test isapprox(p.metrics.wz, w; rtol = 0.02)
        @test abs(p.metrics.cx) < 1e-6 && abs(p.metrics.cz) < 1e-6
        @test p.metrics.peak ≈ maximum(BMO.intensity(pd)[3])
        @test startswith(p.ax.subtitle[], "P = ")
        close(gui)

        # logarithmic color scale with a floor, fixed color range
        m, pd = _fixture()
        # an area of ±5 w, where the floor applies
        area = (; n = 20, x_min = -2.5e-3, x_max = 2.5e-3, z_min = -2.5e-3, z_max = 2.5e-3)
        gui = _live_view(System([m, pd]), _gauss(); detectors = [pd => (:intensity, (; area..., colorscale = :log))])
        p = gui.panels[1]
        _, _, I = BMO.intensity(pd; area...)
        Imax = maximum(I)
        @test p.heat_I[] ≈ Float32.(log10.(max.(I, 1e-4 * Imax)))
        @test minimum(p.heat_I[]) ≈ log10(1e-4 * Imax) rtol = 1e-5
        @test collect(p.heat_plot.colorrange[]) ≈ [log10(1e-4 * Imax), log10(Imax)]
        close(gui)
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), _gauss(); detectors = [pd => (:intensity, (; n = 20, colorrange = (0, 5)))])
        @test collect(gui.panels[1].heat_plot.colorrange[]) == [0, 5]
        @test size(gui.panels[1].heat_I[]) == (20, 20)
        close(gui)
        @test_throws ArgumentError _live_view(System([m, pd]), _gauss(); detectors = [pd => (:intensity, (; colorscale = :sqrt))])

        # history: one point per solve, at most 300
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); throttle = false,
            fine_step = 1e-4, detectors = [pd => (:spot, (; history = true))])
        p = gui.panels[1]
        @test length(p.history_axes) == 2
        @test length(p.history_value[]) == 1
        gui.controls.selected[] = m
        _key!(gui, Keyboard.up)
        _key!(gui, Keyboard.up)
        @test length(p.history_value[]) == 3
        @test last(p.history_value[]) == Point2f(3, 1)
        @test last(p.history_cx[]) ≈ Point2f(3, 1e3 * p.metrics.cx)
        @test p.history_axes[1].ylabel[] == "N"
        for _ in 1:310
            Ext._resolve!(gui, nothing)
        end
        @test length(p.history_value[]) == 300
        @test length(p.history_cz[]) == 300
        @test first(p.history_value[])[1] == 14 && last(p.history_value[])[1] == 313
        close(gui)

        # profiles through the centroid
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), _gauss(); detectors = [pd => (:intensity, (; n = 30, profiles = true, history = true))])
        p = gui.panels[1]
        @test !isnothing(p.profiles_ax)
        @test length(p.profiles_ax.scene.plots) == 2
        x, z, I = BMO.intensity(pd; n = 30)
        i, j = argmin(abs.(x .- p.metrics.cx)), argmin(abs.(z .- p.metrics.cz))
        @test p.profile_x[] ≈ [Point2f(1e3 * x[k], I[k, j]) for k in 1:30]
        @test p.profile_z[] ≈ [Point2f(1e3 * z[k], I[i, k]) for k in 1:30]
        @test p.history_axes[1].ylabel[] == "P [mW]"
        @test last(p.history_value[])[2] ≈ 1e3 * p.metrics.P
        close(gui)
    end

    # Moves the cursor onto the projection of the point `p` [m] in the 3D view
    function _cursor_to!(gui, p)
        scene = gui.ax.scene
        px = Makie.project(scene, :data, :pixel, Point3(p))
        vp = scene.viewport[]
        events(scene).mouseposition[] = (vp.origin[1] + px[1], vp.origin[2] + px[2])
        return nothing
    end

    @testset "beam inspection" begin
        # lens in front of the mirror, such that the refractive index matters
        m, pd = _fixture()
        lens = SphericalLens(0.05, -0.05, 5e-3, 25.4e-3)
        translate3d!(lens, [0, 0.05, 0])
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        gui = _live_view(System([lens, m, pd]), beam; throttle = false)
        rays = collect(BMO.rays(beam))
        @test length(rays) >= 4
        # the segment behind the lens, viewed from above, and the segment towards the detector
        for (k, eye) in ((3, [0.1, 0.06, 0.25]), (4, [0.05, 0.0, 0.25]))
            r = rays[k]
            mid = collect(BMO.position(r)) .+ 0.5 * length(r) .* collect(BMO.direction(r))
            set_view(gui.ax, eye, mid, [0.0, 1, 0])
            _cursor_to!(gui, mid)
            _select!(gui)
            info = gui.inspection
            @test !isnothing(info)
            @test maximum(abs.(info.point .- mid)) < 1e-6
            @test info.direction ≈ collect(BMO.direction(r))
            L = sum(length, rays[1:(k - 1)]) + 0.5 * length(r)
            opl = sum(q -> length(q) * BMO.refractive_index(q), rays[1:(k - 1)]) +
                  0.5 * length(r) * BMO.refractive_index(r)
            @test abs(info.length - L) < 1e-6
            @test abs(info.opl - opl) < 1e-6
            @test info.opl > info.length + 1e-3 # through the lens
            @test isnothing(info.w) && isnothing(info.R)
            @test startswith(gui.status.text[], "beam at (")
            @test occursin("OPL", gui.status.text[])
            @test !isnothing(gui.inspection_plot)
            @test isnothing(gui.controls.selected[])
        end
        # esc removes the marker
        plot = gui.inspection_plot
        _key!(gui, Keyboard.escape)
        @test isnothing(gui.inspection) && isnothing(gui.inspection_plot)
        @test !any(q -> q === plot, gui.ax.scene.plots)
        # a click elsewhere removes it as well
        r = rays[4]
        mid = collect(BMO.position(r)) .+ 0.5 * length(r) .* collect(BMO.direction(r))
        _cursor_to!(gui, mid)
        _select!(gui)
        @test !isnothing(gui.inspection)
        _cursor_to!(gui, mid .+ [0.0, 0.03, 0])
        _select!(gui)
        @test isnothing(gui.inspection) && isnothing(gui.inspection_plot)

        # components win over beams: the beam passes the center of the mirror
        set_view(gui.ax, [0.1, 0.0, 0.25], collect(BMO.position(m)), [0.0, 1, 0])
        _cursor_to!(gui, collect(BMO.position(m)))
        _select!(gui)
        @test gui.controls.selected[] === m
        @test isnothing(gui.inspection)
        close(gui)

        # Gaussian beamlet: w and R of gauss_parameters
        m, pd = _fixture()
        gauss = _gauss()
        gui = _live_view(System([m, pd]), gauss; throttle = false)
        r = BMO.rays(gauss.chief)[1]
        mid = collect(BMO.position(r)) .+ 0.4 * length(r) .* collect(BMO.direction(r))
        set_view(gui.ax, [0.1, 0.04, 0.2], mid, [0.0, 1, 0])
        _cursor_to!(gui, mid)
        _select!(gui)
        info = gui.inspection
        @test maximum(abs.(info.point .- mid)) < 1e-6
        @test abs(info.length - 0.4 * length(r)) < 1e-6
        w, R = BMO.gauss_parameters(gauss, info.length)
        @test info.w == w && info.R == R
        @test occursin("w = ", gui.status.text[]) && occursin("R = ", gui.status.text[])
        close(gui)
    end

    @testset "measuring" begin
        m1 = RoundPlanoMirror(25e-3, 5e-3)
        zrotate3d!(m1, deg2rad(45))
        translate3d!(m1, [0, 0.1, 0])
        m2 = RoundPlanoMirror(25e-3, 5e-3)
        zrotate3d!(m2, deg2rad(-30))
        xrotate3d!(m2, 0.1)
        translate3d!(m2, [0.2, 0.1, 0.05])
        pick_plot = Ref{Any}(nothing)
        gui = _live_view(System([m1, m2]), Beam([0.0, 0, 0], [0.0, 1, 0]); throttle = false,
            labels = Dict(m1 => "M1", m2 => "M2"), pick = ax -> (pick_plot[], 0))
        _plot(obj) = gui.controls.h.handles[findfirst(oh -> oh.obj === obj, gui.controls.h.handles)].plots[1]
        gui.measure_toggle.active[] = true
        @test startswith(gui.status.text[], "measure:")
        pick_plot[] = _plot(m1)
        _select!(gui)
        @test length(gui.measure_points) == 1
        @test occursin("click the second point", gui.status.text[])
        pick_plot[] = _plot(m2)
        _select!(gui)
        d = gui.measurement
        @test abs(d.distance - norm(collect(BMO.position(m2)) - collect(BMO.position(m1)))) < 1e-9
        n1, n2 = BMO.orientation(m1)[:, 2], BMO.orientation(m2)[:, 2]
        @test d.angle ≈ acos(dot(n1, n2))
        @test occursin("M1 to M2: distance", gui.status.text[])
        @test occursin("angle", gui.status.text[])
        @test length(gui.measure_plots) == 2
        @test gui.measure_plots[1] isa Makie.Lines
        # a third point starts a new measurement, the toggle clears it
        _select!(gui)
        @test length(gui.measure_points) == 1 && isnothing(gui.measurement)
        gui.measure_toggle.active[] = false
        @test isempty(gui.measure_points) && isempty(gui.measure_plots)
        close(gui)
    end

    @testset "camera tools" begin
        m, pd = _fixture()
        lens = SphericalLens(0.05, -0.05, 5e-3, 25.4e-3)
        translate3d!(lens, [0, 0.05, 0])
        views = ["top" => ([0.0, 0.05, 0.5], [0.0, 0.05, 0.0], [0.0, 1.0, 0.0])]
        gui = _live_view(System([lens, m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); throttle = false, views)
        cam = cameracontrols(gui.ax.scene)
        _view() = Ext._current_view(gui)
        # the home view is taken at the first tick, e.g. after `set_view` before `display`
        set_view(gui.ax, [0.3, -0.2, 0.25], [0.0, 0.05, 0.0], [0.0, 0.0, 1.0])
        _tick!(gui)
        home = _view()
        @test gui.home_set

        # g: zoom to the selected lens, the view direction is kept
        eye0, lookat0, _ = _view()
        o0 = normalize(eye0 - lookat0)
        gui.controls.selected[] = lens
        _key!(gui, Keyboard.g)
        _tick!(gui)
        bb = Ext._selection_bbox(gui.controls, lens, Ext._object_plots(gui.controls.h, lens))
        center = collect(minimum(bb) .+ Makie.widths(bb) ./ 2)
        d = norm(collect(Makie.widths(bb)))
        eye, lookat, _ = _view()
        @test lookat ≈ center atol = 1e-6
        @test normalize(eye - lookat) ≈ o0 atol = 1e-5
        @test norm(eye - lookat) ≈ (d / 2) / sind(cam.fov[] / 2) rtol = 1e-5
        # animated
        gui.controls.selected[] = nothing
        _key!(gui, Keyboard.g)
        _tick!(gui, 0.05)
        @test !isnothing(gui.camera_animation)
        _tick!(gui)
        @test isnothing(gui.camera_animation)
        eye, lookat, _ = _view()
        bbs = [Makie.boundingbox(p) for h in gui.system_handles for oh in h.handles for p in oh.plots]
        bb = reduce(Ext.GeometryBasics.union, bbs)
        @test lookat ≈ collect(minimum(bb) .+ Makie.widths(bb) ./ 2) atol = 1e-6

        # home
        notify(gui.home_button.clicks)
        _tick!(gui)
        for (a, b) in zip(_view(), home)
            @test a ≈ b atol = 1e-6
        end

        # views menu
        @test first.(gui.views_menu.options[]) == ["top"]
        gui.views_menu.i_selected[] = 1
        _tick!(gui)
        for (a, b) in zip(_view(), views[1].second)
            @test a ≈ b atol = 1e-6
        end
        # save view: appended and printed as Julia code
        set_view(gui.ax, [0.1, 0.2, 0.3], [0.0, 0.01, 0.02], [0.0, 0.0, 1.0])
        out = mktemp() do path, io
            redirect_stdout(io) do
                notify(gui.save_view_button.clicks)
            end
            flush(io)
            read(path, String)
        end
        @test startswith(out, "\"view 2\" => ([")
        entry = eval(Meta.parse(out))
        @test entry.first == "view 2"
        for (a, b) in zip(entry.second, _view())
            @test a == b
        end
        @test first.(gui.views_menu.options[]) == ["top", "view 2"]
        @test gui.views[2] == entry
        close(gui)

        @test_throws ArgumentError _live_view(System([lens]), Beam([0.0, 0, 0], [0.0, 1, 0]); views = ["a" => [1, 2, 3]])
    end
end

end # module
