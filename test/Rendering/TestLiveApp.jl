module TestLiveApp

using BeamletOptics
using Makie
using Test

const BMO = BeamletOptics

@testset "Live view app layout" begin
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

    _key!(gui, key) = (events(gui.ax.scene).keyboardbutton[] = Makie.KeyEvent(key, Keyboard.press))

    # Deterministic solves, see TestLiveView.jl
    _live_app(args...; kwargs...) =
        live_view(args...; merge((; trace_budget = Inf, layout = :app), kwargs)...)

    _width(block) = widths(block.layoutobservables.computedbbox[])[1]
    _height(block) = widths(block.layoutobservables.computedbbox[])[2]
    _marker(gui, src) = only(oh for oh in gui.controls.h.handles if oh.obj === src)

    @testset "construction" begin
        m, pd = _fixture()
        gui = _live_app(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]);
            sliders = ["a" => (0:0.1:1, v -> nothing)])
        @test gui isa Ext.AppView
        @test gui isa Ext.LiveView
        @test gui.layout isa Ext.AppLayout
        @test Tuple(gui.fig.scene.viewport[].widths) == (1600, 950)
        @test length(gui.panels) == 1
        @test gui.layout.dock.shown
        @test isnothing(gui.menu)
        @test gui.sliders isa Makie.SliderGrid
        @test first.(gui.layout.sections[:left]) == ["Objects", "Parameters"]
        @test first.(gui.layout.sections[:right]) == ["Properties"]
        @test first.(gui.layout.dock_panels) == ["Detector 1"]
        @test first.(gui.layout.groups) == [:trace, :camera, :display, :tools, :panels, :help]
        @test occursin("1 ray", gui.layout.info.text[])
        @test occursin("perspective", gui.layout.info.text[])
        @test sprint(show, gui) == "LiveView(1 systems, 1 detector panels)"
        close(gui)

        # without detectors and sliders the dock is collapsed, the 3D view fills the height
        m, pd = _fixture()
        gui = _live_app(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); detectors = [],
            size = (1280, 800))
        @test isempty(gui.panels)
        @test isnothing(gui.sliders)
        @test !gui.layout.dock.shown
        @test !gui.layout.collapse.dock.active[]
        @test first.(gui.layout.sections[:left]) == ["Objects"]
        @test _height(gui.ax) > 650
        close(gui)

        # the compact layout is the default
        gui = live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); trace_budget = Inf)
        @test gui isa Ext.CompactView
        @test gui.menu isa Makie.Menu
        close(gui)

        @test_throws ArgumentError _live_app(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]);
            layout = :fancy)
        @test_throws ArgumentError _live_app(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]);
            theme = :blue)
        @test_throws ArgumentError live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]);
            theme = :blue)
    end

    @testset "dark theme" begin
        m, pd = _fixture()
        gui = _live_app(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); theme = :dark)
        t = Ext._APP_THEMES[:dark]
        @test gui.layout.theme == t
        @test gui.fig.scene.backgroundcolor[] == t.background
        @test gui.ax.scene.backgroundcolor[] == t.view
        @test gui.status.color[] == t.text
        close(gui)
    end

    @testset "toolbar drives the shared logic" begin
        m, pd = _fixture()
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        gui = _live_app(System([m, pd]), beam; throttle = false, clip_planes = [[0, 0.1, 0] => [0, 1, 0]])
        # orthographic
        cam = cameracontrols(gui.ax.scene)
        @test cam.settings.projectiontype[] == Makie.Perspective
        gui.orthographic_toggle.active[] = true
        @test cam.settings.projectiontype[] == Makie.Orthographic
        @test occursin("orthographic", gui.layout.info.text[])
        gui.orthographic_toggle.active[] = false
        @test cam.settings.projectiontype[] == Makie.Perspective
        # auto trace and trace button
        @test gui.auto_trace === gui.auto_trace_toggle.active
        gui.auto_trace_toggle.active[] = false
        gui.controls.selected[] = m
        _key!(gui, Keyboard.left)
        @test gui.stale
        gui.trace_button.clicks[] += 1
        @test !gui.stale
        gui.auto_trace_toggle.active[] = true
        # clipping, in sync with the key `c`
        @test gui.clipping && gui.layout.clip_toggle.active[]
        gui.layout.clip_toggle.active[] = false
        @test !gui.clipping
        @test isempty(gui.ax.scene.theme.clip_planes[])
        _key!(gui, Keyboard.c)
        @test gui.clipping && gui.layout.clip_toggle.active[]
        @test length(gui.ax.scene.theme.clip_planes[]) == 1
        # clip beams
        gui.clip_beams_toggle.active[] = true
        @test gui.clip_beams
        gui.clip_beams_toggle.active[] = false
        @test !gui.clip_beams
        # sources
        marker = _marker(gui, beam)
        gui.sources_toggle.active[] = false
        @test !any(p -> p.visible[], marker.plots)
        gui.sources_toggle.active[] = true
        @test all(p -> p.visible[], marker.plots)
        # measure
        gui.measure_toggle.active[] = true
        @test startswith(gui.status.text[], "measure:")
        gui.measure_toggle.active[] = false
        # help and fit
        help = gui.controls.help_shown
        gui.layout.help_button.clicks[] += 1
        @test gui.controls.help_shown != help
        gui.layout.fit_button.clicks[] += 1
        @test !isnothing(gui.camera_animation)
        # the views menu is in the toolbar
        n = length(gui.views)
        redirect_stdout(() -> (gui.save_view_button.clicks[] += 1), devnull)
        @test length(gui.views) == n + 1
        @test length(gui.views_menu.options[]) == n + 1
        close(gui)
    end

    @testset "inspector" begin
        m, pd = _fixture()
        gui = _live_app(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); throttle = false,
            labels = Dict(m => "M1"))
        @test gui.layout.selection_label.text[] == "No selection"
        gui.controls.selected[] = m
        @test gui.layout.selection_label.text[] == "M1"
        @test gui.layout.type_label.text[] == string(nameof(typeof(m)))
        @test gui.pose_boxes[2].displayed_string[] == "100.0"
        gui.pose_boxes[1].stored_string[] = "5"
        @test BMO.position(m)[1] ≈ 5e-3
        @test gui.layout.mode_label.text[] == "mode: translate"
        gui.step_box.stored_string[] = "1 mrad"
        @test gui.layout.mode_label.text[] == "mode: rotate"
        # hide and show all are in the inspector
        gui.hide_button.clicks[] += 1
        @test m in gui.hidden
        @test isnothing(gui.controls.selected[])
        @test gui.layout.selection_label.text[] == "No selection"
        gui.show_all_button.clicks[] += 1
        @test isempty(gui.hidden)
        close(gui)
    end

    @testset "collapsing" begin
        m, pd = _fixture()
        gui = _live_app(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]);
            sliders = ["a" => (0:0.1:1, v -> nothing)])
        layout = gui.layout
        w0, h0 = _width(gui.ax), _height(gui.ax)
        box = gui.pose_boxes[1]
        layout.collapse.right.active[] = false
        @test !layout.right.shown
        @test _width(gui.ax) ≈ w0 + 260
        # hidden and laid out off-screen, where it can not take clicks
        @test !box.blockscene.visible[]
        @test maximum(box.layoutobservables.computedbbox[])[1] < 0
        layout.collapse.left.active[] = false
        @test _width(gui.ax) ≈ w0 + 500
        @test !gui.sliders.blockscene.visible[]
        layout.collapse.dock.active[] = false
        @test !layout.dock.shown
        @test _height(gui.ax) > h0 + 200
        @test !gui.panels[1].ax.blockscene.visible[]
        # restored
        foreach(t -> t.active[] = true, layout.collapse)
        @test layout.left.shown && layout.right.shown && layout.dock.shown
        @test _width(gui.ax) ≈ w0
        @test _height(gui.ax) ≈ h0
        @test box.blockscene.visible[]
        @test minimum(box.layoutobservables.computedbbox[])[1] > 0
        @test gui.panels[1].ax.scene.visible[]
        close(gui)
    end

    @testset "slots" begin
        m, pd = _fixture()
        gui = _live_app(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); detectors = [])
        # a new toolbar group after the built-in ones
        b = Button(Ext._add_toolbar_entry!(gui, :custom); label = "Mine")
        @test first(gui.layout.groups[end]) == :custom
        @test b in contents(gui.layout.groups[end].second)
        # a sidebar section below the built-in ones
        g = Ext._add_sidebar_section!(gui, :right, "Extra")
        @test g isa GridLayout
        @test first.(gui.layout.sections[:right]) == ["Properties", "Extra"]
        @test_throws ArgumentError Ext._add_sidebar_section!(gui, :top, "Extra")
        # the first dock panel shows the dock
        @test !gui.layout.dock.shown
        d = Ext._add_dock_panel!(gui, "Mine")
        Axis(d[1, 1])
        @test gui.layout.dock.shown
        @test gui.layout.collapse.dock.active[]
        close(gui)

        # the compact layout has no slots yet
        gui = live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); trace_budget = Inf)
        @test_throws ArgumentError Ext._add_toolbar_entry!(gui, :custom)
        @test_throws ArgumentError Ext._add_sidebar_section!(gui, :left, "Extra")
        @test_throws ArgumentError Ext._add_dock_panel!(gui, "Mine")
        close(gui)
    end
end

end
