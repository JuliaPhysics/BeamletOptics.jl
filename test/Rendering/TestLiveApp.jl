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
        @test isnothing(gui.widgets.menu)
        @test gui.sliders isa Makie.SliderGrid
        @test first.(gui.layout.sections[:left]) == ["Objects", "Parameters"]
        @test first.(gui.layout.sections[:right]) == ["Properties"]
        @test first.(gui.layout.dock_panels) == ["Detector 1"]
        @test first.(gui.layout.groups) == [:trace, :camera, :display, :tools, :panels, :help]
        @test occursin("1 ray", gui.widgets.info.text[])
        @test occursin("perspective", gui.widgets.info.text[])
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
        @test gui.widgets.menu isa Makie.Menu
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

        # the colors of the 3D view and of the panels follow the theme
        _rgb(c) = RGBf(Makie.to_color(c))
        _plots(gui, obj) = only(oh for oh in gui.controls.h.handles if oh.obj === obj).plots
        _detector_color(gui, pd) = _rgb(first(p for p in _plots(gui, pd) if p isa Makie.Mesh).color[])
        m, pd = _fixture()
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        dark = _live_app(System([m, pd]), beam; theme = :dark,
            clip_planes = [[0, 0.05, 0] => [0, 1, 0]], detectors = [pd => (:intensity, (; profiles = true))])
        @test dark.beam_handles[1].plot.color[] == t.rays
        @test _detector_color(dark, pd) == t.materials[:detector]
        # the mirror keeps the color of the look
        @test _rgb(first(_plots(dark, m)).color[]) == Ext._materials()[:reflective].color
        plane = only(dark.clip.planes)
        @test only(p for p in _plots(dark, plane) if p isa Makie.Lines).color[] == t.clip_plane
        @test all(p -> p.strokecolor[] == t.marker_stroke,
            filter(p -> p isa Makie.Scatter, [_plots(dark, plane); _plots(dark, beam)]))
        help = only(p for p in dark.controls.plots if p isa Makie.Text && p.parent === dark.ax.blockscene)
        @test help.color[] == t.help
        # the floating cards in the colors of the theme, with a border
        @test dark.cards.selection.background.strokevisible[] && dark.cards.selection.background.strokecolor[] == t.border
        @test dark.cards.selection.background.color[] == t.sidebar
        @test dark.widgets.trace_button.plots[3].backgroundcolor[] == t.tooltip
        px, pz = filter(p -> p isa Makie.Lines, only(dark.panels).profiles_ax.scene.plots)
        @test (px.color[], pz.color[]) == (t.gizmo[1], t.gizmo[3])
        close(dark)
        # the light theme keeps the colors of the compact layout
        m, pd = _fixture()
        light = _live_app(System([m, pd]), beam)
        @test _detector_color(light, pd) == Ext._materials()[:detector].color
        @test _rgb(light.beam_handles[1].plot.color[]) == _rgb(:blue)
        # the cards in the light colors, with a border
        @test light.cards.selection.background.color[] == Ext._app_theme(:light).sidebar
        close(light)
    end

    @testset "toolbar drives the shared logic" begin
        m, pd = _fixture()
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        gui = _live_app(System([m, pd]), beam; throttle = false, clip_planes = [[0, 0.1, 0] => [0, 1, 0]])
        # flat icon buttons and toggles
        @test gui.widgets.trace_button isa Ext._IconButton
        @test gui.widgets.auto_trace_toggle isa Ext._IconToggle
        @test gui.layout.clip_toggle isa Ext._IconToggle
        @test gui.layout.fit_button isa Ext._IconButton
        @test all(t -> t isa Ext._IconToggle, gui.layout.collapse)
        @test gui.widgets.orthographic_toggle.tooltip[] == "Orthographic"
        @test gui.layout.clip_toggle.tooltip[] == "Clipping (c)"
        # orthographic, via a click on the toggle
        cam = cameracontrols(gui.ax.scene)
        @test cam.settings.projectiontype[] == Makie.Perspective
        bb = gui.widgets.orthographic_toggle.box.layoutobservables.computedbbox[]
        events(gui.fig.scene).mouseposition[] = Tuple(Makie.origin(bb) .+ Makie.widths(bb) ./ 2)
        for action in (Mouse.press, Mouse.release)
            events(gui.fig.scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, action)
        end
        @test gui.widgets.orthographic_toggle.active[]
        @test cam.settings.projectiontype[] == Makie.Orthographic
        gui.widgets.orthographic_toggle.active[] = false
        # the views icon opens the views menu
        gui.layout.views_button.clicks[] += 1
        @test gui.widgets.views_menu.is_open[]
        gui.layout.views_button.clicks[] += 1
        @test !gui.widgets.views_menu.is_open[]
        gui.widgets.orthographic_toggle.active[] = true
        @test cam.settings.projectiontype[] == Makie.Orthographic
        @test occursin("orthographic", gui.widgets.info.text[])
        gui.widgets.orthographic_toggle.active[] = false
        @test cam.settings.projectiontype[] == Makie.Perspective
        # auto trace and trace button
        @test gui.trace.auto === gui.widgets.auto_trace_toggle.active
        gui.widgets.auto_trace_toggle.active[] = false
        gui.controls.selected[] = m
        _key!(gui, Keyboard.left)
        @test gui.trace.stale
        gui.widgets.trace_button.clicks[] += 1
        @test !gui.trace.stale
        gui.widgets.auto_trace_toggle.active[] = true
        # clipping, in sync with the key `c`
        @test gui.clip.enabled && gui.layout.clip_toggle.active[]
        gui.layout.clip_toggle.active[] = false
        @test !gui.clip.enabled
        @test isempty(gui.ax.scene.theme.clip_planes[])
        _key!(gui, Keyboard.c)
        @test gui.clip.enabled && gui.layout.clip_toggle.active[]
        @test length(gui.ax.scene.theme.clip_planes[]) == 1
        # clip beams
        gui.widgets.clip_beams_toggle.active[] = true
        @test gui.clip.beams
        gui.widgets.clip_beams_toggle.active[] = false
        @test !gui.clip.beams
        # sources
        marker = _marker(gui, beam)
        gui.widgets.sources_toggle.active[] = false
        @test !any(p -> p.visible[], marker.plots)
        gui.widgets.sources_toggle.active[] = true
        @test all(p -> p.visible[], marker.plots)
        # measure
        gui.widgets.measure_toggle.active[] = true
        @test startswith(gui.status.text[], "measure:")
        gui.widgets.measure_toggle.active[] = false
        # help and fit
        help = gui.controls.help_shown
        gui.layout.help_button.clicks[] += 1
        @test gui.controls.help_shown != help
        gui.layout.fit_button.clicks[] += 1
        @test !isnothing(gui.camera.animation)
        # the views menu is in the toolbar
        n = length(gui.camera.views)
        redirect_stdout(() -> (gui.widgets.save_view_button.clicks[] += 1), devnull)
        @test length(gui.camera.views) == n + 1
        @test length(gui.widgets.views_menu.options[]) == n + 1
        close(gui)
    end

    @testset "inspector" begin
        m, pd = _fixture()
        gui = _live_app(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); throttle = false,
            labels = Dict(m => "M1"))
        @test gui.layout.inspector.name.text[] == "No selection"
        gui.controls.selected[] = m
        @test gui.layout.inspector.name.text[] == "M1"
        @test gui.layout.inspector.type.text[] == string(nameof(typeof(m)))
        # the pose rows of the card, docked in the inspector
        card = gui.layout.inspector.card
        @test Ext._card_widget(card, :y).displayed_string[] == "100.0"
        Ext._card_widget(card, :x).stored_string[] = "5"
        @test BMO.position(m)[1] ≈ 5e-3
        @test gui.layout.inspector.mode.selected[] == :move
        gui.widgets.step_box.stored_string[] = "1 mrad"
        @test gui.layout.inspector.mode.selected[] == :rotate
        # the step box is in the inspector, not on a card
        @test gui.widgets.step_box !== gui.cards.selection.step_box
        # no component menu, the eyes of the tree and "show all" in its title replace it
        @test gui.widgets.show_all_button isa Ext._IconButton
        close(gui)
    end

    # A system with a group of two lenses, a static housing, a detector, a source and a clip plane
    function _tree_fixture(; kwargs...)
        l1 = SphericalLens(0.1, -0.1, 4e-3, 25.4e-3)
        l2 = SphericalLens(-0.1, 0.1, 2e-3, 25.4e-3)
        translate3d!(l2, [0, 0.015, 0])
        group = ObjectGroup([l1, l2])
        translate3d!(group, [0, 0.04, 0])
        m, pd = _fixture()
        cube = BMO.CubeMesh(0.01)
        translate3d!(cube, [0.03, 0.1, 0])
        housing = NonInteractableObject(cube)
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        gui = _live_app(System([group, m, pd, housing]), beam; throttle = false,
            labels = Dict(pd => "PD1"), clip_planes = [[0, 0.1, 0] => [0, 0, 1]], kwargs...)
        return gui, (; l1, l2, group, m, pd, housing, beam)
    end

    _rows(gui) = gui.layout.tree.rows
    _row(gui, key) = only(r for r in _rows(gui) if r.key === key)
    _labels(gui) = [r.label for r in _rows(gui)]

    @testset "object tree rows" begin
        gui, o = _tree_fixture()
        tree = gui.layout.tree
        @test tree isa Ext._ObjectTree
        # systems expanded, groups collapsed; sources and clip planes after the systems
        @test _labels(gui) == ["System 1", "ObjectGroup 1", "Mirror 1", "PD1",
            "NonInteractableObject 1", "Beam 1", "Clip plane 1"]
        @test [r.kind for r in _rows(gui)] ==
              [:system, :group, :mirror, :detector, :mesh, :source, :clip_plane]
        @test [r.depth for r in _rows(gui)] == [0, 1, 1, 1, 1, 0, 0]
        @test _row(gui, o.group).expandable && !_row(gui, o.group).expanded
        @test _row(gui, gui.system_handles[1]).expanded
        # eyes for everything that is rendered, none for clip planes
        @test all(r -> r.visible === true, _rows(gui)[1:6])
        @test isnothing(_row(gui, gui.clip.planes[1]).visible)
        # the names are used in the status line and the inspector
        @test Ext._label(gui, o.m) == "Mirror 1"
        # the kinds of the tree, by dispatch
        @test Ext._tree_kind(o.l1) == :lens
        @test Ext._tree_kind(o.pd) == :detector
        @test Ext._tree_kind(gui.system_handles[1]) == :system
        @test Ext._tree_kind(RoundThinBeamsplitter(0.01)) == :beamsplitter
        # one plot per part, independent of the number of rows
        @test length(tree.scene.plots) == 7
        # the compact layout has no tree and its hooks do nothing
        close(gui)
        gui = live_view(System([o.m, o.pd]), o.beam; trace_budget = Inf)
        @test !hasproperty(gui.layout, :tree)
        @test isnothing(Ext._on_clip_planes_changed!(gui))
        @test Ext._label(gui, o.m) == "Mirror 1"
        _key!(gui, Keyboard.p)
        @test gui.labels[gui.clip.planes[1]] == "Clip plane 1"
        close(gui)
    end

    @testset "object tree interaction" begin
        gui, o = _tree_fixture()
        tree, ctrl = gui.layout.tree, gui.controls
        # a click on a row selects like a click in the 3D view
        tree.clicked[] = o.m
        @test ctrl.selected[] === o.m
        @test tree.selected === o.m
        @test gui.layout.inspector.name.text[] == "Mirror 1"
        @test startswith(gui.status.text[], "Mirror 1")
        # objects that are not movable (e.g. static ones, or not among the `objects` kwarg of the
        # controls) are inspected instead of selected, they are listed and can be hidden
        filter!(x -> x !== o.housing, ctrl.movable)
        tree.clicked[] = o.housing
        @test isnothing(ctrl.selected[]) && gui.objects.inspected === o.housing
        @test gui.layout.inspector.name.text[] == "NonInteractableObject 1"
        # a selection in the 3D view (here: of an object in a collapsed group) expands the group
        # and highlights the row
        ctrl.selected[] = o.l2
        @test _row(gui, o.group).expanded
        @test "Lens 2" in _labels(gui)
        @test _row(gui, o.l2).depth == 2
        @test tree.selected === o.l2
        ctrl.selected[] = nothing
        @test isnothing(tree.selected)
        # expanders keep their state per group
        tree.expand_clicked[] = o.group
        @test !_row(gui, o.group).expanded
        @test !("Lens 1" in _labels(gui))
        tree.expand_clicked[] = gui.system_handles[1]
        @test _labels(gui) == ["System 1", "Beam 1", "Clip plane 1"]
        tree.expand_clicked[] = gui.system_handles[1]
        @test length(_rows(gui)) == 7
        # the eye hides and shows, the row is muted
        handle(obj) = only(oh for oh in ctrl.h.handles if oh.obj === obj)
        tree.eye_clicked[] = o.housing
        @test o.housing in gui.objects.hidden
        @test !any(p -> p.visible[], handle(o.housing).plots)
        @test _row(gui, o.housing).visible === false
        @test tree.plots.labels.color[][5] == tree.muted_color
        tree.eye_clicked[] = o.housing
        @test isempty(gui.objects.hidden)
        @test all(p -> p.visible[], handle(o.housing).plots)
        @test _row(gui, o.housing).visible === true
        # hiding a group hides its objects and clears a selection within it
        ctrl.selected[] = o.l1
        tree.eye_clicked[] = o.group
        @test o.l1 in gui.objects.hidden && o.l2 in gui.objects.hidden
        @test isnothing(ctrl.selected[])
        @test _row(gui, o.group).visible === false
        # the system eye hides everything, "show all" shows everything again
        tree.eye_clicked[] = gui.system_handles[1]
        @test all(r -> r.visible === false, _rows(gui)[1:5])
        gui.widgets.show_all_button.clicks[] += 1
        @test isempty(gui.objects.hidden)
        @test all(r -> r.visible !== false, _rows(gui))
        # the source marker
        tree.eye_clicked[] = o.beam
        @test !any(p -> p.visible[], handle(o.beam).plots)
        tree.eye_clicked[] = o.beam
        @test all(p -> p.visible[], handle(o.beam).plots)
        close(gui)
    end

    @testset "object tree clip planes" begin
        gui, o = _tree_fixture()
        tree = gui.layout.tree
        # `p` adds a numbered plane, selected and shown in the tree
        _key!(gui, Keyboard.p)
        plane = gui.clip.planes[end]
        @test _labels(gui)[end] == "Clip plane 2"
        @test tree.selected === plane
        @test gui.layout.inspector.name.text[] == "Clip plane 2"
        # `Delete` removes it, the numbers are not reused
        _key!(gui, Keyboard.delete)
        @test _labels(gui)[end] == "Clip plane 1"
        @test isnothing(tree.selected)
        _key!(gui, Keyboard.p)
        @test _labels(gui)[end] == "Clip plane 3"
        # a click selects a plane
        ctrl = gui.controls
        ctrl.selected[] = nothing
        tree.clicked[] = gui.clip.planes[1]
        @test ctrl.selected[] === gui.clip.planes[1]
        close(gui)
    end

    @testset "collapsing" begin
        m, pd = _fixture()
        gui = _live_app(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]);
            sliders = ["a" => (0:0.1:1, v -> nothing)])
        layout = gui.layout
        w0, h0 = _width(gui.ax), _height(gui.ax)
        box = gui.widgets.step_box
        layout.collapse.right.active[] = false
        @test !layout.right.shown
        @test _width(gui.ax) ≈ w0 + 300
        # hidden and laid out off-screen, where it can not take clicks
        @test !box.blockscene.visible[]
        @test maximum(box.layoutobservables.computedbbox[])[1] < 0
        layout.collapse.left.active[] = false
        @test _width(gui.ax) ≈ w0 + 540
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
        # a new toolbar group after the built-in ones, but before "Help", which stays last
        b = Button(Ext._add_toolbar_entry!(gui, :custom); label = "Mine")
        @test first.(gui.layout.groups[(end - 1):end]) == [:custom, :help]
        @test b in contents(gui.layout.groups[end - 1].second)
        # the groups and separators alternate in the columns of the toolbar
        cols(x) = Makie.GridLayoutBase.gridcontent(x).span.cols
        @test [cols(g.second) for g in gui.layout.groups] == [1:1, 3:3, 5:5, 7:7, 9:9, 11:11, 13:13]
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
