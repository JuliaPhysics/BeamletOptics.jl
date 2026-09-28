module TestLiveExtras

using BeamletOptics
using Makie
using Test

const BMO = BeamletOptics

@testset "Live view extras and opacity" begin
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
    # A housing next to the beam path, as a mesh like a `MeshDummy` of an STL file
    function _housing(p = [0.05, 0.05, 0])
        housing = NonInteractableObject(BMO.CubeMesh(0.02))
        translate3d!(housing, p)
        return housing
    end
    _live_view(args...; kwargs...) = live_view(args...; merge((; trace_budget = Inf, throttle = false), kwargs)...)
    _tick!(gui) = (events(gui.ax.scene).tick[] = Makie.Tick(Makie.RegularRenderTick, 0, 0.0, 1.0))
    _handle(gui, obj) = only(oh for oh in gui.controls.h.handles if oh.obj === obj)
    _alphas(gui, obj) = [p.alpha[] for p in _handle(gui, obj).plots]
    _card_names(c) = [w.name for (_, w) in c.widgets]

    @testset "opacity row by type" begin
        m, pd = _fixture()
        housing = _housing()
        lens = SphericalLens(0.1, -0.1, 4e-3, 25.4e-3)
        opacity(obj) = any(r -> any(c -> c isa CardWidget && c.name === :opacity, r.cells), card_rows(obj))
        @test opacity(housing)
        @test opacity(IntersectableObject(BMO.CubeMesh(0.01)))
        @test !opacity(lens)
        @test !opacity(m)
        @test !opacity(pd)
        # the pose rows come first
        @test length(card_rows(housing)) == length(pose_card_rows(housing)) + 1
    end

    @testset "opacity via the card" begin
        m, pd = _fixture()
        housing = _housing()
        gui = _live_view(System([m, pd, housing]), Beam([0.0, 0, 0], [0.0, 1, 0]))
        ctrl = gui.controls
        # the card of the selection shows the slider at the opacity as rendered
        Ext._select!(gui, housing)
        _tick!(gui)
        slider = Ext._card_widget(gui.card, :opacity)
        @test slider isa Slider
        @test slider.value[] == 100
        @test Ext._card_widget(gui.card, :opacity_value).text[] == "100 %"
        @test Ext._opacity(gui, housing) == 1
        plots = _handle(gui, housing).plots
        @test all(p -> !p.transparency[], plots)
        # the slider scales the alpha of the plots, which become transparent
        Makie.set_close_to!(slider, 40)
        @test Ext._opacity(gui, housing) ≈ 0.4
        @test all(a -> a ≈ 0.4f0, _alphas(gui, housing))
        @test all(p -> p.transparency[], plots)
        @test Ext._card_widget(gui.card, :opacity_value).text[] == "40 %"
        @test occursin("opacity 40 %", gui.status.text[])
        # no solve: the optics are unchanged, the beams are not outdated
        @test !gui.stale
        # the other objects keep their look
        @test all(a -> a == 1, _alphas(gui, m))
        # 0 % hides the object like "hide", the selection and its card stay
        Makie.set_close_to!(slider, 0)
        @test housing in gui.hidden
        @test all(p -> !p.visible[], plots)
        @test ctrl.selected[] === housing
        @test Ext._card_widget(gui.card, :hide).label[] == "show"
        # raising the slider shows it again
        Makie.set_close_to!(slider, 60)
        @test !(housing in gui.hidden)
        @test all(p -> p.visible[], plots)
        @test all(a -> a ≈ 0.6f0, _alphas(gui, housing))
        # the opacity survives hiding and showing
        Ext._toggle_hidden!(gui, housing)
        @test housing in gui.hidden
        @test Ext._opacity(gui, housing) ≈ 0.6
        Ext._toggle_hidden!(gui, housing)
        @test !(housing in gui.hidden)
        @test all(a -> a ≈ 0.6f0, _alphas(gui, housing))
        # showing an object hidden at 0 % restores its initial opacity
        Ext._set_opacity!(gui, housing, 0)
        @test housing in gui.hidden
        Ext._show_all!(gui)
        @test !(housing in gui.hidden)
        @test Ext._opacity(gui, housing) == 1
        # 100 % restores the plots as rendered, opaque
        @test all(a -> a == 1, _alphas(gui, housing))
        @test all(p -> !p.transparency[], plots)
        close(gui)
    end

    @testset "opacity of transparent mechanics and edges" begin
        m, pd = _fixture()
        # a mechanical part with an analytic shape, rendered transparent as in
        # `render!(ax, housing; transparency = true, color = ...)`, with feature edges (a mesh of an
        # STL file has none)
        housing = NonInteractableObject(BMO.BoxSDF(0.02, 0.01, 0.01))
        translate3d!(housing, [0.05, 0.05, 0])
        gui = _live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]);
            extras = [housing => (; transparency = true, color = RGBAf(0.7, 0.8, 0.9, 0.2), edges = true)])
        plots = _handle(gui, housing).plots
        @test length(plots) == 2
        @test any(p -> p isa Makie.Lines, plots)
        @test Ext._opacity(gui, housing) ≈ 0.2
        # 100 % makes it opaque: the opaque rendering, the edges scaled with it (at most opaque)
        Ext._set_opacity!(gui, housing, 1)
        @test Ext._opacity(gui, housing) == 1
        faces = only(p for p in plots if !(p isa Makie.Lines))
        edges = only(p for p in plots if p isa Makie.Lines)
        @test Ext._plot_opacity(faces) ≈ 1
        @test !faces.transparency[]
        @test Ext._plot_opacity(edges) <= 1
        # back to 20 %, as rendered
        Ext._set_opacity!(gui, housing, 0.2)
        @test faces.transparency[] && edges.transparency[]
        @test Ext._plot_opacity(faces) ≈ 0.2
        e0 = Ext._plot_opacity(edges)
        # the edges follow the opacity
        Ext._set_opacity!(gui, housing, 0.1)
        @test Ext._plot_opacity(edges) ≈ e0 / 2
        close(gui)
    end

    @testset "extras" begin
        m, pd = _fixture()
        housing = _housing()
        mount = _housing([-0.05, 0.05, 0])
        sys = System([m, pd])
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        gui = _live_view(sys, beam; extras = [housing => (; color = :lightblue), mount],
            labels = Dict(housing => "housing"))
        ctrl = gui.controls
        # rendered, with the render kwargs, but not part of any system
        @test Ext._leaves(gui.extras) == [housing, mount]
        @test !any(o -> o === housing, BMO.objects(sys))
        @test Makie.to_color(_handle(gui, housing).plots[1].color[]) == Makie.to_color(:lightblue)
        @test Ext._is_extra(gui, housing) && !Ext._is_extra(gui, m)
        # selectable like the objects of the systems: via the plot under the cursor, the menu
        @test Ext._pick_leaf(ctrl.h, _handle(gui, housing).plots[1]) === housing
        @test Ext._is_movable(ctrl, housing)
        @test any(o -> o === housing, gui.menu_objects)
        Ext._select!(gui, housing)
        _tick!(gui)
        @test Ext._card_widget(gui.card, :opacity) isa Slider
        # moving an extra does not solve, the systems keep their beams
        hits = length(BMO.hits(pd))
        solve_time = gui.solve_time
        BMO.translate3d!(housing, [0, 0, 0.01])
        ctrl.on_change(housing)
        @test !gui.stale
        @test startswith(gui.status.text[], "housing")
        @test length(BMO.hits(pd)) == hits
        # the move is exported
        code = export_changes(gui; io = devnull)
        @test occursin("# housing (NonInteractableObject)", code)
        @test occursin("translate_to3d!(housing", code)
        # hide and show
        Ext._toggle_hidden!(gui, mount)
        @test mount in gui.hidden
        close(gui)

        # invalid entries
        @test_throws ArgumentError _live_view(sys, beam; extras = [housing => :blue])
        @test_throws ArgumentError _live_view(sys, beam; extras = ["housing"])
        @test_throws ArgumentError _live_view(sys, beam; extras = housing)
        # an object of a system, or twice
        @test_throws ArgumentError _live_view(sys, beam; extras = [m])
        @test_throws ArgumentError _live_view(sys, beam; extras = [housing, housing])
        # no extras: an empty handle, nothing in the tree
        gui = _live_view(sys, beam)
        @test isempty(gui.extras.handles)
        close(gui)
    end

    @testset "extras in the app layout" begin
        m, pd = _fixture()
        housing = _housing()
        group = ObjectGroup([_housing([0, 0.2, 0]), _housing([0, 0.25, 0])])
        gui = _live_view(System([m, pd]), Beam([0.0, 0, 0], [0.0, 1, 0]); layout = :app,
            extras = [housing, group])
        tree = gui.layout.tree
        rows = tree.rows
        labels = [r.label for r in rows]
        @test labels[1:4] == ["System 1", "Mirror 1", "Detector 1", "Extras"]
        extras_row = rows[4]
        @test extras_row.key === gui.extras
        @test extras_row.kind == :group && extras_row.depth == 0 && extras_row.expanded
        @test [(r.label, r.depth) for r in rows[5:6]] == [("NonInteractableObject 1", 1), ("ObjectGroup 1", 1)]
        # a click in the tree selects the extra, its card is docked in the inspector
        tree.clicked[] = housing
        @test gui.controls.selected[] === housing
        card = gui.layout.inspector.card
        slider = Ext._card_widget(card, :opacity)
        @test slider isa Slider
        Makie.set_close_to!(slider, 0)
        @test housing in gui.hidden
        # the eye of the row shows the object hidden at 0 %
        row(key) = only(r for r in tree.rows if r.key === key)
        @test row(housing).visible === false
        tree.eye_clicked[] = housing
        @test !(housing in gui.hidden)
        @test row(housing).visible === true
        @test Ext._opacity(gui, housing) == 1
        # the eye of "Extras" hides all extras
        tree.eye_clicked[] = gui.extras
        @test all(leaf -> leaf in gui.hidden, Ext._leaves(gui.extras))
        @test row(gui.extras).visible === false
        # the group is selected in the tree like one of a system, its objects are revealed
        Ext._select!(gui, BMO.shape(group)[1])
        @test row(group).expanded
        close(gui)
    end
end

end
