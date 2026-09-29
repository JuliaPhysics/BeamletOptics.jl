module TestLiveSelection

using BeamletOptics
using Makie
using Test

const BMO = BeamletOptics

@testset "Live view selection" begin
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

    # Each change is solved immediately, see TestLiveView.jl
    _live_view(args...; kwargs...) = live_view(args...; merge((; trace_budget = Inf), kwargs)...)

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
        @test first.(gui.widgets.menu.options[]) == ["M1", "Detector 1", "ObjectGroup 1", "  G1", "  Mirror 1", "Beam 1"]
        @test all(gui.objects.menu .=== Any[m, pd, g, g.objects[1], g.objects[2], beam])
        @test gui.widgets.menu.i_selected[] == 0

        # the menu selects, the selection in the 3D view updates the menu
        gui.widgets.menu.i_selected[] = 4
        @test gui.controls.selected[] === g.objects[1]
        @test !isempty(gui.controls.box_obs[])
        @test startswith(gui.status.text[], "G1 at (")
        _select!(gui)
        @test gui.controls.selected[] === m
        @test gui.widgets.menu.i_selected[] == 1
        _key!(gui, Keyboard.escape)
        @test gui.widgets.menu.i_selected[] == 0

        # keys are ignored while the menu is open, e.g. for its search
        gui.widgets.menu.is_open[] = true
        _key!(gui, Keyboard.m)
        @test gui.controls.mode[] == :move
        gui.widgets.menu.is_open[] = false

        # ray picking of the mirror in the center of the view
        set_view(gui.ax, [0.3, -0.2, 0.3], [0.0, 0.1, 0.0], [0.0, 0, 1])
        scene = gui.ax.scene
        vp = scene.viewport[]
        events(scene).mouseposition[] = (vp.origin[1] + vp.widths[1] / 2, vp.origin[2] + vp.widths[2] / 2)
        @test first(Ext._ray_pick(gui.controls, scene)) === m

        # hide
        Ext._toggle_hidden!(gui, nothing)
        @test startswith(gui.status.text[], "select a component")
        gui.widgets.menu.i_selected[] = 1
        notify(Ext._card_widget(gui.cards.selection, :hide).clicks)
        @test isnothing(gui.controls.selected[])
        @test gui.widgets.menu.i_selected[] == 0
        @test m in gui.objects.hidden
        @test all(p -> !p.visible[], gui.controls.h.handles[1].plots)
        @test isnothing(first(Ext._ray_pick(gui.controls, scene)))
        # also not via the `pick` function
        _select!(gui)
        @test isnothing(gui.controls.selected[])
        # still traced
        Ext._trace!(gui)
        @test length(BMO.hits(pd)) == 1

        # a hidden object can be selected in the menu and shown again
        gui.widgets.menu.i_selected[] = 1
        @test gui.controls.selected[] === m
        notify(Ext._card_widget(gui.cards.selection, :hide).clicks)
        @test !(m in gui.objects.hidden)
        @test all(p -> p.visible[], gui.controls.h.handles[1].plots)
        @test gui.controls.selected[] === m

        # groups: all objects, show all
        gui.widgets.menu.i_selected[] = 3
        notify(Ext._card_widget(gui.cards.selection, :hide).clicks)
        gui.widgets.menu.i_selected[] = 1
        notify(Ext._card_widget(gui.cards.selection, :hide).clicks)
        @test length(gui.objects.hidden) == 3
        leaf_plots = [p for oh in gui.controls.h.handles if oh.obj in (m, g.objects...) for p in oh.plots]
        @test all(p -> !p.visible[], leaf_plots)
        notify(gui.widgets.show_all_button.clicks)
        @test isempty(gui.objects.hidden)
        @test all(p -> p.visible[], leaf_plots)
        @test first(Ext._ray_pick(gui.controls, scene)) === m
        close(gui)
    end

    @testset "automatic names and clip planes" begin
        m1, pd = _fixture()
        m2 = RoundPlanoMirror(25e-3, 5e-3)
        translate3d!(m2, [0.2, 0.1, 0])
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        gui = _live_view(System([m1, m2, pd]), beam; throttle = false,
            labels = Dict(m2 => "Own name"), detectors = [])
        gui.export_clipboard = false
        # objects without a label get "Type i" like in the tree of the app, labelled ones keep it
        @test Ext._label(gui, m1) == "Mirror 1"
        @test Ext._label(gui, m2) == "Own name"
        @test first.(gui.widgets.menu.options[]) == ["Mirror 1", "Own name", "Detector 1", "Beam 1"]
        app = _live_view(System([m1, m2, pd]), beam; throttle = false, layout = :app,
            labels = Dict(m2 => "Own name"), detectors = [])
        @test [Ext._label(app, o) for o in (m1, m2, pd, beam)] == ["Mirror 1", "Own name", "Detector 1", "Beam 1"]
        close(app)

        # clip planes are numbered
        p1 = Ext._add_clip_plane!(gui, [0, 0.05, 0], [0, 1, 0])
        p2 = Ext._add_clip_plane!(gui, [0, 0.06, 0], [0, 1, 0])
        @test Ext._label(gui, p1) == "Clip plane 1"
        @test Ext._label(gui, p2) == "Clip plane 2"

        # automatic names are no variable names in the export
        translate3d!(m1, [0, 1e-3, 0])
        names = Ext._export_names(gui, Any[m1, m2])
        @test names[m1] == "obj1" && names[m2] == "obj2"
        code = export_changes(gui; io = devnull)
        @test occursin("# Mirror\n", code) || occursin("(Mirror)", code)
        @test !occursin("Mirror 1", code)
        close(gui)
    end

end

end # module
