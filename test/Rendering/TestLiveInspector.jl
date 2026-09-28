module TestLiveInspector

using BeamletOptics
using Makie
using Test

const BMO = BeamletOptics

@testset "Live view inspector" begin
    Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)
    @test !isnothing(Ext)

    # A lens group, a mirror that reflects the beam onto a detector and a beamsplitter off the path
    function _fixture(; kwargs...)
        l1 = SphericalLens(0.1, -0.1, 4e-3, 25.4e-3)
        l2 = SphericalLens(-0.1, 0.1, 2e-3, 25.4e-3)
        translate3d!(l2, [0, 0.015, 0])
        group = ObjectGroup([l1, l2])
        translate3d!(group, [0, 0.04, 0])
        m = RoundPlanoMirror(25e-3, 5e-3)
        zrotate3d!(m, deg2rad(45))
        translate3d!(m, [0, 0.1, 0])
        pd = Detector(5e-3)
        zrotate3d!(pd, -π / 2)
        translate3d!(pd, [0.1, 0.1, 0])
        bs = RoundThinBeamsplitter(10e-3; reflectance = 0.3)
        translate3d!(bs, [0, 0, 0.1])
        pd2 = Detector(5e-3)
        translate3d!(pd2, [0, 0.3, 0.2])
        beam = Beam([0.0, 0, 0], [0.0, 1, 0], 633e-9)
        gui = live_view(System([group, m, pd, bs, pd2]), beam; trace_budget = Inf, layout = :app,
            throttle = false, labels = Dict(pd => "PD1"), detectors = [pd],
            clip_planes = [[0, 0.1, 0] => [0, 0, 1]], kwargs...)
        return gui, (; l1, l2, group, m, pd, bs, pd2, beam)
    end

    _rows(gui) = gui.layout.inspector.list.rows
    _value(gui, label) = (i = findfirst(r -> r[1] == label, _rows(gui));
        isnothing(i) ? nothing : _rows(gui)[i][2])
    # Declared widgets of the docked card by name, see `card_rows`
    _w(gui, name) = Ext._card_widget(gui.layout.inspector.card, name)
    _rect(x) = x.layoutobservables.computedbbox[]

    @testset "formatting" begin
        @test Ext._property_row("Diameter [m]", 25.4e-3) == ("Diameter", "25.4 mm")
        @test Ext._property_row("Wavelength [m]", 633e-9) == ("Wavelength", "633 nm")
        @test Ext._property_row("Offset [m]", -2e-6) == ("Offset", "-2 µm")
        @test Ext._property_row("Length [m]", 1.5) == ("Length", "1.5 m")
        @test Ext._property_row("Angle [rad]", deg2rad(5)) == ("Angle", "87.3 mrad")
        @test Ext._property_row("Power [W]", 1e-3) == ("Power", "0.001 W")
        @test Ext._property_row("Size [m]", [10e-3, 2e-3]) == ("Size", "(10, 2) mm")
        @test Ext._property_row("Optical axis", [0.0, -1.0, 1e-17]) == ("Optical axis", "(0, -1, 0)")
        @test Ext._property_row("n(λ₀)", 1.51680) == ("n(λ₀)", "1.517")
        @test Ext._property_row("Hits", 12) == ("Hits", "12")
        @test Ext._property_row("Stops beams", true) == ("Stops beams", "yes")
        @test Ext._property_row("Shape", "Mesh") == ("Shape", "Mesh")
    end

    @testset "selection" begin
        gui, o = _fixture()
        insp = gui.layout.inspector
        ctrl = gui.controls
        # nothing selected: a summary
        @test insp.name.text[] == "No selection"
        @test _value(gui, "Systems") == "1"
        @test _value(gui, "Objects") == "6"
        @test _value(gui, "Sources") == "1"
        @test _value(gui, "Detector panels") == "1"
        @test _value(gui, "Clip planes") == "1"
        @test endswith(_value(gui, "Last trace"), "ms")
        # without a selection, the docked card is empty and has no height
        @test isempty(insp.card.widgets) && isempty(insp.card.rows.content)
        @test !insp.pin.box.visible[]
        # a lens: header, the rows of its card (the pose) and its properties below them
        ctrl.selected[] = o.l1
        @test insp.name.text[] == "Lens 1"
        @test insp.type.text[] == "Lens"
        @test insp.icon[] === Ext._icon(:lens)
        @test _w(gui, :y).displayed_string[] == "40.0"
        @test _w(gui, :hide) isa Button && _w(gui, :hide).label[] == "hide"
        @test all(k -> _w(gui, k) isa Textbox, (:x, :y, :z, :rx, :ry, :rv))
        @test isnothing(_w(gui, :panel_mode))
        @test insp.pin.box.visible[] && !insp.pin.active[]
        @test _value(gui, "Thickness") == "4 mm"
        @test _value(gui, "n(λ₀)") == "1.5"
        @test isnothing(_value(gui, "Position"))   # in the pose rows
        @test isnothing(_value(gui, "Type"))       # in the header
        # the actions in the header, the rows below it, the properties below the rows
        @test minimum(_rect(insp.card.actions))[1] > maximum(_rect(insp.name))[1] - 1
        @test maximum(_rect(insp.card.rows))[2] < minimum(_rect(insp.type))[2]
        @test maximum(_rect(insp.list.box))[2] < minimum(_rect(insp.card.rows))[2]
        # the boxes fill the width of the sidebar, inside of it
        sidebar = _rect(gui.layout.right.box)
        @test maximum(_rect(_w(gui, :z)))[1] <= maximum(sidebar)[1]
        @test minimum(_rect(_w(gui, :x)))[1] >= minimum(sidebar)[1]
        # a group
        ctrl.selected[] = o.group
        @test insp.name.text[] == "ObjectGroup 1"
        @test insp.icon[] === Ext._icon(:group)
        @test _value(gui, "Parts") == "2"
        # the source
        ctrl.selected[] = o.beam
        @test insp.icon[] === Ext._icon(:source)
        @test _value(gui, "Wavelength") == "633 nm"
        @test _value(gui, "Direction") == "(0, 1, 0)"
        # a beamsplitter, with its own icon
        ctrl.selected[] = o.bs
        @test insp.icon[] === Ext._icon(:beamsplitter)
        @test Ext._tree_kind(o.bs) == :beamsplitter
        @test _value(gui, "Reflectance") == "0.3"
        @test Ext._tree_kind(RoundLinearPolarizer(25e-3, 1e-3, 1e-3, λ -> 1.5)) == :polarizer
        @test Ext._tree_kind(PolarizationFilter(1e-2)) == :polarizer
        # a clip plane
        ctrl.selected[] = gui.clip_planes[1]
        @test insp.name.text[] == "Clip plane 1"
        @test insp.icon[] === Ext._icon(:clip_plane)
        @test _value(gui, "Normal") == "(0, 0, 1)"
        @test isnothing(_w(gui, :hide)) && _w(gui, :flip) isa Button && _w(gui, :remove) isa Button
        # the detector: hits and the options of its panel
        ctrl.selected[] = o.pd
        @test insp.name.text[] == "PD1"
        @test _value(gui, "Hits") == "1"
        @test _value(gui, "Size") == "(5, 5) mm"
        @test _w(gui, :panel_mode).label[] == "auto" && !_w(gui, :panel_log).active[]
        # deselected: the summary again, the widgets of the card are removed
        blocks = copy(insp.card.blocks)
        ctrl.selected[] = nothing
        @test insp.name.text[] == "No selection"
        @test isempty(insp.card.blocks) && isempty(insp.card.widgets)
        @test all(b -> b.parent === nothing, blocks)
        close(gui)
    end

    @testset "updates in place" begin
        gui, o = _fixture()
        insp, ctrl = gui.layout.inspector, gui.controls
        plots = (insp.list.labels, insp.list.values, insp.list.lines)
        n = length(gui.fig.scene.children)
        ctrl.selected[] = o.pd
        blocks = copy(insp.card.blocks)
        @test !isempty(blocks)
        # a move updates the values, the widgets and plots stay
        _w(gui, :z).stored_string[] = "10"
        @test _value(gui, "Hits") == "0"
        @test _w(gui, :z).displayed_string[] == "10.0"
        @test insp.card.blocks == blocks
        @test (insp.list.labels, insp.list.values, insp.list.lines) === plots
        _w(gui, :z).stored_string[] = "0"
        @test _value(gui, "Hits") == "1"
        # a move in the 3D view, too
        translate3d!(o.pd, [0, 0, 0.002])
        Ext._update_inspector!(gui)
        @test _w(gui, :z).displayed_string[] == "2.0"
        translate3d!(o.pd, [0, 0, -0.002])
        # another detector keeps the widgets, too
        ctrl.selected[] = o.pd2
        @test insp.card.blocks == blocks
        @test _w(gui, :panel_mode).label[] == "no panel"
        # the list has a constant number of plots
        ctrl.selected[] = o.l1
        @test (insp.list.labels, insp.list.values, insp.list.lines) === plots
        @test length(insp.list.box.blockscene.plots) == 4   # with the (invisible) box
        # collapsed, the inspector is not updated; it is refreshed when shown again
        gui.layout.collapse.right.active[] = false
        ctrl.selected[] = o.m
        @test insp.shown === o.l1
        gui.layout.collapse.right.active[] = true
        @test insp.shown === o.m
        @test insp.name.text[] == "Mirror 1"
        close(gui)
    end

    @testset "long lists" begin
        gui, o = _fixture()
        list = gui.layout.inspector.list
        rows = [("row $i", string(i)) for i in 1:30]
        Ext._set_rows!(list, rows)
        @test length(list.rows) == Ext._PROPERTY_MAX_ROWS
        @test list.rows[end][1] == "… $(30 - Ext._PROPERTY_MAX_ROWS + 1) more"
        # long values are ellipsized
        Ext._set_rows!(list, [("Shape", "x"^200)])
        @test endswith(list.values.text[][1], "…")
        close(gui)
    end

    @testset "mode" begin
        gui, _ = _fixture()
        ctrl = gui.controls
        sel = gui.layout.inspector.mode.selected
        @test ctrl.mode[] == :move && sel[] == :move
        # the control drives the controls
        gui.layout.inspector.mode.buttons[2].clicks[] += 1
        @test sel[] == :rotate
        @test ctrl.mode[] == :rotate
        # and follows the key `m` and the step box
        events(gui.ax.scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.m, Keyboard.press)
        @test ctrl.mode[] == :move && sel[] == :move
        gui.step_box.stored_string[] = "1 mrad"
        @test sel[] == :rotate
        close(gui)
    end

    @testset "detector panel rows" begin
        gui, o = _fixture()
        ctrl = gui.controls
        p = only(gui.panels)
        ctrl.selected[] = o.pd
        mode, log = _w(gui, :panel_mode), _w(gui, :panel_log)
        @test p.mode == :auto && mode.label[] == "auto"
        # the panel mode cycles, the color scale is switched, both in place
        notify(mode.clicks)
        @test p.mode == :spot && mode.label[] == "spot"
        @test p.scatter_plot.visible[]
        notify(mode.clicks)
        @test p.mode == :intensity && mode.label[] == "intensity"
        @test p.heat_plot.visible[]
        log.active[] = true
        @test p.colorscale == :log
        notify(mode.clicks)
        @test p.mode == :auto
        # the rows show the state of the panel when selected again
        ctrl.selected[] = nothing
        p.colorscale = :linear
        ctrl.selected[] = o.pd
        @test !_w(gui, :panel_log).active[]
        # a detector without a panel: the inputs only show a message
        ctrl.selected[] = o.pd2
        notify(_w(gui, :panel_mode).clicks)
        @test occursin("no panel", gui.status.text[])
        # only in the inspector: the floating cards show the rows of `card_rows`
        @test isempty(Ext._docked_rows(o.m))
        @test length(Ext._declarations(gui.card, o.pd)[2]) == length(card_rows(o.pd))
        @test length(Ext._declarations(gui.layout.inspector.card, o.pd)[2]) == length(card_rows(o.pd)) + 1
        close(gui)
    end

    @testset "docked card" begin
        gui, o = _fixture()
        insp, ctrl = gui.layout.inspector, gui.controls
        # no floating card for the selection, its card is docked in the inspector
        ctrl.selected[] = o.m
        events(gui.ax.scene).tick[] = Makie.Tick(Makie.RegularRenderTick, 0, 0.0, 1.0)
        @test !Ext._selection_card_shown(gui)
        @test !any(c -> c.scene.visible[], gui.cards)
        @test insp.name.text[] == "Mirror 1"
        # the actions of the card: hide shows the hint of the tree and clears the selection
        notify(_w(gui, :hide).clicks)
        @test o.m in gui.hidden && isnothing(ctrl.selected[])
        @test occursin("eye", gui.status.text[])
        Ext._toggle_hidden!(gui, o.m)
        @test !(o.m in gui.hidden)
        # a clip plane: flip and remove
        plane = gui.clip_planes[1]
        ctrl.selected[] = plane
        n = plane.dir[:, 2]
        notify(_w(gui, :flip).clicks)
        @test plane.dir[:, 2] ≈ -n
        notify(_w(gui, :remove).clicks)
        @test isempty(gui.clip_planes) && isnothing(ctrl.selected[])
        @test isempty(insp.card.widgets)
        close(gui)
    end

    @testset "ray slider of sources" begin
        m = RoundPlanoMirror(25e-3, 5e-3)
        zrotate3d!(m, deg2rad(45))
        translate3d!(m, [0, 0.1, 0])
        pd = Detector(5e-3)
        zrotate3d!(pd, -π / 2)
        translate3d!(pd, [0.1, 0.1, 0])
        src = CollimatedSource([0.0, 0, 0], [0.0, 1, 0], 2e-3, 1e-6; num_rings = 2, num_rays = 40)
        gui = live_view(System([m, pd]), src; trace_budget = Inf, layout = :app, preview = false)
        gui.controls.selected[] = src
        rays = _w(gui, :rays)
        @test rays isa Slider && rays.value[] == 40
        @test _w(gui, :ray_count).text[] == "40 rays"
        # the slider fills the sidebar
        @test maximum(_rect(rays))[1] <= maximum(_rect(gui.layout.right.box))[1]
        Makie.set_close_to!(rays, 200)
        @test length(src) == 200 && _w(gui, :ray_count).text[] == "200 rays"
        @test length(BMO.hits(pd)) == 200
        close(gui)
    end

    @testset "pin" begin
        gui, o = _fixture()
        insp, ctrl = gui.layout.inspector, gui.controls
        n = length(gui.cards)
        ctrl.selected[] = o.m
        # the pin in the header pins a floating card to the selected object
        insp.pin.active[] = true
        c = only(filter(c -> c.pinned, gui.cards))
        @test c.obj === o.m && c !== gui.card && length(gui.cards) == n + 1
        events(gui.ax.scene).tick[] = Makie.Tick(Makie.RegularRenderTick, 0, 0.0, 1.0)
        @test c.scene.visible[] && Ext._card_widget(c, :hide) isa Button
        # it stays when the selection changes, the pin follows the selection
        ctrl.selected[] = o.pd
        @test !insp.pin.active[] && c.scene.visible[]
        ctrl.selected[] = o.m
        @test insp.pin.active[]
        # "unpin" on the floating card unpins it, the pin of the inspector follows
        notify(c.pin_button.clicks)
        @test !c.pinned && !c.scene.visible[] && !insp.pin.active[]
        # pin and unpin from the inspector
        insp.pin.active[] = true
        @test Ext._is_pinned(gui, o.m)
        insp.pin.active[] = false
        @test !Ext._is_pinned(gui, o.m)
        # the keyboard step stays in the inspector
        @test gui.step_box !== gui.card.step_box
        close(gui)
    end

    @testset "keyboard" begin
        gui, o = _fixture()
        ctrl = gui.controls
        ev = events(gui.ax.scene)
        cam = cameracontrols(gui.ax.scene)
        ctrl.selected[] = o.m
        # a box of the docked card takes the keyboard: the camera and the controls ignore the keys
        box = _w(gui, :x)
        box.focused[] = true
        @test Ext._typing(gui) && !cam.selected[]
        P0 = Vector{Float64}(BMO.position(o.m))
        ev.keyboardbutton[] = Makie.KeyEvent(Keyboard.left, Keyboard.press)
        @test Vector{Float64}(BMO.position(o.m)) == P0
        box.focused[] = false
        @test !Ext._typing(gui) && cam.selected[]
        # the step box as well
        gui.step_box.focused[] = true
        @test Ext._typing(gui) && !cam.selected[]
        gui.step_box.focused[] = false
        @test !Ext._typing(gui)
        # another selection ends the input into the old box
        box = _w(gui, :x)
        box.focused[] = true
        ctrl.selected[] = o.pd
        @test !box.focused[] && !Ext._typing(gui)
        close(gui)
    end

    @testset "a failing properties method" begin
        gui, o = _fixture()
        m = Module()
        Core.eval(m, quote
            import BeamletOptics
            struct Broken end
            BeamletOptics.properties(::Broken) = error("boom")
        end)
        # errors of `properties` are shown, not thrown
        rows = Base.invokelatest(Ext._inspector_rows, gui, Base.invokelatest(m.Broken))
        @test only(rows)[1] == "Error"
        @test occursin("boom", only(rows)[2])
        close(gui)
    end
end

end
