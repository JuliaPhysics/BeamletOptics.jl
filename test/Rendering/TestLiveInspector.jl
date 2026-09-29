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
    _center(r) = Point2f(minimum(r) .+ Makie.widths(r) ./ 2)
    _tick!(gui) = (events(gui.ax.scene).tick[] = Makie.Tick(Makie.RegularRenderTick, 0, 0.0, 1.0))
    # A left click at the figure pixel `xy`: mouse events like those of a window
    function _click!(gui, xy)
        ev = events(gui.fig.scene)
        ev.mouseposition[] = (Float64(xy[1]), Float64(xy[2]))
        for action in (Mouse.press, Mouse.release)
            ev.mousebutton[] = Makie.MouseButtonEvent(Mouse.left, action)
        end
        return nothing
    end

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
        ctrl.selected[] = gui.clip.planes[1]
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
        gui.widgets.step_box.stored_string[] = "1 mrad"
        @test sel[] == :rotate
        close(gui)
    end

    @testset "part of the selection shared with the floating card" begin
        gui, o = _fixture()
        insp = gui.layout.inspector
        # the step box of the inspector is the keyboard step
        @test gui.widgets.step_box === insp.step_box
        gui.controls.selected[] = o.m
        Ext._refresh_inspector!(gui)
        rows = _rows(gui)
        @test rows == Ext._inspector_rows(gui, o.m) && !isempty(rows)
        close(gui)
        # the floating card of the compact layout shows the same rows for the same object
        gui, o = _fixture(layout = :compact)
        c = gui.cards.selection
        gui.controls.selected[] = o.m
        Ext._update_selection_box!(gui.controls)
        notify(c.properties_button.clicks)
        _tick!(gui)
        @test c.list.rows == rows
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
        # the rows of `card_rows`, the same as on the floating cards
        key(c) = map(Ext._layout_key, Ext._declarations(c, o.pd))
        @test length(Ext._declarations(gui.layout.inspector.card, o.pd)[2]) == length(card_rows(o.pd))
        @test key(gui.layout.inspector.card) == key(gui.cards.selection)
        close(gui)
    end

    @testset "docked card" begin
        gui, o = _fixture()
        insp, ctrl = gui.layout.inspector, gui.controls
        # no floating card for the selection, its card is docked in the inspector
        ctrl.selected[] = o.m
        events(gui.ax.scene).tick[] = Makie.Tick(Makie.RegularRenderTick, 0, 0.0, 1.0)
        @test !Ext._selection_card_shown(gui)
        @test !any(c -> c.scene.visible[], gui.cards.all)
        @test insp.name.text[] == "Mirror 1"
        # the actions of the card: hide shows the hint of the tree and clears the selection
        notify(_w(gui, :hide).clicks)
        @test o.m in gui.objects.hidden && isnothing(ctrl.selected[])
        @test occursin("eye", gui.status.text[])
        Ext._toggle_hidden!(gui, o.m)
        @test !(o.m in gui.objects.hidden)
        # a clip plane: flip and remove
        plane = gui.clip.planes[1]
        ctrl.selected[] = plane
        n = plane.dir[:, 2]
        notify(_w(gui, :flip).clicks)
        @test plane.dir[:, 2] ≈ -n
        notify(_w(gui, :remove).clicks)
        @test isempty(gui.clip.planes) && isnothing(ctrl.selected[])
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
        n = length(gui.cards.all)
        ctrl.selected[] = o.m
        # the pin in the header docks a card of the selected object below the inspector, no
        # floating card
        insp.pin.active[] = true
        c = only(insp.pinned)
        @test c.pinned && c.obj === o.m && Ext._is_pinned(gui, o.m) && length(gui.cards.all) == n
        _tick!(gui)
        @test !any(c -> c.scene.visible[], gui.cards.all)
        @test Ext._card_widget(c, :hide) isa Button && Ext._card_widget(c, :x) isa Textbox
        @test c.head.title.text[] == "Mirror 1" && c.head.icon[] === Ext._icon(:mirror)
        @test c.head.pin.active[]
        # it stays when the selection changes, the pin of the inspector follows the selection
        ctrl.selected[] = o.pd
        @test !insp.pin.active[] && only(insp.pinned) === c
        # a second pinned card below the first
        insp.pin.active[] = true
        c2 = insp.pinned[2]
        @test c2.obj === o.pd && Ext._card_widget(c2, :panel_mode) isa Button
        @test maximum(_rect(c2.parent))[2] <= minimum(_rect(c.parent))[2]
        # the first card may have collapsed to make room, see `_fit_pinned!`; expanded again
        c.collapsed && notify(c.head.collapse.clicks)
        @test !c.collapsed
        # the widgets of a pinned card act on its object, not on the selection
        ctrl.selected[] = o.l1
        Ext._card_widget(c, :z).stored_string[] = "5"
        @test BMO.position(o.m)[3] ≈ 5e-3
        @test Ext._card_widget(c, :z).displayed_string[] == "5.0"
        # and take the keyboard
        Ext._card_widget(c, :y).focused[] = true
        @test Ext._typing(gui)
        Ext._card_widget(c, :y).focused[] = false
        # collapsed to the head and the actions, expanded again
        notify(c.head.collapse.clicks)
        @test c.collapsed && isempty(c.rows.content) && Ext._card_widget(c, :hide) isa Button
        @test c.head.collapse.icon[] === Ext._icon(:expand)
        notify(c.head.collapse.clicks)
        @test !c.collapsed && Ext._card_widget(c, :x) isa Textbox
        # the pin of a pinned card unpins it, the next card moves up
        c.head.pin.active[] = false
        @test only(insp.pinned) === c2 && !Ext._is_pinned(gui, o.m)
        @test isempty(c.widgets) && c.head.title.parent === nothing
        # pin and unpin from the inspector
        ctrl.selected[] = o.m
        insp.pin.active[] = true
        @test Ext._is_pinned(gui, o.m) && length(insp.pinned) == 2
        insp.pin.active[] = false
        @test !Ext._is_pinned(gui, o.m) && only(insp.pinned) === c2
        # the keyboard step stays in the inspector
        @test gui.widgets.step_box !== gui.cards.selection.step_box
        close(gui)
    end

    @testset "pinned cards fit into the sidebar" begin
        gui, o = _fixture()
        insp, ctrl = gui.layout.inspector, gui.controls
        ctrl.selected[] = o.m
        # more pinned cards than fit: the older ones collapse, the newest stays expanded
        for obj in (o.m, o.pd, o.bs, o.pd2, o.l1)
            Ext._toggle_pin!(gui, obj)
            # the older cards collapse, then the property list is shortened; only beyond that
            # room, the sidebar overflows
            @test !Ext._overflows(gui) ||
                  (all(c -> c.collapsed, insp.pinned[1:(end - 1)]) && isempty(insp.list.rows))
            @test length(insp.list.rows) <= length(Ext._inspector_rows(gui, ctrl.selected[]))
            @test !insp.pinned[end].collapsed
        end
        @test count(c -> c.collapsed, insp.pinned) >= 1
        # an expanded card stays expanded, others collapse instead
        c = first(insp.pinned)
        @test c.collapsed
        notify(c.head.collapse.clicks)
        @test !c.collapsed && Ext._card_widget(c, :x) isa Textbox
        @test !Ext._overflows(gui) || all(d -> d.collapsed, filter(d -> d !== c, insp.pinned))
        # a collapsed card is not expanded automatically, e.g. after a solve
        collapsed = [d.collapsed for d in insp.pinned]
        Ext._update_inspector!(gui; force = true)
        @test [d.collapsed for d in insp.pinned] == collapsed
        close(gui)
    end

    @testset "pin and unpin by clicks" begin
        # Regression: the press and release of a click on a widget of the sidebar reached the
        # controls of the 3D view, whose release on the "background" cleared the selection, which
        # emptied the inspector before the pin could act
        gui, o = _fixture()
        insp, ctrl = gui.layout.inspector, gui.controls
        ctrl.selected[] = o.m
        _tick!(gui)
        _click!(gui, _center(_rect(insp.pin.box)))
        # the selection and its docked card stay, the pin is active, a card is pinned below
        @test ctrl.selected[] === o.m
        @test insp.shown === o.m && insp.name.text[] == "Mirror 1"
        @test _w(gui, :x) isa Textbox && !isempty(insp.card.rows.content)
        @test insp.pin.active[] && Ext._is_pinned(gui, o.m)
        @test gui.layout.right.shown
        _tick!(gui)
        c = only(insp.pinned)
        @test c.obj === o.m && !any(c -> c.scene.visible[], gui.cards.all)
        # a click on a widget of the docked card keeps the selection, too
        _click!(gui, _center(_rect(_w(gui, :y))))
        @test ctrl.selected[] === o.m && _w(gui, :y).focused[]
        _w(gui, :y).focused[] = false
        # unpinned via the pin of the inspector: the inspector stays
        _click!(gui, _center(_rect(insp.pin.box)))
        @test !Ext._is_pinned(gui, o.m) && !insp.pin.active[]
        @test ctrl.selected[] === o.m && insp.name.text[] == "Mirror 1" && _w(gui, :x) isa Textbox
        @test isempty(insp.pinned)
        # pinned again, then unpinned via the pin of the pinned card; a click on its widgets keeps
        # the selection
        _click!(gui, _center(_rect(insp.pin.box)))
        c = only(insp.pinned)
        _click!(gui, _center(_rect(Ext._card_widget(c, :y))))
        @test ctrl.selected[] === o.m && Ext._card_widget(c, :y).focused[]
        Ext._card_widget(c, :y).focused[] = false
        _click!(gui, _center(_rect(c.head.pin.box)))
        @test isempty(insp.pinned) && !insp.pin.active[]
        @test ctrl.selected[] === o.m && insp.name.text[] == "Mirror 1" && _w(gui, :x) isa Textbox
        # a click on the background of the 3D view still clears the selection
        vp = gui.ax.scene.viewport[]
        _click!(gui, Point2f(minimum(vp) .+ (5, 5)))
        @test isnothing(ctrl.selected[]) && insp.name.text[] == "No selection"
        close(gui)
    end

    @testset "float and dock pinned cards" begin
        gui, o = _fixture()
        insp, ctrl = gui.layout.inspector, gui.controls
        _floating(obj) = Ext._floating_cards(gui, obj)
        _inside(r, v) = all(minimum(r) .>= minimum(v) .- 1.0f-3) && all(maximum(r) .<= maximum(v) .+ 1.0f-3)
        ctrl.selected[] = o.m
        # pinned: docked by default, with the float button in its head
        Ext._toggle_pin!(gui, o.m)
        d = only(insp.pinned)
        @test Ext._is_pinned(gui, o.m) && !Ext._is_floating(gui, o.m) && isempty(_floating(o.m))
        @test d.head.float.icon[] === Ext._icon(:float) && d.head.float.tooltip[] == "Float in the 3D view"
        # floated: out of the sidebar, a floating card next to the object with the dock button
        notify(d.head.float.clicks)
        _tick!(gui)
        @test isempty(insp.pinned) && Ext._is_pinned(gui, o.m) && Ext._is_floating(gui, o.m)
        @test insp.pin.active[] && isempty(d.widgets) && d.head.title.parent === nothing
        c = only(_floating(o.m))
        @test c !== gui.cards.selection && c.scene.visible[] && c.title.text[] == "Mirror 1"
        @test c.icon[] === Ext._icon(:mirror) && Ext._card_widget(c, :x) isa Textbox
        @test c.dock_button.icon[] === Ext._icon(:dock) && c.dock_button.tooltip[] == "Dock in the sidebar"
        # the card of the selection stays docked in the inspector, whose header has no float button
        @test insp.name.text[] == "Mirror 1" && _w(gui, :x) isa Textbox && !gui.cards.selection.scene.visible[]
        _icons(g) = [p.marker[] for gc in g.content if gc.content isa Box
                     for p in gc.content.blockscene.plots if p isa Scatter]
        @test Ext._icon(:pin) in _icons(insp.card.header) || Ext._icon(:pinned) in _icons(insp.card.header)
        @test !(Ext._icon(:float) in _icons(insp.card.header))        # in the 3D view, linked to its object, the dock button in its head
        bg = _rect(c.background)
        @test _inside(bg, Rect2f(gui.ax.scene.viewport[])) && length(c.link[]) == 2
        @test _inside(_rect(c.dock_button.box), bg) && _inside(_rect(c.pin_button.box), bg)
        # in the colors of the theme, like the floating cards of the compact layout
        @test RGBAf(Makie.to_color(c.background.color[])) == RGBAf(Makie.to_color(gui.layout.theme.sidebar))
        # it stays with its object when the selection changes, its widgets act on its object
        ctrl.selected[] = o.pd
        _tick!(gui)
        @test c.scene.visible[] && c.obj === o.m && !insp.pin.active[]
        Ext._card_widget(c, :z).stored_string[] = "3"
        @test BMO.position(o.m)[3] ≈ 3e-3
        @test BMO.position(o.pd)[3] ≈ 0 atol = 1e-12
        # typing into its textbox blocks the keys of the 3D view
        box = Ext._card_widget(c, :y)
        box.focused[] = true
        @test Ext._typing(gui)
        P0 = Vector{Float64}(BMO.position(o.pd))
        events(gui.ax.scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.left, Keyboard.press)
        @test Vector{Float64}(BMO.position(o.pd)) == P0
        box.focused[] = false
        @test !Ext._typing(gui)
        # hidden and shown again: the card stays floating, like in the compact layout
        notify(Ext._card_widget(c, :hide).clicks)
        _tick!(gui)
        @test o.m in gui.objects.hidden && Ext._is_floating(gui, o.m) && c.scene.visible[]
        notify(Ext._card_widget(c, :hide).clicks)
        @test !(o.m in gui.objects.hidden) && only(_floating(o.m)) === c
        # docked again: at the end of the pinned cards, in its collapsed state
        Ext._toggle_pin!(gui, o.pd)
        notify(c.collapse_button.clicks)
        @test c.collapsed
        notify(c.dock_button.clicks)
        _tick!(gui)
        @test !Ext._is_floating(gui, o.m) && !c.pinned && !c.scene.visible[]
        @test [d.obj for d in insp.pinned] == [o.pd, o.m] && insp.pinned[2].collapsed
        @test insp.pinned[2].head.collapse.icon[] === Ext._icon(:expand)
        # floated again, still collapsed; unpinned by the pin of the floating card
        notify(insp.pinned[2].head.float.clicks)
        _tick!(gui)
        c = only(_floating(o.m))
        @test c.collapsed && c.scene.visible[] && [d.obj for d in insp.pinned] == [o.pd]
        c.pin_button.active[] = false
        _tick!(gui)
        @test !Ext._is_pinned(gui, o.m) && isempty(_floating(o.m)) && !c.scene.visible[]
        # pinned again: docked
        Ext._toggle_pin!(gui, o.m)
        @test Ext._is_pinned(gui, o.m) && !Ext._is_floating(gui, o.m)
        @test [d.obj for d in insp.pinned] == [o.pd, o.m]
        # a floating card is unpinned by the pin of the inspector, too
        Ext._float!(gui, o.m)
        c = only(_floating(o.m))
        ctrl.selected[] = o.m
        @test insp.pin.active[]
        insp.pin.active[] = false
        _tick!(gui)
        @test !Ext._is_pinned(gui, o.m) && !c.pinned && !c.scene.visible[]
        @test [d.obj for d in insp.pinned] == [o.pd]
        # nothing to float or dock without a pinned card
        Ext._float!(gui, o.m)
        Ext._dock!(gui, o.pd)
        @test !Ext._is_pinned(gui, o.m) && [d.obj for d in insp.pinned] == [o.pd]
        close(gui)
    end

    @testset "floating cards do not count for the sidebar" begin
        gui, o = _fixture()
        insp, ctrl = gui.layout.inspector, gui.controls
        ctrl.selected[] = o.m
        # an expanded floating card, then more docked cards than fit into the sidebar
        Ext._toggle_pin!(gui, o.m)
        Ext._float!(gui, o.m)
        c = only(Ext._floating_cards(gui, o.m))
        for obj in (o.pd, o.bs, o.pd2, o.l1)
            Ext._toggle_pin!(gui, obj)
            @test !Ext._overflows(gui) ||
                  (all(d -> d.collapsed, insp.pinned[1:(end - 1)]) && isempty(insp.list.rows))
            @test !insp.pinned[end].collapsed
        end
        # the floating card is neither in the stack nor collapsed by the fitting
        @test all(d -> d.obj !== o.m, insp.pinned) && length(insp.pinned) == 4
        @test !c.collapsed && c.pinned
        # floating the docked cards makes room: no docked card is collapsed for them
        for obj in (o.pd, o.bs, o.pd2)
            Ext._float!(gui, obj)
        end
        @test only(insp.pinned).obj === o.l1 && !Ext._overflows(gui)
        @test count(c -> c.pinned, gui.cards.all) == 4
        _tick!(gui)
        # the floating cards cover neither each other nor the view cube
        rects = [_rect(c.background) for c in gui.cards.all if c.scene.visible[]]
        @test length(rects) == 4
        @test !any(Ext._overlaps(rects[i], rects[j]) for i in eachindex(rects) for j in (i + 1):length(rects))
        close(gui)
    end

    @testset "float and dock by clicks" begin
        # The clicks on the float button of a docked card and on the dock button of a floating
        # card keep the selection, like those on the pins
        gui, o = _fixture()
        insp, ctrl = gui.layout.inspector, gui.controls
        ctrl.selected[] = o.m
        _tick!(gui)
        _click!(gui, _center(_rect(insp.pin.box)))
        d = only(insp.pinned)
        _tick!(gui)
        sleep(0.3)  # later than a double click
        _click!(gui, _center(_rect(d.head.float.box)))
        @test ctrl.selected[] === o.m && insp.name.text[] == "Mirror 1" && _w(gui, :x) isa Textbox
        @test isempty(insp.pinned) && Ext._is_floating(gui, o.m) && insp.pin.active[]
        _tick!(gui)
        c = only(Ext._floating_cards(gui, o.m))
        @test c.scene.visible[]
        # a click on a widget of the floating card keeps the selection, too
        sleep(0.3)
        _click!(gui, _center(_rect(Ext._card_widget(c, :y))))
        @test ctrl.selected[] === o.m && Ext._card_widget(c, :y).focused[] && Ext._typing(gui)
        Ext._card_widget(c, :y).focused[] = false
        sleep(0.3)
        _click!(gui, _center(_rect(c.dock_button.box)))
        @test ctrl.selected[] === o.m && insp.name.text[] == "Mirror 1"
        @test !Ext._is_floating(gui, o.m) && only(insp.pinned).obj === o.m
        _tick!(gui)
        @test !c.scene.visible[]
        # floated again and unpinned by a click on the pin of the floating card
        sleep(0.3)
        _click!(gui, _center(_rect(only(insp.pinned).head.float.box)))
        _tick!(gui)
        c = only(Ext._floating_cards(gui, o.m))
        sleep(0.3)
        _click!(gui, _center(_rect(c.pin_button.box)))
        @test ctrl.selected[] === o.m && !Ext._is_pinned(gui, o.m) && !insp.pin.active[]
        @test isempty(insp.pinned) && _w(gui, :x) isa Textbox
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
        gui.widgets.step_box.focused[] = true
        @test Ext._typing(gui) && !cam.selected[]
        gui.widgets.step_box.focused[] = false
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
