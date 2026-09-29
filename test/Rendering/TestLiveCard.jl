module TestLiveCard

using BeamletOptics
using Makie
using Test

const BMO = BeamletOptics

# A component with its own card rows, see `card_rows`: its height as a text and a slider that lifts
# it, which changes the optics (`solve = true`); no actions in the head
struct CardTestObject{T, S <: BMO.AbstractShape{T}} <: BMO.AbstractObject{T}
    shape::S
end
BMO.intersect3d(::CardTestObject, ::BMO.AbstractRay) = nothing
BMO.interact3d(::BMO.AbstractSystem, ::CardTestObject, ::BMO.AbstractBeam, ::BMO.AbstractRay) = nothing
_height(o) = round(Int, 1e3 * BMO.position(o)[3])
BMO.card_rows(o::CardTestObject) = (pose_card_rows(o)...,
    CardRow("height", CardWidget(Label; name = :height, value = (gui, o) -> "$(_height(o)) mm")),
    CardRow("lift", CardWidget(Slider; name = :lift, range = 0:5, width = 120, solve = true,
        value = (gui, o) -> _height(o), on = (gui, o, v) -> translate3d!(o, [0, 0, 1e-3 * (v - _height(o))]))))
BMO.card_actions(::CardTestObject) = ()

@testset "Live component card" begin
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
    _gauss() = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 1e-6, 0.5e-3)
    _live_view(args...; kwargs...) = live_view(args...; merge((; trace_budget = Inf), kwargs)...)
    _tick!(gui) = (events(gui.ax.scene).tick[] = Makie.Tick(Makie.RegularRenderTick, 0, 0.0, 1.0))
    # a click on the pin toggle of the card `c`
    _pin!(c) = (c.pin_button.active[] = !c.pin_button.active[])

    _rect(x) = x.layoutobservables.computedbbox[]
    _center(r) = Point2f(minimum(r) .+ Makie.widths(r) ./ 2)
    _away(x) = x.layoutobservables.suggestedbbox[] == Ext._CARD_AWAY
    _move!(gui, xy) = (events(gui.ax.scene).mouseposition[] = (Float64(xy[1]), Float64(xy[2])))
    _press!(gui) = (events(gui.ax.scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press))
    _release!(gui) = (events(gui.ax.scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release))
    _click!(gui, xy) = (_move!(gui, xy); _press!(gui); _release!(gui))
    function _drag!(gui, a, b)
        _move!(gui, a)
        _press!(gui)
        for t in range(0, 1; length = 10)
            _move!(gui, a .+ t .* (b .- a))
        end
        _release!(gui)
    end
    _eye(gui) = Vector{Float64}(cameracontrols(gui.ax.scene).eyeposition[])
    # The card lies inside the 3D view, with the margin
    function _inside(gui)
        r, v = _rect(gui.cards.selection.background), Rect2f(gui.ax.scene.viewport[])
        return all(minimum(r) .>= minimum(v) .+ Ext._CARD_MARGIN .- 1.0f-3) &&
               all(maximum(r) .<= maximum(v) .- Ext._CARD_MARGIN .+ 1.0f-3)
    end
    function _select!(gui, obj)
        gui.controls.selected[] = obj
        Ext._update_selection_box!(gui.controls)
        return nothing
    end
    # Declared widgets by name, the pose boxes by their index
    _w(c, name) = Ext._card_widget(c, name)
    _pose(c, k) = _w(c, (:x, :y, :z, :rx, :ry, :rv)[k])

    @testset "widgets in the card" begin
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), _gauss(); labels = Dict(m => "M1"))
        c = gui.cards.selection
        @test gui.widgets.step_box === c.step_box
        # no textboxes and no hide button below the 3D view
        @test !any(b -> b isa Textbox, gui.fig.content)
        @test !any(b -> b isa Button && b.label[] == "hide", gui.fig.content)
        # drawn after the 3D scene
        @test c.scene.transformation.translation[][3] == Ext._card_z(1)
        # hidden without a selection, without widgets and with all parts away
        @test !c.scene.visible[] && isempty(c.widgets)
        @test all(_away, (c.head, c.rows, c.step, c.actions, c.background))

        _select!(gui, m)
        @test c.scene.visible[]
        @test c.title.text[] == "M1"
        # the default declarations: "hide" and the pose rows
        @test _w(c, :hide).label[] == "hide"
        @test all(k -> _pose(c, k) isa Textbox, 1:6)
        @test !_away(c.rows) && !_away(c.actions) && !_away(c.step)
        @test _inside(gui)
        # the head and the actions in the first line, the rows below, the step below them
        @test minimum(_rect(c.actions))[1] > maximum(_rect(c.head))[1]
        @test maximum(_rect(c.rows))[2] < minimum(_rect(c.head))[2]
        @test maximum(_rect(c.step))[2] < minimum(_rect(c.rows))[2]
        # the pose of the object in the boxes, the boxes of both rows line up
        @test _pose(c, 2).displayed_string[] == "100.0"
        @test minimum(_rect(_pose(c, 1)))[1] ≈ minimum(_rect(_pose(c, 4)))[1]
        # at the bounding box of the object (not at the gizmo), at one of its sides, see "placement"
        set_view(gui.ax, [0.3, -0.2, 0.3], [0.0, 0.1, 0.0], [0.0, 0, 1])
        _tick!(gui)
        sel = Ext._screen_rect(gui.ax.scene, gui.controls.box_obs[], m)
        r = _rect(c.background)
        gap = Ext._CARD_GAP
        @test minimum(r)[1] ≈ maximum(sel)[1] + gap || maximum(r)[1] ≈ minimum(sel)[1] - gap ||
              maximum(r)[2] ≈ minimum(sel)[2] - gap || minimum(r)[2] ≈ maximum(sel)[2] + gap
        # the line from the center of the bounding box to the card, which covers its end
        a, b = c.link[]
        @test a ≈ Ext._link_anchor(gui.ax.scene, gui.controls.box_obs[], m)
        @test all(minimum(sel) .<= a .<= maximum(sel)) && b in r

        # follows the camera
        r0 = _rect(c.background)
        set_view(gui.ax, [0.2, -0.1, 0.1], [0.0, 0.1, 0.0], [0.0, 0, 1])
        _tick!(gui)
        @test _rect(c.background) != r0
        @test _inside(gui)
        # a selection behind the camera: at the edge of the view
        set_view(gui.ax, [0.0, 0.3, 0.0], [0.0, 0.6, 0.0], [0.0, 0, 1])
        _tick!(gui)
        @test all(iszero, Makie.widths(Ext._screen_rect(gui.ax.scene, gui.controls.box_obs[], m)))
        @test c.scene.visible[] && _inside(gui)

        # an open menu would be covered
        gui.widgets.menu.is_open[] = true
        _tick!(gui)
        @test !c.scene.visible[]
        gui.widgets.menu.is_open[] = false
        _tick!(gui)
        @test c.scene.visible[]

        # collapsed to the head and the actions
        h0 = Makie.widths(_rect(c.background))[2]
        notify(c.collapse_button.clicks)
        @test c.collapsed && c.collapse_button.icon[] === Ext._icon(:expand)
        @test _away(c.rows) && _away(c.step) && !_away(c.actions)
        @test Makie.widths(_rect(c.background))[2] < h0 / 2
        notify(c.collapse_button.clicks)
        @test !c.collapsed && c.collapse_button.icon[] === Ext._icon(:collapse)
        @test !_away(c.rows)

        # deselected: hidden
        _select!(gui, nothing)
        @test !c.scene.visible[]
        @test all(_away, (c.head, c.rows, c.step, c.actions, c.background))
        _select!(gui, m)
        close(gui)
        @test !c.scene.visible[]
    end

    @testset "part of the selection: step, mode and properties" begin
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), _gauss(); labels = Dict(m => "M1"), size = (1600, 1000))
        ctrl = gui.controls
        c = gui.cards.selection
        _select!(gui, m)
        _tick!(gui)
        # step and mode below the rows, the properties collapsed by default
        @test c.mode.selected[] == :move && c.mode.keys == [:move, :rotate]
        @test !_away(c.step) && !c.properties_shown && _away(c.properties)
        @test c.properties_button.icon[] === Ext._icon(:expand)
        h0 = Makie.widths(_rect(c.background))[2]
        # expanded: the rows of the inspector below the step, the card grows
        notify(c.properties_button.clicks)
        @test c.properties_shown && !_away(c.properties) && c.properties_button.icon[] === Ext._icon(:collapse)
        @test c.list.rows == Ext._inspector_rows(gui, m) && !isempty(c.list.rows)
        @test maximum(_rect(c.list.box))[2] < minimum(_rect(c.step))[2]
        @test Makie.widths(_rect(c.background))[2] > h0
        # the list fills the card, which contains it
        @test all(minimum(_rect(c.list.box)) .>= minimum(_rect(c.background)))
        @test all(maximum(_rect(c.list.box)) .<= maximum(_rect(c.background)))
        @test _inside(gui)
        # the state stays with the card when another object is selected
        _select!(gui, pd)
        _tick!(gui)
        @test gui.cards.selection === c && c.properties_shown
        @test c.list.rows == Ext._inspector_rows(gui, pd)
        # a click on "Rotate" sets the mode of the controls, the key `m` sets it back
        _click!(gui, _center(_rect(c.mode.buttons[2])))
        @test ctrl.mode[] == :rotate && c.mode.selected[] == :rotate
        @test ctrl.selected[] === pd
        events(gui.ax.scene).keyboardbutton[] = Makie.KeyEvent(Keyboard.m, Keyboard.press)
        @test ctrl.mode[] == :move && c.mode.selected[] == :move
        # a step in an angle unit switches to the rotate mode
        c.step_box.stored_string[] = "1 mrad"
        @test ctrl.mode[] == :rotate && c.mode.selected[] == :rotate
        # collapsed again
        notify(c.properties_button.clicks)
        @test !c.properties_shown && _away(c.properties) && !_away(c.step)

        # a pinned card shows neither step, mode nor properties, also if they were expanded
        notify(c.properties_button.clicks)
        _pin!(c)
        @test c.pinned && c.scene.visible[]
        @test _away(c.step) && _away(c.properties) && !_away(c.rows)
        # the new card of the selection has its own state, collapsed
        _select!(gui, m)
        _tick!(gui)
        c2 = gui.cards.selection
        @test c2 !== c && !c2.properties_shown && !_away(c2.step) && _away(c2.properties)
        @test c2.mode.selected[] == ctrl.mode[]
        @test _away(c.step) && _away(c.properties)
        close(gui)
    end

    @testset "placement" begin
        view = Rect2f(0, 0, 800, 600)
        size = Vec2f(300, 150)
        margin, gap = Ext._CARD_MARGIN, Ext._CARD_GAP
        # right of the selection, top-aligned
        @test Ext._card_position(Rect2f(100, 200, 50, 80), size, view) == Point2f(150 + gap, 280)
        # left of it if the right side has no room
        @test Ext._card_position(Rect2f(600, 200, 50, 80), size, view) == Point2f(600 - gap - 300, 280)
        # inside the view: below its top edge, at the edges for a selection outside of the view
        @test Ext._card_position(Rect2f(100, 500, 50, 200), size, view)[2] == 600 - margin
        @test Ext._card_position(Rect2f(-500, -400, 10, 10), size, view) == Point2f(margin, margin + 150)
        # neither side has room: below, left-aligned, or above
        @test Ext._card_position(Rect2f(50, 200, 700, 80), size, view) == Point2f(50, 200 - gap)
        @test Ext._card_position(Rect2f(50, 100, 700, 80), size, view) == Point2f(50, 180 + gap + 150)
        # no room at all, e.g. an object that fills the view: right, moved into the view
        @test Ext._card_position(Rect2f(0, 0, 800, 600), size, view) == Point2f(800 - margin - 300, 600 - margin)

        # off the view cube
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), _gauss())
        cube = gui.widgets.view_cube
        C = Rect2f(Makie.viewport(cube.scene)[])
        view = Rect2f(gui.ax.scene.viewport[])
        p = Point2f(minimum(C)[1] - 50, maximum(C)[2])
        @test Ext._overlaps(Ext._card_rect(p, size), C)
        q = Ext._avoid(p, size, view, Ext._obstacles(cube))
        r = Ext._card_rect(q, size)
        @test !Ext._overlaps(r, C)
        @test all(minimum(r) .>= minimum(view) .+ margin) && all(maximum(r) .<= maximum(view) .- margin)
        # beside or below the cube, at the margin
        @test minimum(C)[2] - maximum(r)[2] ≈ margin || minimum(C)[1] - maximum(r)[1] ≈ margin
        far = Point2f(minimum(view)[1] + 20, minimum(view)[2] + 200)
        @test Ext._avoid(far, size, view, Ext._obstacles(cube)) == far
        @test Ext._avoid(p, size, view, Rect2f[]) == p
        close(gui)
    end

    @testset "mouse and keyboard" begin
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), _gauss())
        c, ctrl = gui.cards.selection, gui.controls
        ev = events(gui.ax.scene)
        _select!(gui, m)
        @test c.scene.visible[]
        vp = Rect2f(gui.ax.scene.viewport[])
        # a point of the 3D view beside the card and the objects
        beside = Point2f(minimum(vp) .+ 30)
        _move!(gui, beside)
        @test !Ext._over_card(c, ev)
        _move!(gui, _center(_rect(c.title)))
        @test Ext._over_card(c, ev) && ctrl.ignore_mouse()

        # a click on the card keeps the selection
        _click!(gui, _center(_rect(c.title)))
        @test ctrl.selected[] === m
        # a drag from the card moves the card, not the camera, a drag beside it moves the camera
        eye = _eye(gui)
        _drag!(gui, _center(_rect(c.title)), _center(_rect(c.title)) .+ Point2f(150, -100))
        @test _eye(gui) ≈ eye
        @test ctrl.selected[] === m && !isnothing(c.spot)
        ev.mouseposition[] = Tuple(Float64.(_center(_rect(c.title))))
        ev.scroll[] = (0.0, 3.0)
        @test _eye(gui) ≈ eye
        _drag!(gui, beside, beside .+ Point2f(150, 100))
        @test !(_eye(gui) ≈ eye)
        _tick!(gui)

        # the declared widgets get the clicks
        hide = _w(c, :hide)
        xy = _center(_rect(hide))
        _click!(gui, xy)
        @test hide.clicks[] == 1
        @test m in gui.objects.hidden && isnothing(ctrl.selected[])
        @test !c.scene.visible[]
        # hidden widgets are away and take no clicks at their former position
        _click!(gui, xy)
        @test hide.clicks[] == 1 && c.collapse_button.clicks[] == 0
        # selected in the menu, the hidden object shows "show"; the same declarations keep the widgets
        gui.widgets.menu.i_selected[] = findfirst(o -> o === m, gui.objects.menu)
        @test ctrl.selected[] === m && _w(c, :hide) === hide && hide.label[] == "show"
        notify(hide.clicks)
        @test !(m in gui.objects.hidden) && hide.label[] == "hide"

        # a textbox gets the keyboard, a press elsewhere ends the input
        box = _pose(c, 1)
        _click!(gui, _center(_rect(box)))
        @test box.focused[] && Ext._typing(gui)
        P0 = Vector{Float64}(BMO.position(m))
        ev.keyboardbutton[] = Makie.KeyEvent(Keyboard.left, Keyboard.press)
        @test Vector{Float64}(BMO.position(m)) == P0
        _click!(gui, beside)
        @test !box.focused[] && !Ext._typing(gui)
        close(gui)
    end

    @testset "clip planes" begin
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), _gauss())
        c, ctrl = gui.cards.selection, gui.controls
        plane = Ext._add_clip_plane!(gui, [0, 0.05, 0], [0, 1, 0])
        @test ctrl.selected[] === plane
        @test c.title.text[] == "Clip plane 1"
        @test isnothing(_w(c, :hide)) && !isnothing(_w(c, :flip))
        n = plane.dir[:, 2]
        notify(_w(c, :flip).clicks)
        @test plane.dir[:, 2] ≈ -n
        notify(_w(c, :remove).clicks)
        @test !(plane in gui.clip.planes)
        @test isnothing(ctrl.selected[]) && !c.scene.visible[]
        # another object gets its own actions
        _select!(gui, m)
        @test isnothing(_w(c, :flip)) && !isnothing(_w(c, :hide))
        # a pinned plane is unpinned when it is removed
        plane = Ext._add_clip_plane!(gui, [0, 0.05, 0], [0, 1, 0])
        pinned = gui.cards.selection
        _pin!(pinned)
        @test pinned.pinned && pinned.obj === plane
        notify(_w(pinned, :remove).clicks)
        @test !pinned.pinned && !pinned.scene.visible[]
        close(gui)
    end

    @testset "pinned cards" begin
        m, pd = _fixture()
        beam = _gauss()
        # room for three cards next to each other, one of them with the rows of the detector
        gui = _live_view(System([m, pd]), beam; labels = Dict(m => "M1", pd => "PD"), size = (1800, 1100))
        ctrl = gui.controls
        c1 = gui.cards.selection
        # nothing to pin without a selection
        _pin!(c1)
        @test !c1.pinned && !c1.pin_button.active[]

        _select!(gui, m)
        _pin!(c1)
        @test c1.pinned && c1.obj === m && c1.pin_button.active[]
        # the selection gets another card
        @test gui.cards.selection !== c1 && length(gui.cards.all) == 2 && gui.widgets.step_box === gui.cards.selection.step_box
        # the pinned object shows its pinned card only, without the keyboard step
        @test c1.scene.visible[] && !gui.cards.selection.scene.visible[]
        @test _away(c1.step) && !_away(c1.rows)
        @test length(c1.link[]) == 2

        # the pinned card stays when the selection changes
        _select!(gui, pd)
        @test c1.scene.visible[] && gui.cards.selection.scene.visible[]
        @test c1.title.text[] == "M1" && gui.cards.selection.title.text[] == "PD"
        _select!(gui, nothing)
        @test c1.scene.visible[] && !gui.cards.selection.scene.visible[]
        # and follows its object, also when it is moved elsewhere
        set_view(gui.ax, [0.3, -0.2, 0.3], [0.0, 0.1, 0.0], [0.0, 0, 1])
        _tick!(gui)
        r0 = _rect(c1.background)
        translate3d!(m, [0.0, 0.0, 0.02])
        _tick!(gui)
        @test _rect(c1.background) != r0
        @test _pose(c1, 3).displayed_string[] == "20.0"

        # its widgets act on its object, not on the selection
        _select!(gui, pd)
        _pose(c1, 3).stored_string[] = "5"
        @test BMO.position(m)[3] ≈ 5e-3
        @test BMO.position(pd)[3] ≈ 0 atol = 1e-12
        notify(_w(c1, :hide).clicks)
        @test m in gui.objects.hidden && ctrl.selected[] === pd
        @test c1.scene.visible[] && _w(c1, :hide).label[] == "show"
        notify(_w(c1, :hide).clicks)
        @test !(m in gui.objects.hidden) && _w(c1, :hide).label[] == "hide"
        _pose(c1, 1).focused[] = true
        @test Ext._typing(gui)
        _pose(c1, 1).focused[] = false
        @test !Ext._typing(gui)

        # a second pinned card; the widgets of the newest card take the clicks
        _pin!(gui.cards.selection)
        @test length(gui.cards.all) == 3 && count(c -> c.pinned, gui.cards.all) == 2
        _select!(gui, beam)
        c3 = gui.cards.selection
        @test c3 === gui.cards.all[3] && c3.scene.visible[]
        # no card covers another
        rects = [_rect(c.background) for c in gui.cards.all]
        @test !any(Ext._overlaps(rects[i], rects[j]) for i in 1:3 for j in (i + 1):3)
        _click!(gui, _center(_rect(c3.collapse_button.box)))
        @test c3.collapsed
        sleep(0.3)  # later than a double click
        _click!(gui, _center(_rect(c3.collapse_button.box)))
        @test !c3.collapsed
        # the menu hides all cards
        gui.widgets.menu.is_open[] = true
        _tick!(gui)
        @test !any(c -> c.scene.visible[], gui.cards.all)
        gui.widgets.menu.is_open[] = false
        _tick!(gui)
        @test count(c -> c.scene.visible[], gui.cards.all) == 3

        # unpinned: hidden, and reused for the next pin
        _pin!(c1)
        @test !c1.pinned && !c1.scene.visible[] && !c1.pin_button.active[]
        _pin!(c3)
        @test gui.cards.selection === c1 && length(gui.cards.all) == 3
        close(gui)
        @test !any(c -> c.scene.visible[], gui.cards.all)
    end

    @testset "cards moved by the mouse" begin
        # the spot of a card is relative to the nearest corner of the view, inside the view
        view, size = Rect2f(0, 0, 1000, 800), Vec2f(200, 100)
        spot = Ext._card_spot(Point2f(780, 750), size, view)
        @test spot == (true, true, Vec2f(20, 50))
        @test Ext._spot_position(spot, size, view) == Point2f(780, 750)
        @test Ext._spot_position(spot, size, Rect2f(0, 0, 1200, 600)) == Point2f(980, 550)
        @test Ext._spot_position(spot, size, Rect2f(0, 0, 150, 90)) == Point2f(8, 82)
        spot = Ext._card_spot(Point2f(30, 140), size, view)
        @test spot == (false, false, Vec2f(30, 40))
        @test Ext._spot_position(spot, Vec2f(200, 30), view) == Point2f(30, 70)

        m, pd = _fixture()
        gui = _live_view(System([m, pd]), _gauss(); size = (1400, 900))
        c, ctrl = gui.cards.selection, gui.controls
        ev = events(gui.ax.scene)
        vp = Rect2f(gui.ax.scene.viewport[])
        _select!(gui, m)
        _tick!(gui)
        @test isnothing(c.spot)
        top_left(c) = Point2f(minimum(_rect(c.background))[1], maximum(_rect(c.background))[2])
        # a drag at a widget does not move the card
        box = _pose(c, 1)
        p0 = top_left(c)
        _drag!(gui, _center(_rect(box)), _center(_rect(box)) .+ Point2f(-200, 100))
        @test isnothing(c.spot) && top_left(c) == p0
        box.focused[] && Makie.defocus!(box)
        # a drag at the title moves the card to the top left corner of the view, off the objects
        target = Point2f(minimum(vp)[1] + 40, maximum(vp)[2] - 40)
        a = _center(_rect(c.title))
        eye = _eye(gui)
        _drag!(gui, a, a .+ (target - p0))
        @test _eye(gui) ≈ eye && ctrl.selected[] === m
        @test top_left(c) ≈ target
        @test c.spot[1:2] == (false, true)
        # it stays there when the camera moves, with its line to the object
        set_view(gui.ax, [0.3, -0.2, 0.3], [0.0, 0.1, 0.0], [0.0, 0, 1])
        _tick!(gui)
        @test top_left(c) ≈ target && length(c.link[]) == 2
        # and for another selected object
        _select!(gui, pd)
        _tick!(gui)
        @test top_left(c) ≈ target && c.pose[1] === pd
        # a drag beyond the view keeps the card inside the view
        a = _center(_rect(c.title))
        _drag!(gui, a, a .+ Point2f(-500, 500))
        r = _rect(c.background)
        @test minimum(r)[1] ≈ minimum(vp)[1] + Ext._CARD_MARGIN
        @test maximum(r)[2] ≈ maximum(vp)[2] - Ext._CARD_MARGIN
        _drag!(gui, _center(_rect(c.title)), _center(_rect(c.title)) .+ (target - top_left(c)))
        @test top_left(c) ≈ target

        # pinned, the card stays at its spot; the new card of the selection keeps off it
        _pin!(c)
        @test c.pinned && top_left(c) ≈ target
        _select!(gui, m)
        _tick!(gui)
        s = gui.cards.selection
        @test s !== c && isnothing(s.spot) && s.scene.visible[]
        @test !Ext._overlaps(_rect(s.background), _rect(c.background))
        # a double click at the head places the card next to its object again
        a = _center(_rect(c.title))
        _click!(gui, a)
        _click!(gui, a)
        @test isnothing(c.spot) && c.pinned
        @test !(top_left(c) ≈ target)
        # two slow clicks do not
        a = _center(_rect(c.title))
        _drag!(gui, a, a .+ (target - top_left(c)))
        @test !isnothing(c.spot)
        a = _center(_rect(c.title))
        _click!(gui, a)
        sleep(1.2 * Ext._CARD_DOUBLE_CLICK)
        _click!(gui, a)
        @test !isnothing(c.spot)
        # the buttons at the head are no handle: a click on the collapse button collapses the card
        _click!(gui, _center(_rect(c.collapse_button.box)))
        @test c.collapsed && !isnothing(c.spot) && top_left(c) ≈ target
        # an unpinned card forgets its spot
        _pin!(c)
        @test !c.pinned && isnothing(c.spot)
        close(gui)
    end

    @testset "no dock toggle in the compact layout" begin
        # the compact layout has no sidebar: its pinned cards always float, without a dock button
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), _gauss())
        _select!(gui, m)
        c = gui.cards.selection
        _pin!(c)
        _tick!(gui)
        @test c.pinned && c.scene.visible[] && isnothing(c.dock_button)
        @test all(c -> isnothing(c.dock_button), gui.cards.all)
        @test Ext._is_floating(gui, m) && Ext._is_pinned(gui, m)
        # the head holds the pin and the chevron only
        @test length(c.tools.content) == 2
        # floating and docking do nothing
        Ext._dock!(gui, m)
        Ext._float!(gui, m)
        _tick!(gui)
        @test c.pinned && c.obj === m && c.scene.visible[] && count(c -> c.pinned, gui.cards.all) == 1
        close(gui)
    end

    @testset "declared rows" begin
        m, pd = _fixture()
        t = CardTestObject(BMO.shape(RoundPlanoMirror(0.01, 0.002)))
        translate3d!(t, [-0.05, 0.05, 0.0])
        # a second mirror off the beam path, with the same declarations as `m`
        m2 = RoundPlanoMirror(25e-3, 5e-3)
        translate3d!(m2, [0.1, 0.25, 0])
        solves = Ref(0)
        gui = _live_view(System([m, pd, t, m2]), _gauss(); on_change = (gui, obj) -> (solves[] += 1))
        c, ctrl = gui.cards.selection, gui.controls
        _select!(gui, t)
        # the pose rows, then its own rows; no actions
        @test _pose(c, 1) isa Textbox && _w(c, :height).text[] == "0 mm"
        @test isnothing(_w(c, :hide)) && isempty(c.actions.content)
        lift = _w(c, :lift)
        @test lift isa Slider && lift.value[] == 0
        # an input calls `on`, solves again (`solve = true`) and shows the new values on the card
        n = solves[]
        Makie.set_close_to!(lift, 3)
        @test BMO.position(t)[3] ≈ 3e-3
        @test _w(c, :height).text[] == "3 mm" && _pose(c, 3).displayed_string[] == "3.0"
        @test solves[] > n

        # the declared widgets take the clicks like the others: a click at the right end of the
        # slider sets its last value, and neither deselects nor moves the camera
        set_view(gui.ax, [0.3, -0.2, 0.3], [0.0, 0.1, 0.0], [0.0, 0, 1])
        _tick!(gui)
        eye = _eye(gui)
        r = _rect(lift)
        _click!(gui, Point2f(maximum(r)[1] - 2, _center(r)[2]))
        @test lift.value[] == 5 && BMO.position(t)[3] ≈ 5e-3
        @test ctrl.selected[] === t && _eye(gui) ≈ eye

        # another object with other declarations: new widgets, the old ones are removed
        old = copy(c.blocks)
        _select!(gui, m)
        @test isnothing(_w(c, :lift)) && _w(c, :hide) isa Button
        @test all(b -> b.parent === nothing, filter(b -> !any(x -> x === b, c.blocks), old))
        # an object with the same declarations (another mirror) keeps the widgets and shows its values
        kept = copy(c.blocks)
        _select!(gui, m2)
        @test length(c.blocks) == length(kept) && all(c.blocks .=== kept)
        @test _pose(c, 2).displayed_string[] == "250.0"
        close(gui)
    end

    @testset "ray slider of sources" begin
        m, pd = _fixture()
        src = CollimatedSource([0.0, 0, 0], [0.0, 1, 0], 2e-3, 1e-6; num_rings = 2, num_rays = 40)
        # without the preview solve of the rendered rays only, see `preview` of `live_view`
        gui = _live_view(System([m, pd]), src; preview = false)
        c = gui.cards.selection
        _select!(gui, src)
        rays = _w(c, :rays)
        @test rays.value[] == 40 && _w(c, :ray_count).text[] == "40 rays"
        # steps 1-2-5 from the minimum of the rings, 20 × 2, up to 20 000
        @test first(rays.range[]) == 40 && last(rays.range[]) == 20_000 && 200 in rays.range[]
        Makie.set_close_to!(rays, 200)
        @test length(src) == 200 && _w(c, :ray_count).text[] == "200 rays"
        # solved again: the new rays reach the detector
        @test length(BMO.hits(pd)) == 200
        close(gui)
        # a source that wraps given beams has no slider
        wrapped = CollimatedSource(BMO.beams(src), 2e-3, [0.0, 0, 0], [0.0, 1, 0])
        gui = _live_view(System([m, pd]), wrapped)
        _select!(gui, wrapped)
        @test isnothing(_w(gui.cards.selection, :rays)) && _pose(gui.cards.selection, 1) isa Textbox
        close(gui)
    end

    @testset "colors of the theme" begin
        _rgba(x) = RGBAf(Makie.to_color(x))
        _label(c, text) = only(b for b in c.blocks if b isa Label && b.text[] == text)
        for layout in (:compact, :app), theme in (:light, :dark)
            m, pd = _fixture()
            gui = _live_view(System([m, pd]), _gauss(); layout, theme, labels = Dict(m => "M1"))
            t = Ext._app_theme(theme)
            @test gui.layout.theme === t
            # a pinned card: floating in the compact layout, docked below the inspector in the app
            # layout, whose floating cards keep the colors of the theme
            Ext._toggle_pin!(gui, m)
            @test _rgba(gui.cards.selection.background.color[]) == _rgba(t.sidebar)
            @test _rgba(gui.cards.selection.background.strokecolor[]) == _rgba(t.border)
            if layout === :app
                c = only(gui.layout.inspector.pinned)
                @test c.theme === t && !any(c -> c.scene.visible[], gui.cards.all)
                @test _rgba(c.head.title.color[]) == _rgba(t.text) && c.head.title.text[] == "M1"
                @test c.head.pin.active[] && c.head.icon[] === Ext._icon(:mirror)
            else
                c = only(filter(c -> c.pinned, gui.cards.all))
                @test c.theme === t && c.scene.visible[]
                @test _rgba(c.title.color[]) == _rgba(t.text) && c.title.text[] == "M1"
                @test c.pin_button.active[] && c.icon[] === Ext._icon(:mirror)
                # the line to the object
                @test _rgba(c.scene.plots[1].color[]) == _rgba(t.accent)
            end
            # the widgets in the style of the app layout, the rotations in the gizmo colors
            @test _rgba(_pose(c, 1).boxcolor[]) == _rgba(t.field)
            @test _rgba(_pose(c, 1).textcolor[]) == _rgba(t.text)
            @test _rgba(_label(c, "x").color[]) == _rgba(t.text)
            @test [_rgba(_label(c, k).color[]) for k in ("rx", "ry", "rv")] == _rgba.(collect(t.gizmo))
            # the progress window: panel and text from the theme, the bar orange
            panel, track, bar, text = gui.trace.progress.plots
            @test _rgba(panel.color[]) == _rgba(t.sidebar) && _rgba(text.color[]) == _rgba(t.text)
            @test _rgba(track.color[]) == _rgba(t.border) && bar.color[] == Ext._PROGRESS_FILL_COLOR
            close(gui)
        end
        # the icon toggle pins the card of the selection with a click
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), _gauss())
        _select!(gui, m)
        _tick!(gui)
        c = gui.cards.selection
        _click!(gui, _center(_rect(c.pin_button.box)))
        @test c.pinned && c.pin_button.active[] && gui.cards.selection !== c
        sleep(0.3)  # later than a double click
        _click!(gui, _center(_rect(c.pin_button.box)))
        @test !c.pinned && !c.pin_button.active[]
        @test_throws ArgumentError _live_view(System([m, pd]), _gauss(); theme = :blue)
        close(gui)
    end

    @testset "panel rows of detectors" begin
        m, pd = _fixture()
        pd2 = Detector(5e-3)
        translate3d!(pd2, [0, 0.3, 0.2])
        gui = _live_view(System([m, pd, pd2]), _gauss(); detectors = [pd])
        c = gui.cards.selection
        p = only(gui.panels)
        # the pose rows, the beam, the signal, then the mode and the color scale of the panel
        @test length(card_rows(pd)) == length(pose_card_rows(pd)) + 3
        _select!(gui, pd)
        _tick!(gui)
        mode, log = _w(c, :panel_mode), _w(c, :panel_log)
        @test _pose(c, 1) isa Textbox && mode isa Button && log isa Toggle
        @test mode.label[] == "auto" && !log.active[]
        # a click on the button cycles the mode, the toggle switches the color scale
        _click!(gui, _center(_rect(mode)))
        @test p.mode == :spot && mode.label[] == "spot"
        @test gui.controls.selected[] === pd
        log.active[] = true
        @test p.colorscale == :log
        # a detector without a panel: the same widgets, the inputs only show a message
        _select!(gui, pd2)
        _tick!(gui)
        @test _w(c, :panel_mode) === mode && mode.label[] == "no panel" && !log.active[]
        notify(mode.clicks)
        @test occursin("no panel", gui.status.text[]) && p.mode == :spot
        close(gui)
    end

    @testset "info cards of the inspection and the measurement" begin
        _info(y; w = nothing, R = nothing) = (; point = [0.0, y, 0.0], direction = [0.0, 1.0, 0.0],
            length = y, opl = y, w, R)
        _text(c, name) = _w(c, name).text[]
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), _gauss(); labels = Dict(m => "M1", pd => "PD"))
        # an inspected point: a transient card at the point, its pin off, no step box
        Ext._show_inspection!(gui, _info(0.05))
        c = Ext._info_card(gui)
        @test c.transient && c.obj isa Ext._BeamPoint && c.scene.visible[] && !c.pin_button.active[]
        @test c.title.text[] == "Beam" && c.icon[] === Ext._icon(:trace)
        @test _text(c, :at) == "(0.000, 50.000, 0.000) mm" && _text(c, :path) == "50.000 mm"
        @test isnothing(_w(c, :w)) && isnothing(_w(c, :hide)) && _away(c.step)
        @test length(c.link[]) == 2
        # the next inspection replaces it, with the rows of a Gaussian beamlet
        Ext._show_inspection!(gui, _info(0.06; w = 1e-3, R = 0.0))
        @test Ext._info_card(gui) === c && _text(c, :at) == "(0.000, 60.000, 0.000) mm"
        @test _text(c, :w) == "1 mm" && _text(c, :R) == "∞"
        # its pin keeps it, the next inspection gets another card
        c.pin_button.active[] = true
        @test !c.transient && c.pinned && isnothing(Ext._info_card(gui))
        # the card keeps the marker of the point, which the next inspection does not remove
        marker = only(c.obj.plots)
        @test isnothing(gui.measure.inspection_plot) && marker in gui.ax.scene.plots
        Ext._show_inspection!(gui, _info(0.07))
        c2 = Ext._info_card(gui)
        @test c2 !== c && c.scene.visible[] && c2.scene.visible[]
        # clearing the inspection (a click elsewhere, Esc) removes the transient card only
        Ext._clear_inspection!(gui)
        @test isnothing(Ext._info_card(gui)) && !c2.scene.visible[] && c.scene.visible[]
        # unpinned: removed
        c.pin_button.active[] = false
        @test !c.pinned && !c.scene.visible[]
        @test !(marker in gui.ax.scene.plots)

        # a measurement: the first point, then the distance, its components and the angle
        gui.widgets.measure_toggle.active[] = true
        Ext._add_measure_point!(gui, BMO.position(m), m)
        c = Ext._info_card(gui)
        @test c.obj isa Ext._Measurement && c.title.text[] == "Measurement"
        @test _text(c, :at) == "M1, (0.000, 100.000, 0.000) mm"
        Ext._add_measure_point!(gui, BMO.position(pd), pd)
        @test Ext._info_card(gui) === c
        @test _text(c, :from) == "M1" && _text(c, :to) == "PD"
        @test _text(c, :distance) == "100.0 mm" && _text(c, :delta) == "(100.000, 0.000, 0.000) mm"
        @test _text(c, :angle) == "785 mrad"
        # its pin keeps the card with the line and the markers; switching measuring off removes
        # only the transient card of a new measurement
        c.pin_button.active[] = true
        kept = copy(c.obj.plots)
        @test length(kept) == 2 && isempty(gui.measure.plots)
        Ext._add_measure_point!(gui, BMO.position(m), m)
        c3 = Ext._info_card(gui)
        @test c3 !== c && all(p -> p in gui.ax.scene.plots, kept)
        # unpinned: its plots are removed
        c.pin_button.active[] = false
        @test !c.pinned && !any(p -> p in gui.ax.scene.plots, kept)
        c = c3
        gui.widgets.measure_toggle.active[] = false
        @test isnothing(Ext._info_card(gui)) && !c.scene.visible[]
        close(gui)

        # app layout: the transient card floats at the point, its pin docks the item below the
        # inspector
        m, pd = _fixture()
        gui = _live_view(System([m, pd]), _gauss(); layout = :app)
        Ext._show_inspection!(gui, _info(0.05))
        c = Ext._info_card(gui)
        @test c.scene.visible[]
        c.pin_button.active[] = true
        d = only(gui.layout.inspector.pinned)
        @test d.obj isa Ext._BeamPoint && d.head.title.text[] == "Beam"
        @test Ext._card_widget(d, :at).text[] == "(0.000, 50.000, 0.000) mm"
        @test isnothing(Ext._info_card(gui)) && !c.scene.visible[]
        # the kept marker moves with the card between the sidebar and the 3D view
        x = d.obj
        marker = only(x.plots)
        Ext._float!(gui, x)
        f = only(Ext._floating_cards(gui, x))
        @test !f.transient && isempty(gui.layout.inspector.pinned) && marker in gui.ax.scene.plots
        Ext._dock!(gui, x)
        @test Ext._is_docked(gui, x) && isempty(Ext._floating_cards(gui, x))
        @test marker in gui.ax.scene.plots
        # the dock button of a transient card keeps it, docked once
        Ext._show_inspection!(gui, _info(0.06))
        c = Ext._info_card(gui)
        notify(c.dock_button.clicks)
        @test length(gui.layout.inspector.pinned) == 2 && isnothing(Ext._info_card(gui))
        # unpinned: its marker is removed
        Ext._unpin!(gui, x)
        @test !Ext._is_pinned(gui, x) && !(marker in gui.ax.scene.plots)
        close(gui)
    end

    @testset "widget protocol of the blocks of Makie" begin
        # the cards use the public verbs only, see `card_input` and `card_show!`
        @test !isdefined(Ext, :_widget_input) && !isdefined(Ext, :_show!)
        fig = Figure()
        slider = Slider(fig[1, 1]; range = 0:10)
        toggle = Toggle(fig[2, 1])
        box = Textbox(fig[3, 1])
        button = Button(fig[4, 1])
        menu = Menu(fig[5, 1]; options = ["a", "b"])
        label = Label(fig[6, 1], "")
        @test BMO.card_input(slider) === slider.value
        @test BMO.card_input(toggle) === toggle.active
        @test BMO.card_input(box) === box.stored_string
        @test BMO.card_input(button) === button.clicks
        @test BMO.card_input(menu) === menu.selection
        # no input: labels and any other type
        @test isnothing(BMO.card_input(label))
        @test isnothing(BMO.card_input(42))
        BMO.card_show!(slider, 7)
        @test slider.value[] == 7
        BMO.card_show!(toggle, true)
        @test toggle.active[]
        BMO.card_show!(button, "hide")
        @test button.label[] == "hide"
        BMO.card_show!(label, 1.5)
        @test label.text[] == "1.5"
        # a textbox shows the value without an input, but keeps what is typed, unless `force`
        n = Ref(0)
        on(_ -> (n[] += 1), box.stored_string)
        BMO.card_show!(box, 12)
        @test box.displayed_string[] == "12" && box.stored_string[] == "12"
        box.focused[] = true
        BMO.card_show!(box, 13)
        @test box.displayed_string[] == "12"
        BMO.card_show!(box, 13; force = true)
        @test box.displayed_string[] == "13"
        @test n[] == 0
        # shows nothing for any other type
        @test isnothing(BMO.card_show!(menu, "b"))
        @test isnothing(BMO.card_show!(42, 1))
    end
end

end
