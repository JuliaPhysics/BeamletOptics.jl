module TestLiveCustom

using BeamletOptics
using Makie
using Test

const BMO = BeamletOptics

@testset "Live view customization" begin
    Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)
    @test !isnothing(Ext)
    GLB = Makie.GridLayoutBase

    # Rays along +y, the mirror at 45° reflects them along +x onto the detector. The source is a
    # beam group, which is solved as a preview while moving, see `preview`
    function _fixture(; kwargs...)
        m = RoundPlanoMirror(25e-3, 5e-3)
        zrotate3d!(m, deg2rad(45))
        translate3d!(m, [0, 0.1, 0])
        pd = Detector(5e-3)
        zrotate3d!(pd, -π / 2)
        translate3d!(pd, [0.1, 0.1, 0])
        cs = CollimatedSource([0.0, 0, 0], [0.0, 1, 0], 2e-3, 1e-6; num_rings = 2, num_rays = 40)
        gui = live_view(System([m, pd]) => cs; trace_budget = Inf, throttle = false,
            fine_step = 1e-4, idle_delay = 10.0, kwargs...)
        return gui, m, pd
    end

    _key!(gui, key) = (events(gui.ax.scene).keyboardbutton[] = Makie.KeyEvent(key, Keyboard.press))
    _tick!(gui) = (events(gui.ax.scene).tick[] = Makie.Tick(Makie.RegularRenderTick, 0, 0.0, 1.0))

    # A preview solve after a key step of the mirror, then the full solve once the movement pauses
    function _move!(gui, m)
        gui.controls.selected[] = m
        _key!(gui, Keyboard.up)
        @test gui.trace.preview
        return nothing
    end
    function _pause!(gui)
        gui.trace.last_change -= 100
        _tick!(gui)
        @test !gui.trace.preview
        return nothing
    end

    # A panel with an axis whose `update` counts its calls
    function _counted_panel!(gui, title; kwargs...)
        n = Ref(0)
        layout = add_panel!(gui, title; kwargs...) do layout
            Axis(layout[1, 1])
            return g -> (n[] += 1)
        end
        return layout, n
    end

    @testset "panels: compact" begin
        gui, m, _ = _fixture()
        layout, n = _counted_panel!(gui, "Power")
        # below the detector panel, in the grid in fig[1, 2], with its title above
        box = GLB.gridcontent(layout).parent
        @test GLB.gridcontent(box).parent === gui.layout.panels
        @test GLB.gridcontent(box).span.rows == 2:2
        @test any(c -> c.content isa Label && c.content.text[] == "Power", box.content)
        @test only(gui.custom.panels).layout === layout
        # updated right away, after full solves only
        @test n[] == 1
        _move!(gui, m)
        @test n[] == 1
        _pause!(gui)
        @test n[] == 2
        _key!(gui, Keyboard.t)
        @test n[] == 3
        # an own axis below the detector panel still works
        @test Axis(gui.fig[1, 2][3, 1]) isa Axis
        # a `do` block that returns a plot has no update
        layout2 = add_panel!(gui, "Plot") do layout
            lines!(Axis(layout[1, 1]), [1, 2], [1, 2])
        end
        @test isnothing(gui.custom.panels[2].update)
        @test GLB.gridcontent(GLB.gridcontent(layout2).parent).span.rows == 4:4
        # errors in `update` are logged once
        @test_logs (:error, r"panel \"Broken\"") add_panel!(_ -> (_ -> error("broken panel")), gui, "Broken")
        @test gui.custom.panels[3].last_error == "broken panel"
        @test_logs _key!(gui, Keyboard.t)
        close(gui)

        # without detector panels: a new column, the rows below the 3D view span both columns
        gui, _, _ = _fixture(; detectors = [])
        @test isnothing(gui.layout.panels)
        layout, n = _counted_panel!(gui, "Power")
        grid = gui.layout.panels
        @test GLB.gridcontent(grid).span.cols == 2:2
        @test GLB.gridcontent(GLB.gridcontent(layout).parent).span.rows == 1:1
        @test GLB.gridcontent(gui.layout.status_row).span.cols == 1:2
        @test GLB.gridcontent(gui.layout.tool_row).span.cols == 1:2
        @test n[] == 1
        close(gui)
    end

    @testset "panels: app" begin
        gui, m, _ = _fixture(; layout = :app)
        tabs = gui.layout.tabs
        layout, n = _counted_panel!(gui, "Power")
        # a tab behind the detector tab, which stays active; hidden, i.e. not updated
        @test tabs.bar.titles == ["Detector 1", "Power"]
        @test last(gui.layout.dock_panels) == ("Power" => layout)
        @test tabs.active == 1
        @test tabs.panels[2] === only(gui.custom.panels)
        @test n[] == 0
        @test tabs.panels[2] in tabs.stale
        _key!(gui, Keyboard.t)
        @test n[] == 0
        # updated once its tab is opened
        Ext._activate_tab!(gui, 2)
        @test n[] == 1
        @test !(tabs.panels[2] in tabs.stale)
        # shown: after full solves only
        _move!(gui, m)
        @test n[] == 1
        _pause!(gui)
        @test n[] == 2
        # hidden again: stale after a solve, updated when the tab is opened
        Ext._activate_tab!(gui, 1)
        _key!(gui, Keyboard.t)
        @test n[] == 2
        Ext._activate_tab!(gui, 2)
        @test n[] == 3
        # a collapsed dock updates nothing, the active panel once it is expanded
        gui.layout.collapse.dock.active[] = false
        _key!(gui, Keyboard.t)
        @test n[] == 3
        gui.layout.collapse.dock.active[] = true
        @test n[] == 4
        # `select = true` shows the new tab at once
        layout2, n2 = _counted_panel!(gui, "Spectrum"; select = true)
        @test tabs.active == 3
        @test n2[] == 1
        @test layout2.parent !== nothing
        close(gui)

        # without detector panels, the first panel is the active tab of the shown dock
        gui, _, _ = _fixture(; layout = :app, detectors = [])
        @test !gui.layout.dock.shown
        layout, n = _counted_panel!(gui, "Power")
        @test gui.layout.dock.shown
        @test gui.layout.tabs.active == 1
        @test n[] == 1
        close(gui)
    end

    @testset "controls and retrace! ($layout)" for layout in (:compact, :app)
        gui, m, pd = _fixture(; layout, auto_trace = false)
        clicked = Ref(0)
        tb = Ref{Any}(nothing)
        content = add_controls!(gui, "Mine") do l
            tb[] = Textbox(l[1, 1]; placeholder = "x")
            b = Button(l[1, 2]; label = "go")
            on(_ -> (clicked[] += 1), b.clicks)
        end
        @test content isa GridLayout
        if layout == :compact
            # a row above the status row, over the width of the window
            box = GLB.gridcontent(content).parent
            @test GLB.gridcontent(box).span.rows.stop + 1 ==
                  GLB.gridcontent(gui.layout.status_row).span.rows.start
            @test GLB.gridcontent(box).span.cols == 1:2
        else
            # a section of the left sidebar
            @test last(gui.layout.sections[:left]) == ("Mine" => content)
        end
        # a focused textbox of the controls takes the keyboard
        @test gui.custom.boxes == [tb[]]
        cam = cameracontrols(gui.ax.scene)
        @test !Ext._typing(gui)
        tb[].focused[] = true
        @test Ext._typing(gui)
        @test gui.controls.ignore_keys()
        @test !cam.selected[]
        pos = position(m)
        gui.controls.selected[] = m
        _key!(gui, Keyboard.up)
        @test position(m) == pos
        tb[].focused[] = false
        @test !Ext._typing(gui)
        @test cam.selected[]
        # retrace! with auto tracing off: the change is made, the beams are outdated
        hits0 = length(BMO.hits(pd))
        @test hits0 == 40
        retrace!(() -> zrotate3d!(m, deg2rad(20)), gui)
        @test gui.trace.stale
        @test length(BMO.hits(pd)) == hits0
        # with auto tracing on, the change is solved: the reflected beam misses the detector
        gui.trace.auto[] = true
        @test !gui.trace.stale
        @test isnothing(BMO.hits(pd))
        zrotate3d!(m, deg2rad(-20))
        retrace!(gui)
        _pause!(gui)
        @test length(BMO.hits(pd)) == 40
        close(gui)
    end

    @testset "update of controls ($layout)" for layout in (:compact, :app)
        gui, m, _ = _fixture(; layout)
        n = Ref(0)
        add_controls!(gui, "Counted") do l
            Label(l[1, 1], "x")
            return g -> (n[] += 1)
        end
        @test only(gui.custom.controls).title == "Counted"
        # called right away, then after full solves only, like the update of a panel
        @test n[] == 1
        _move!(gui, m)
        @test n[] == 1
        _pause!(gui)
        @test n[] == 2
        _key!(gui, Keyboard.t)
        @test n[] == 3
        # a `do` block that ends with a listener has no update
        add_controls!(gui, "Listener") do l
            b = Button(l[1, 1]; label = "go")
            on(_ -> nothing, b.clicks)
        end
        @test isnothing(gui.custom.controls[2].update)
        # errors in `update` are logged once
        @test_logs (:error, r"controls \"Broken\"") add_controls!(_ -> (_ -> error("broken controls")), gui, "Broken")
        @test gui.custom.controls[3].last_error == "broken controls"
        @test_logs _key!(gui, Keyboard.t)
        @test n[] == 4
        close(gui)
    end

    @testset "tools ($layout)" for layout in (:compact, :app)
        gui, _, _ = _fixture(; layout)
        calls = Any[]
        b = add_tool!(g -> push!(calls, g), gui, "Tool"; key = Keyboard._2, tooltip = "My tool")
        t = add_tool!((g, a) -> push!(calls, a), gui, "Toggle"; toggle = true, icon = :measure)
        # a click calls `f(gui)`, a toggle passes its state
        b.clicks[] += 1
        @test calls == [gui]
        t.active[] = true
        t.active[] = false
        @test calls == [gui, true, false]
        # the key clicks the button, not while typing
        _key!(gui, Keyboard._2)
        @test length(calls) == 4
        tb = Ref{Any}(nothing)
        add_controls!(l -> (tb[] = Textbox(l[1, 1])), gui, "Box")
        tb[].focused[] = true
        _key!(gui, Keyboard._2)
        @test length(calls) == 4
        tb[].focused[] = false
        # the key of a toggle switches it
        add_tool!((g, a) -> push!(calls, a), gui, "Keyed"; toggle = true, key = Keyboard._3)
        _key!(gui, Keyboard._3)
        @test calls[end] === true
        # taken keys and unknown icons
        err(f) = try
            f()
            ""
        catch e
            e isa ArgumentError ? e.msg : rethrow()
        end
        @test occursin("trace (t)", err(() -> add_tool!(identity, gui, "x"; key = Keyboard.t)))
        @test occursin("Camera3D (forward_key)", err(() -> add_tool!(identity, gui, "x"; key = Keyboard.w)))
        @test occursin("Camera3D (fix_x_key)", err(() -> add_tool!(identity, gui, "x"; key = Keyboard.x)))
        @test occursin("source markers", err(() -> add_tool!(identity, gui, "x"; key = Keyboard._1)))
        @test occursin("key steps", err(() -> add_tool!(identity, gui, "x"; key = Keyboard.left)))
        @test occursin("the tool \"Tool\"", err(() -> add_tool!(identity, gui, "x"; key = Keyboard._2)))
        @test occursin("measure", err(() -> add_tool!(identity, gui, "x"; icon = :nonsense)))
        # an own icon
        path = BezierPath("M -0.3 -0.3 L 0.3 -0.3 L 0 0.3 Z")
        own = add_tool!(identity, gui, "Own"; icon = path)
        own_toggle = add_tool!((g, a) -> nothing, gui, "Own toggle"; icon = path, toggle = true)
        # errors of `f` are logged
        bad = add_tool!(_ -> error("broken tool"), gui, "Bad")
        @test_logs (:error, r"tool \"Bad\"") (bad.clicks[] += 1)
        if layout == :compact
            # text widgets in the tool row, which ends with its filler
            @test b isa Button
            @test b.label[] == "Tool (2)"
            @test t isa Toggle
            row = gui.layout.tool_row
            filler = only(GLB.contents(row[1, GLB.ncols(row)]))
            @test filler isa Label && filler.text[] == ""
            # the compact layout shows the name, not the icon
            @test own isa Button && own.label[] == "Own"
        else
            # icon buttons in the toolbar group `:user` before "Help"
            @test b isa Ext._IconButton
            @test b.tooltip[] == "My tool (2)"
            @test t isa Ext._IconToggle
            @test first.(gui.layout.groups)[(end - 1):end] == [:user, :help]
            # the own icon in the toolbar, also on a toggle, on and off
            @test own.icon[] === path
            @test own_toggle.icon[] === path
            own_toggle.active[] = true
            @test own_toggle.icon[] === path
        end
        close(gui)
    end
end

end
