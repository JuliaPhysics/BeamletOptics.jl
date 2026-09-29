module TestLiveIcons

using BeamletOptics
using Makie
using Test

@testset "Live view icons" begin
    Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)
    @test !isnothing(Ext)

    @testset "icon set" begin
        names = (:trace, :auto_trace, :home, :fit, :views, :save_view, :orthographic, :clip,
            :clip_beams, :sources, :measure, :export, :panel_left, :panel_right, :panel_bottom,
            :help, :eye, :eye_off, :expand, :collapse, :lens, :mirror, :detector, :source, :group,
            :clip_plane, :mesh, :object, :system, :beamsplitter, :polarizer, :pin, :pinned, :chart,
            :float, :dock)
        @test Set(keys(Ext._ICONS)) == Set(names)
        for name in names
            icon = Ext._icon(name)
            @test icon isa Makie.BezierPath
            @test count(c -> c isa Makie.ClosePath, icon.commands) >= 1
            bb = Makie.bbox(icon)
            @test all(minimum(bb) .>= -0.5) && all(maximum(bb) .<= 0.5)
            # Not a degenerate path
            @test all(Makie.widths(bb) .> 0.2)
            # Cached: parsed once
            @test Ext._icon(name) === icon
        end
        @test_throws ArgumentError Ext._icon(:no_such_icon)
        # an own icon is used as it is, e.g. on an icon button
        path = Makie.BezierPath("M -0.3 -0.3 L 0.3 -0.3 L 0 0.3 Z")
        @test Ext._icon(path) === path
        fig = Figure()
        @test Ext._IconButton(fig[1, 1]; icon = path).icon[] === path
        t = Ext._IconToggle(fig[1, 2]; icon = path, icon_off = :eye_off)
        @test t.icon[] === Ext._icon(:eye_off)
        t.active[] = true
        @test t.icon[] === path
    end

    @testset "SVG path parser" begin
        # Absolute and relative commands give the same path; the view box maps to the unit square
        a = Ext._svg_path("M0-960H960V0H0Z")
        b = Ext._svg_path("m0-960h960v960h-960z")
        @test a == b
        @test Makie.bbox(a) == Rect2d(-0.5, -0.5, 1, 1)
        # y down in SVG, y up in the path
        p = Ext._svg_path("M0-960L480-480L0 0Z")
        @test p.commands[1].p ≈ Point2d(-0.5, 0.5)
        # Implicit lines after M, degenerate subpaths are dropped
        @test Ext._svg_path("M0-960 960-960 960 0Zm10-10Z") ==
              Ext._svg_path("M0-960L960-960L960 0Z")
        # Smooth quadratic curve: the control point is reflected
        q = Ext._svg_path("M0-480Q240-960 480-480T960-480")
        @test q.commands[3].c1[2] < 0 # second arc bends downwards
        @test_throws ArgumentError Ext._svg_path("M0 0A10 10 0 0 1 20 20")
        # A subpath closed by Z ends with an explicit line to its start, such that Makie's
        # bounding box, which scales the marker, includes the start (e.g. the tip of an arrow)
        t = Ext._svg_path("M0-960L480-960L480-480ZM960-480L720-240L720-720Z")
        @test t.commands[end - 1] == Makie.LineTo(Point2d(0.5, 0.0))
        @test maximum(Makie.bbox(t))[1] ≈ 0.5
    end

    fig = Figure(size = (400, 200))
    b = Ext._IconButton(fig[1, 1]; icon = :trace, tooltip = "Trace (t)", tooltip_delay = 0)
    t = Ext._IconToggle(fig[1, 2]; icon = :eye, icon_off = :eye_off, tooltip = "Visible")
    Box(fig[2, 1:3])
    Makie.update_state_before_display!(fig)
    scene = fig.scene
    center(w) = (bb = w.box.layoutobservables.computedbbox[];
        Tuple(Makie.origin(bb) .+ Makie.widths(bb) ./ 2))
    function click!(w)
        events(scene).mouseposition[] = center(w)
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release)
        return nothing
    end

    @testset "construction" begin
        @test b.clicks[] == 0 && !b.hovered[]
        @test length(b.plots) == 3 && all(p -> p in b.box.blockscene.plots, b.plots)
        @test b.icon[] === Ext._icon(:trace)
        @test b.background[] == Ext._TRANSPARENT
        @test b.icon_color[] == Ext._ICON_COLOR
        # The tooltip is hidden and drawn on top
        @test !b.plots[3].visible[]
        @test Makie.transformationmatrix(b.plots[3])[][3, 4] == Ext._TOOLTIP_Z
        @test b.tooltip[] == "Trace (t)"
        @test (w = b.box.layoutobservables.computedbbox[].widths; w[1] == w[2] == 28)
    end

    @testset "hover" begin
        events(scene).mouseposition[] = (1, 1)
        @test !b.hovered[]
        n = Ref(0)
        on(_ -> n[] += 1, b.background)
        events(scene).mouseposition[] = center(b)
        @test b.hovered[] && b.background[] == Ext._ICON_HOVER_COLOR
        @test b.plots[3].visible[] # no delay
        # Below the button, extending to the right at the left edge of the window
        @test b.plots[3].placement[] === :below && b.plots[3].align[] == 0.15f0
        @test b.plots[3][1][] ≈ Point2f(center(b)[1], b.box.layoutobservables.computedbbox[].origin[2])
        # Moving within the button changes nothing
        events(scene).mouseposition[] = center(b) .+ 3
        @test n[] == 1
        events(scene).mouseposition[] = center(t)
        @test !b.hovered[] && t.hovered[]
        @test b.background[] == Ext._TRANSPARENT && !b.plots[3].visible[]
        @test n[] == 2
        # The toggle waits for the default delay
        @test !t.plots[3].visible[]
        events(scene).mouseposition[] = (1, 1)
        @test !t.hovered[]
    end

    @testset "button clicks" begin
        consumed = Ref(0)
        # Clicks on the button do not reach listeners of lower priority
        on(_ -> (consumed[] += 1; Consume(false)), events(scene).mousebutton; priority = -1)
        click!(b)
        @test b.clicks[] == 1
        click!(b)
        @test b.clicks[] == 2
        @test consumed[] == 2 # the press events, the release events are consumed
        # A click elsewhere does not count
        events(scene).mouseposition[] = (1, 1)
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release)
        @test b.clicks[] == 2
    end

    @testset "toggle" begin
        @test !t.active[]
        @test t.icon[] === Ext._icon(:eye_off) && t.icon_color[] == Ext._ICON_COLOR
        click!(t)
        @test t.active[]
        @test t.icon[] === Ext._icon(:eye)
        @test t.icon_color[] == Ext._ICON_ACTIVE_ICON_COLOR
        @test t.background[] == Ext._ICON_ACTIVE_COLOR
        click!(t)
        @test !t.active[]
        # Set from code
        events(scene).mouseposition[] = (1, 1)
        t.active[] = true
        @test t.background[] == Ext._ICON_ACTIVE_COLOR && t.icon[] === Ext._icon(:eye)
        t.active[] = false
        @test t.background[] == Ext._TRANSPARENT && t.icon[] === Ext._icon(:eye_off)
        # A given observable is used, colors are taken from the keyword arguments
        active = Observable(true)
        t2 = Ext._IconToggle(fig[1, 3]; icon = :clip, active, active_color = :red,
            active_icon_color = :white, size = 36)
        @test t2.active === active
        @test t2.background[] == RGBAf(1, 0, 0, 1) && t2.icon_color[] == RGBAf(1, 1, 1, 1)
        active[] = false
        @test t2.background[] == Ext._TRANSPARENT
        Makie.update_state_before_display!(fig)
        @test t2.box.layoutobservables.computedbbox[].widths[1] == 36
        # delete! removes the plots
        delete!(t2)
        @test !(t2.box.blockscene in fig.scene.children)
    end
end

end
