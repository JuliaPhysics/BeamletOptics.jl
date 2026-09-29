module TestLiveDock

using BeamletOptics
using Makie
using Test

const BMO = BeamletOptics

@testset "Live view analysis dock" begin
    Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)
    @test !isnothing(Ext)

    # Rays along +y through two transparent detectors onto a third one; the first panel records
    # its history, the Gaussian beam of a second system hits the intensity panel PD4
    function _fixture(; layout = :app, kwargs...)
        pds = [Detector(25e-3, false), Detector(25e-3, false), Detector(25e-3)]
        foreach(((k, pd),) -> translate3d!(pd, [0, 0.05k, 0]), enumerate(pds))
        pd4 = Detector(10e-3)
        translate3d!(pd4, [0.06, 0.1, 0])
        src = CollimatedSource([0.0, 0, 0], [0.0, 1, 0], 4e-3, 1e-6; num_rings = 2, num_rays = 40)
        g = GaussianBeamlet([0.06, 0, 0], [0.0, 1, 0], 1e-6, 0.5e-3)
        labels = Dict(pds[1] => "PD1", pds[2] => "PD2", pds[3] => "PD3", pd4 => "PD4")
        detectors = [pds[1] => (:spot, (; history = true)), pds[2], pds[3],
            pd4 => (:intensity, (; n = 20, profiles = true))]
        gui = live_view(System(pds) => src, System([pd4]) => g; layout, labels, detectors,
            trace_budget = Inf, kwargs...)
        return gui
    end

    # Counts the updates of the plots of each panel, i.e. its computations
    function _counters(gui)
        counts = zeros(Int, length(gui.panels))
        for (k, p) in enumerate(gui.panels)
            on(_ -> (counts[k] += 1), p.xy)
            on(_ -> (counts[k] += 1), p.heat_I)
        end
        return counts
    end

    _tabs(gui) = gui.layout.tabs
    _shown(p) = p.ax.blockscene.visible[]

    # A click on the tab `i` of the tab bar
    function _click_tab!(gui, i)
        bar = _tabs(gui).bar
        s = only(s for s in Ext._tab_slots(bar).slots if s.i == i)
        o = Makie.origin(bar.scene.viewport[])
        events(gui.fig.scene).mouseposition[] = (o[1] + (s.x0 + s.x1) / 2, o[2] + Ext._TAB_HEIGHT / 2)
        for action in (Mouse.press, Mouse.release)
            events(gui.fig.scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, action)
        end
        return nothing
    end

    @testset "tabs" begin
        gui = _fixture()
        tabs = _tabs(gui)
        @test tabs.bar isa Ext._TabBar
        @test tabs.bar.titles == ["PD1", "PD2", "PD3", "PD4"]
        @test first.(gui.layout.dock_panels) == tabs.bar.titles
        @test tabs.panels == gui.panels
        @test tabs.active == 1 && tabs.bar.active == 1
        # only the active panel is laid out and drawn
        @test _shown(gui.panels[1])
        @test !any(_shown, gui.panels[2:end])
        @test maximum(gui.panels[2].ax.layoutobservables.computedbbox[])[1] < 0
        # a constant number of plots
        @test length(tabs.bar.scene.plots) == 3
        # a click on a tab switches the panel
        _click_tab!(gui, 3)
        @test tabs.active == 3 && tabs.bar.active == 3
        @test _shown(gui.panels[3])
        @test !_shown(gui.panels[1])
        @test minimum(gui.panels[3].ax.layoutobservables.computedbbox[])[1] > 0
        # the history and the profiles beside the axis
        p1, p4 = gui.panels[1], gui.panels[4]
        @test Ext._GLB.gridcontent(p1.history_axes[1]).span.cols == 2:2
        @test Ext._GLB.gridcontent(p4.profiles_ax).span.cols == 2:2
        close(gui)
    end

    @testset "lazy panels" begin
        gui = _fixture()
        tabs = _tabs(gui)
        # initially only the active panel (and the history panel, which is the active one) is
        # computed, the other panels are stale
        @test occursin("PD1: 40 rays", gui.panels[1].ax.title[])
        @test gui.panels[2].ax.title[] == "PD2: no hits"
        @test Set(tabs.stale) == Set(gui.panels[2:4])
        n = _counters(gui)
        # after a solve only the active panel is computed
        Ext._resolve!(gui, nothing)
        @test n == [1, 0, 0, 0]
        @test Ext._computed_panels(gui, false) == gui.panels[1:1]
        # a stale tab is computed once when it becomes active
        Ext._activate_tab!(gui, 2)
        @test n == [1, 1, 0, 0]
        @test occursin("PD2: 40 rays", gui.panels[2].ax.title[])
        @test !(gui.panels[2] in tabs.stale)
        Ext._activate_tab!(gui, 1)
        Ext._activate_tab!(gui, 2)
        @test n == [1, 1, 0, 0]
        # the intensity panel in a job, like the refinement of a coarse preview
        Ext._activate_tab!(gui, 4)
        @test n[4] == 1
        @test gui.panels[4].heat_plot.visible[]
        @test !isempty(gui.panels[4].profile_x[])
        @test isnothing(gui.trace.job)
        @test !gui.trace.stale
        # the history panel records every full solve, also while it is not shown, but draws it
        # only when shown
        p1 = gui.panels[1]
        drawn = Ref(0)
        on(_ -> (drawn[] += 1), p1.history_value)
        k = p1.history_count
        Ext._resolve!(gui, nothing)
        Ext._resolve!(gui, nothing)
        @test p1.history_count == k + 2
        @test length(p1.history_value[]) == k + 2
        @test drawn[] == 0
        @test n == [1, 1, 0, 3]
        # its field is kept: showing it needs no second computation
        @test haskey(tabs.fields, p1)
        Ext._activate_tab!(gui, 1)
        @test n == [2, 1, 0, 3]
        @test drawn[] == 1
        @test p1.history_count == k + 2
        # previews are not recorded, the history panel is not computed
        @test Ext._computed_panels(gui, true) == [p1]
        Ext._activate_tab!(gui, 2)
        @test Ext._computed_panels(gui, true) == [gui.panels[2]]
        close(gui)
    end

    @testset "collapsed dock" begin
        gui = _fixture()
        tabs = _tabs(gui)
        Ext._activate_tab!(gui, 2)
        n = _counters(gui)
        gui.layout.collapse.dock.active[] = false
        @test !gui.layout.dock.shown
        @test isempty(Ext._shown_panels(gui))
        # only the history panel is computed, for its metrics
        Ext._resolve!(gui, nothing)
        @test n == [0, 0, 0, 0]
        @test gui.panels[2] in tabs.stale
        @test Ext._computed_panels(gui, false) == gui.panels[1:1]
        # expanding computes the active panel
        gui.layout.collapse.dock.active[] = true
        @test n == [0, 1, 0, 0]
        @test !(gui.panels[2] in tabs.stale)
        @test _shown(gui.panels[2])
        @test !_shown(gui.panels[1])
        close(gui)
    end

    @testset "stale hits" begin
        gui = _fixture()
        tabs = _tabs(gui)
        n = _counters(gui)
        # a solve that did not complete (e.g. cancelled) leaves incomplete hits: stale tabs are
        # not computed until the next solve
        Ext._on_solve_started!(gui)
        Ext._activate_tab!(gui, 3)
        @test n == [0, 0, 0, 0]
        @test gui.panels[3] in tabs.stale
        Ext._resolve!(gui, nothing)
        @test n == [0, 0, 1, 0]
        # the app layout computes after each solve
        @test !(gui.panels[3] in tabs.stale)
        close(gui)
    end

    @testset "coarse and preview panels" begin
        gui = _fixture()
        tabs = _tabs(gui)
        Ext._activate_tab!(gui, 4)
        # a panel that shows a coarse preview is computed again when it is shown next
        Ext._resolve!(gui, nothing; coarse = true)
        @test gui.trace.coarse
        @test size(gui.panels[4].heat_I[]) == (16, 16)
        @test Ext._shown_panels(gui) == gui.panels[4:4]
        Ext._activate_tab!(gui, 1)
        @test gui.panels[4] in tabs.stale
        close(gui)
    end

    @testset "slots" begin
        gui = _fixture()
        tabs = _tabs(gui)
        n = _counters(gui)
        # a dock panel of the user is a tab, which becomes active
        d = Ext._add_dock_panel!(gui, "Mine")
        ax = Axis(d[1, 1])
        @test tabs.bar.titles[end] == "Mine"
        @test tabs.active == 5
        @test isnothing(Ext._active_panel(tabs))
        @test !any(_shown, gui.panels)
        # no detector panel is shown, only the history panel is computed
        Ext._resolve!(gui, nothing)
        @test n == [0, 0, 0, 0]
        @test Ext._computed_panels(gui, false) == gui.panels[1:1]
        Ext._activate_tab!(gui, 2)
        @test !ax.blockscene.visible[]
        @test n == [0, 1, 0, 0]
        close(gui)
    end

    @testset "overflow" begin
        gui = _fixture(; size = (700, 700), view_cube = false)
        bar = _tabs(gui).bar
        for k in 1:8
            Ext._add_dock_panel!(gui, "Panel with a long name $k")
        end
        (; slots, overflow) = Ext._tab_slots(bar)
        @test overflow
        # the active (last) tab is scrolled into view
        @test slots[end].i == 12
        @test bar.first > 1
        # the wheel over the bar scrolls by whole tabs
        f = bar.first
        vp = bar.scene.viewport[]
        events(gui.fig.scene).mouseposition[] = Tuple(Makie.origin(vp) .+ Makie.widths(vp) ./ 2)
        events(gui.fig.scene).scroll[] = (0.0, 1.0)
        @test bar.first == f - 1
        Ext._activate_tab!(gui, 1)
        @test bar.first == 1
        close(gui)
    end

    @testset "compact layout computes all panels" begin
        gui = _fixture(; layout = :compact)
        n = _counters(gui)
        Ext._resolve!(gui, nothing)
        @test n == [1, 1, 1, 1]
        @test Ext._computed_panels(gui, false) == gui.panels
        @test Ext._shown_panels(gui) == gui.panels
        @test gui.panels[1].history_count == 2
        close(gui)
    end
end

end
