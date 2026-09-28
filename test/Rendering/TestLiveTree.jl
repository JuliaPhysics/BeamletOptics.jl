module TestLiveTree

using BeamletOptics
using Makie
using Test

@testset "Live object tree" begin
    Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)
    @test !isnothing(Ext)
    Row = Ext._TreeRow

    # A system with a group of objects, a source and a clip plane; keys are strings
    function _rows(n = 4)
        rows = [Row("sys", "System 1", 0, :system, true, true, nothing),
            Row("grp", "Mount", 1, :group, true, true, true)]
        for i in 1:n
            push!(rows, Row("obj$i", "Lens L$i", 2, :lens, false, false, isodd(i)))
        end
        push!(rows, Row("src", "Source 1", 0, :source, false, false, true))
        return rows
    end

    function _fixture(; size = (300, 400))
        fig = Figure(; size)
        tree = Ext._ObjectTree(fig[1, 1])
        return fig, tree
    end

    # Moves the mouse to the pixel `(x, y)` of the tree scene, `x` from the left, `y` from the top
    function _mouse!(tree, x, y)
        vp = tree.scene.viewport[]
        o, w = Makie.origin(vp), Makie.widths(vp)
        events(tree.scene).mouseposition[] = (o[1] + x, o[2] + w[2] - y)
        return nothing
    end

    # Clicks the part (:label, :eye or :expander) of the i-th row
    function _click!(tree, i, part)
        c = Ext._row_columns(tree, tree.rows[i])
        x = part === :label ? c.label + 5 : getfield(c, part)
        vp = tree.scene.viewport[]
        y = Makie.widths(vp)[2] - Ext._row_y(tree, i)
        _mouse!(tree, x, y)
        ev = events(tree.scene)
        ev.mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        ev.mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release)
        return nothing
    end

    _labels(tree) = tree.plots.labels.text[]

    @testset "rows are displayed" begin
        fig, tree = _fixture()
        @test isempty(_labels(tree))
        rows = _rows()
        Ext._set_rows!(tree, rows)
        @test _labels(tree) == [r.label for r in rows]
        # expanders for the system and the group, eyes for all but the system
        @test length(tree.plots.expanders[1][]) == 2
        @test length(tree.plots.eyes[1][]) == length(rows) - 1
        @test length(tree.plots.markers[1][]) == length(rows)
        # hidden rows are muted
        colors = tree.plots.labels.color[]
        @test colors[4] == tree.muted_color
        @test colors[3] == tree.text_color
        # indentation by depth
        x = [p[1] for p in tree.plots.labels[1][]]
        @test x[1] < x[2] < x[3]
        @test x[1] == x[end]
        # no scroll bar if the rows fit
        @test !tree.plots.scrollbar.visible[]
        # long labels are ellipsized to the width
        long = Row("long", "A very long name of an object " ^ 5, 0, :lens, false, false, true)
        Ext._set_rows!(tree, [long])
        label = only(_labels(tree))
        @test endswith(label, "…")
        @test length(label) < length(long.label)
        w = Makie.widths(tree.scene.viewport[])[1]
        @test Ext._row_columns(tree, long).label + Ext._label_width(tree, label) <= w
    end

    @testset "clicks" begin
        fig, tree = _fixture()
        Ext._set_rows!(tree, _rows())
        hits = Dict(:clicked => Any[], :eye => Any[], :expand => Any[])
        on(k -> push!(hits[:clicked], k), tree.clicked)
        on(k -> push!(hits[:eye], k), tree.eye_clicked)
        on(k -> push!(hits[:expand], k), tree.expand_clicked)
        _click!(tree, 3, :label)
        @test hits[:clicked] == ["obj1"]
        _click!(tree, 4, :eye)
        @test hits[:eye] == ["obj2"]
        _click!(tree, 2, :expander)
        @test hits[:expand] == ["grp"]
        # the type marker selects, too
        _click!(tree, 7, :marker)
        @test hits[:clicked] == ["obj1", "src"]
        # no eye and no expander: the columns select the row
        _click!(tree, 3, :expander)
        _click!(tree, 1, :eye)
        @test hits[:clicked] == ["obj1", "src", "obj1", "sys"]
        @test length(hits[:eye]) == 1 && length(hits[:expand]) == 1
        # below the last row: nothing
        _mouse!(tree, 50, 390)
        events(tree.scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        events(tree.scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release)
        @test length(hits[:clicked]) == 4
        # press on one row, release on another: nothing
        c = Ext._row_columns(tree, tree.rows[3])
        h = Makie.widths(tree.scene.viewport[])[2]
        _mouse!(tree, c.label + 5, h - Ext._row_y(tree, 3))
        events(tree.scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        _mouse!(tree, c.label + 5, h - Ext._row_y(tree, 5))
        events(tree.scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release)
        @test length(hits[:clicked]) == 4
        # presses outside of the tree are not consumed
        other = Ref(0)
        on(events(tree.scene).mousebutton, priority = 0) do _
            other[] += 1
            return Consume(false)
        end
        _click!(tree, 3, :label)
        @test other[] == 0
        events(tree.scene).mouseposition[] = (1000.0, 1000.0)
        events(tree.scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press)
        events(tree.scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release)
        @test other[] == 2
    end

    @testset "selection and scrolling" begin
        fig, tree = _fixture()
        rows = _rows(100)
        Ext._set_rows!(tree, rows)
        h = Makie.widths(tree.scene.viewport[])[2]
        @test tree.plots.scrollbar.visible[]
        # only the rows in view are drawn
        nview = length(_labels(tree))
        @test nview < 25
        @test first(_labels(tree)) == "System 1"
        @test !tree.plots.selection.visible[]
        # select a row far down: it is scrolled into view and highlighted
        Ext._set_selected!(tree, "obj80")
        @test tree.plots.selection.visible[]
        y = Ext._row_y(tree, 82)
        @test tree.row_height / 2 <= y <= h - tree.row_height / 2
        @test "Lens L80" in _labels(tree)
        rect = first(tree.plots.selection[1][])
        @test Makie.origin(rect)[2] ≈ y - tree.row_height / 2
        # the selection survives new rows, and is cleared by nothing
        Ext._set_rows!(tree, rows)
        @test tree.plots.selection.visible[]
        Ext._set_selected!(tree, nothing)
        @test !tree.plots.selection.visible[]
        # a key that is not shown highlights nothing
        Ext._set_selected!(tree, "missing")
        @test !tree.plots.selection.visible[]

        # the wheel scrolls while the mouse is over the tree and clamps at the ends
        ev = events(tree.scene)
        other = Ref(0)
        on(ev.scroll, priority = 0) do _
            other[] += 1
            return Consume(false)
        end
        _mouse!(tree, 50, 50)
        for _ in 1:100
            ev.scroll[] = (0.0, 1.0)
        end
        @test tree.offset == 0
        @test first(_labels(tree)) == "System 1"
        ev.scroll[] = (0.0, -1.0)
        @test tree.offset == 3 * tree.row_height
        for _ in 1:100
            ev.scroll[] = (0.0, -1.0)
        end
        @test tree.offset == Ext._max_offset(tree)
        @test last(_labels(tree)) == "Source 1"
        # the last row ends above the bottom padding
        @test Ext._row_y(tree, length(rows)) ≈ Ext._TREE_PAD + tree.row_height / 2
        @test other[] == 0
        # outside of the tree, the wheel is not consumed and does not scroll
        offset = tree.offset
        events(tree.scene).mouseposition[] = (1000.0, 1000.0)
        ev.scroll[] = (0.0, 1.0)
        @test other[] == 1
        @test tree.offset == offset
        # fewer rows clamp the offset
        Ext._set_rows!(tree, _rows())
        @test tree.offset == 0
        @test !tree.plots.scrollbar.visible[]
    end

    @testset "constant number of plots" begin
        _, small = _fixture()
        Ext._set_rows!(small, _rows(10))
        _, large = _fixture()
        Ext._set_rows!(large, _rows(5000))
        @test length(small.scene.plots) == length(large.scene.plots)
        h = Makie.widths(large.scene.viewport[])[2]
        @test length(_labels(large)) <= ceil(Int, h / large.row_height) + 1
        Ext._set_selected!(large, "obj4000")
        @test length(small.scene.plots) == length(large.scene.plots)
        @test "Lens L4000" in _labels(large)
        @test length(_labels(large)) <= 20
    end

    @testset "follows the layout" begin
        fig = Figure(; size = (400, 300))
        tree = Ext._ObjectTree(fig[1, 1])
        Box(fig[1, 2])
        colsize!(fig.layout, 1, Fixed(150))
        Ext._set_rows!(tree, _rows(40))
        vp = tree.scene.viewport[]
        @test Makie.widths(vp)[1] == 150
        n = length(_labels(tree))
        resize!(fig.scene, 400, 600)
        vp = tree.scene.viewport[]
        @test Makie.widths(vp)[2] > 500
        @test length(_labels(tree)) > n
        colsize!(fig.layout, 1, Fixed(200))
        @test Makie.widths(tree.scene.viewport[])[1] == 200
        # requested width
        fig = Figure(; size = (400, 300))
        tree = Ext._ObjectTree(fig[1, 1]; width = 120)
        Box(fig[1, 2])
        @test Makie.widths(tree.scene.viewport[])[1] == 120
        close(tree)
        @test isempty(tree.listeners)
    end

    @testset "timings" begin
        _, tree = _fixture()
        rows = _rows(5000)
        Ext._set_rows!(tree, rows)
        Ext._set_selected!(tree, "obj2500")
        _mouse!(tree, 50, 50)
        t_rows = @elapsed for _ in 1:20
            Ext._set_rows!(tree, rows)
        end
        t_select = @elapsed for i in 1:20
            Ext._set_selected!(tree, "obj$(250 * i)")
        end
        t_scroll = @elapsed for i in 1:100
            events(tree.scene).scroll[] = (0.0, isodd(i ÷ 10) ? 1.0 : -1.0)
        end
        @info "Object tree with 5000 rows" set_rows_ms = 1000 * t_rows / 20 select_ms = 1000 *
            t_select / 20 scroll_ms = 1000 * t_scroll / 100
        # generous limits for slow runners, typical values are well below 1 ms
        @test t_rows / 20 < 0.05
        @test t_scroll / 100 < 0.05
    end
end

end
