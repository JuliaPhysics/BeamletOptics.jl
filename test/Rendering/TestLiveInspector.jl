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
        # a lens: header, pose and properties, no type-dependent sections
        ctrl.selected[] = o.l1
        @test insp.name.text[] == "Lens 1"
        @test insp.type.text[] == "Lens"
        @test insp.icon[] === Ext._icon(:lens)
        @test gui.pose_boxes[2].displayed_string[] == "40.0"
        @test _value(gui, "Thickness") == "4 mm"
        @test _value(gui, "n(λ₀)") == "1.5"
        @test isnothing(_value(gui, "Position"))   # in the pose boxes
        @test isnothing(_value(gui, "Type"))       # in the header
        @test isempty(insp.sections)
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
        # the detector: hits and the panel section
        ctrl.selected[] = o.pd
        @test insp.name.text[] == "PD1"
        @test _value(gui, "Hits") == "1"
        @test _value(gui, "Size") == "(5, 5) mm"
        @test first.(insp.sections) == [BMO.Detector]
        # deselected: the summary again, the sections are removed
        ctrl.selected[] = nothing
        @test insp.name.text[] == "No selection"
        @test isempty(insp.sections)
        @test isempty(Ext._blocks!(Any[], insp.context))
        close(gui)
    end

    @testset "updates in place" begin
        gui, o = _fixture()
        insp, ctrl = gui.layout.inspector, gui.controls
        plots = (insp.list.labels, insp.list.values, insp.list.lines)
        n = length(gui.fig.scene.children)
        ctrl.selected[] = o.pd
        blocks = Ext._blocks!(Any[], insp.context)
        @test !isempty(blocks)
        # a move updates the values, the sections and plots stay
        gui.pose_boxes[3].stored_string[] = "10"
        @test _value(gui, "Hits") == "0"
        @test Ext._blocks!(Any[], insp.context) == blocks
        @test (insp.list.labels, insp.list.values, insp.list.lines) === plots
        gui.pose_boxes[3].stored_string[] = "0"
        @test _value(gui, "Hits") == "1"
        # another detector keeps the sections, too
        ctrl.selected[] = o.pd2
        @test Ext._blocks!(Any[], insp.context) == blocks
        @test occursin("no panel", insp.context.content[1].content.content[2].content.text[])
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

    @testset "detector section" begin
        gui, o = _fixture()
        insp, ctrl = gui.layout.inspector, gui.controls
        p = only(gui.panels)
        ctrl.selected[] = o.pd
        section = insp.context.content[1].content
        segmented(i) = only(c.content for c in section.content
                            if c.content isa GridLayout && c.span.rows == i:i)
        mode_buttons = [c.content for c in segmented(3).content]
        scale_buttons = [c.content for c in segmented(4).content]
        @test p.mode == :auto
        # the panel mode and the color scale are set in place
        mode_buttons[2].clicks[] += 1
        @test p.mode == :spot
        @test p.scatter_plot.visible[]
        mode_buttons[3].clicks[] += 1
        @test p.mode == :intensity
        @test p.heat_plot.visible[]
        scale_buttons[2].clicks[] += 1
        @test p.colorscale == :log
        # the section shows the state of the panel when selected again
        ctrl.selected[] = nothing
        p.colorscale = :linear
        ctrl.selected[] = o.pd
        section = insp.context.content[1].content
        @test [c.content for c in segmented(4).content][1].buttoncolor[] ==
              gui.layout.theme.accent_soft
        close(gui)
    end

    @testset "registry" begin
        gui, o = _fixture()
        insp, ctrl = gui.layout.inspector, gui.controls
        built = Ref(0)
        updated = Any[]
        Ext._add_inspector!(gui, BMO.AbstractReflectiveOptic, function (gui, grid)
            built[] += 1
            Label(grid[1, 1], "Mine")
            return obj -> push!(updated, obj)
        end)
        ctrl.selected[] = o.m
        @test built[] == 1
        @test updated == [o.m]
        @test first.(insp.sections) == [BMO.AbstractReflectiveOptic]
        # other types, other sections
        ctrl.selected[] = o.pd
        @test first.(insp.sections) == [BMO.Detector]
        close(gui)
        # the compact layout has no inspector sections
        gui = live_view(System([o.m]), o.beam; trace_budget = Inf)
        @test_throws ArgumentError Ext._add_inspector!(gui, BMO.Mirror, (g, grid) -> identity)
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
