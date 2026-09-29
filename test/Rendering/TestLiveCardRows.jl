module TestLiveCardRows

using BeamletOptics
using Makie
using Test

const BMO = BeamletOptics

@testset "Live card rows" begin
    Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)
    @test !isnothing(Ext)

    _live_view(args...; kwargs...) =
        live_view(args...; merge((; trace_budget = Inf, throttle = false), kwargs)...)
    function _select!(gui, obj)
        gui.controls.selected[] = obj
        Ext._update_selection_box!(gui.controls)
        return nothing
    end
    # Text of the declared label `name` on the card of the selection
    _text(gui, name) = Ext._card_widget(gui.cards.selection, name).text[]

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

    @testset "mirror, detector, missed component" begin
        m, pd = _fixture()
        far = RoundPlanoMirror(25e-3, 5e-3)
        translate3d!(far, [0.3, 0.3, 0.3])
        gui = _live_view(System([m, pd, far]), Beam([0.0, 0, 0], [0.0, 1, 0]))
        _select!(gui, m)
        @test _text(gui, :beam) == "1 ray, AOI 45.0°"
        _select!(gui, far)
        @test _text(gui, :beam) == "not hit"
        _select!(gui, pd)
        @test _text(gui, :beam) == "1 ray, AOI 0.0°"
        # the number of rays of the spot panel
        @test _text(gui, :signal) == "N = 1"
        # refreshed after the solve of a change: the mirror rotated by 1° about z
        _select!(gui, m)
        Ext._card_widget(gui.cards.selection, :rv).stored_string[] = string(1e3 * deg2rad(1))
        @test _text(gui, :beam) in ("1 ray, AOI 44.0°", "1 ray, AOI 46.0°")
        close(gui)
    end

    @testset "lens, beamsplitters, polarizer" begin
        lens = SphericalLens(0.05, -0.05, 5e-3, 25.4e-3, 1.5)
        translate3d!(lens, [0, 0.05, 0])
        bs = RoundThinBeamsplitter(25e-3; reflectance = 0.3)
        translate3d!(bs, [0, 0.1, 0])
        cube = CubeBeamsplitter(10e-3, λ -> 1.5)
        translate3d!(cube, [0, 0.2, 0])
        pol = RoundLinearPolarizer(25e-3, 1e-3, 1e-3, λ -> 1.5)
        translate3d!(pol, [0, 0.3, 0])
        gui = _live_view(System([lens, bs, cube, pol]), Beam([0.0, 0, 0], [0.0, 1, 0], 633e-9))
        _select!(gui, lens)
        @test _text(gui, :index) == "1.5 at 633 nm"
        @test _text(gui, :thickness) == "5 mm"
        # the front and the back surface on the axis (and reflections of the splitters behind)
        hits = match(r"^(\d+) rays, AOI 0.0°", _text(gui, :beam))
        @test !isnothing(hits) && parse(Int, hits[1]) >= 2
        _select!(gui, bs)
        @test _text(gui, :split) == "R 30.0 %, T 70.0 %"
        @test _text(gui, :beam) != "not hit"
        # a multi-shape component: the hits of its prisms and its coating, behind the split
        _select!(gui, cube)
        @test _text(gui, :split) == "R 50.0 %, T 50.0 %"
        @test _text(gui, :beam) != "not hit"
        # the transmission axis along the local x-axis, i.e. horizontal for an optical axis along y
        _select!(gui, pol)
        @test _text(gui, :axis) == "0.0° from horizontal"
        rotate3d!(pol, [0.0, 1, 0], deg2rad(30))
        Ext._update_inspector!(gui; force = true)
        @test _text(gui, :axis) in ("30.0° from horizontal", "150.0° from horizontal")
        close(gui)
    end

    @testset "sources" begin
        m, pd = _fixture()
        g = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 1e-6, 0.5e-3)
        gui = _live_view(System([m, pd]), g)
        _select!(gui, g)
        @test _text(gui, :gauss) == "1000 nm, w0 500 µm, zR $(Ext._length_string(rayleigh_range(g)))"
        # the power of the intensity panel of the detector
        _select!(gui, pd)
        @test startswith(_text(gui, :signal), "P = ") && endswith(_text(gui, :signal), " mW")
        close(gui)

        src = CollimatedSource([0.0, 0, 0], [0.0, 1, 0], 4e-3, 633e-9; num_rings = 2, num_rays = 40)
        gui = _live_view(System([m, pd]), src; preview = false)
        _select!(gui, src)
        @test _text(gui, :source) == "633 nm, ⌀ 4 mm"
        # the rays of the source hit the mirror: the center ray first, the range of all of them
        _select!(gui, m)
        @test startswith(_text(gui, :beam), "40 rays, AOI 45.0° (")
        close(gui)

        ps = PointSource([0.0, 0, 0], [0.0, 1, 0], 0.01, 633e-9; num_rings = 2, num_rays = 40)
        gui = _live_view(System([m, pd]), ps)
        _select!(gui, ps)
        @test _text(gui, :source) == "633 nm, NA 0.01"
        close(gui)
    end

    @testset "systems" begin
        m, pd = _fixture()
        src = CollimatedSource([0.0, 0, 0], [0.0, 1, 0], 4e-3, 633e-9; num_rings = 2, num_rays = 40)
        sys = System([m, pd])
        gui = _live_view(sys => src, sys => Beam([0.0, 0, 0], [0.0, 1, 0]); preview = false)
        # no pose rows, the number of objects, the rays of both sources and the last solve
        rows = BMO.card_rows(sys)
        @test length(rows) == 3
        Ext._inspect!(gui, sys)
        @test _text(gui, :objects) == "2"
        @test _text(gui, :rays) == "41"
        @test _text(gui, :solve) == Ext._ms_string(gui.trace.solve_time)
        @test isnothing(Ext._card_widget(gui.cards.selection, :x))
        close(gui)
    end

    @testset "refresh of 1000 rays" begin
        m, pd = _fixture()
        src = UniformDiscSource([0.0, 0, 0], [0.0, 1, 0], 4e-3, 633e-9; num_rays = 1000)
        gui = _live_view(System([m, pd]), src; preview = false)
        _select!(gui, m)
        @test startswith(_text(gui, :beam), "1000 rays")
        Ext._update_inspector!(gui)
        t = @elapsed Ext._update_inspector!(gui)
        @test t < 5e-3
        close(gui)
    end
end

end
