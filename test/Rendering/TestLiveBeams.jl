module TestLiveBeams

using BeamletOptics
using Makie
using Test
using LinearAlgebra
using GeometryBasics: coordinates, faces, Point3f

const BMO = BeamletOptics

# Segment geometry of the static renderer
function _expected_segment(ray; flen)
    isect = BMO.intersection(ray)
    len = isnothing(isect) ? flen : length(isect)
    p0 = position(ray)
    p1 = p0 + len * BMO.direction(ray)
    return (Point3f(p0), Point3f(p1))
end

function _expected_beam_points(beam; flen)
    pts = Point3f[]
    for child in BMO.PreOrderDFS(beam)
        for ray in BMO.rays(child)
            p0, p1 = _expected_segment(ray; flen)
            push!(pts, p0, p1)
        end
    end
    return pts
end

# Lens + mirror system
function _lens_mirror_fixture()
    lens = SphericalLens(50e-3, -50e-3, 10e-3, 25e-3, 1.5168)
    mir = RoundPlanoMirror(25.4e-3, 5e-3)
    translate3d!(mir, [0, 60e-3, 0])
    zrotate3d!(mir, deg2rad(45))
    sys = System([lens, mir])
    beam = Beam([0.0, -20e-3, 0.0], [0.0, 1.0, 0.0])
    solve_system!(sys, beam)
    return sys, lens, mir, beam
end

@testset "Live beam rendering" begin
    Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)
    @test !isnothing(Ext)

    @testset "Beam: segment endpoints match static geometry" begin
        sys, lens, mir, beam = _lens_mirror_fixture()
        fig = Figure()
        ax = LScene(fig[1, 1])
        flen = 1.0
        h = Ext.live_render!(ax, beam; flen)
        @test h isa Ext.BeamRenderHandle
        @test h.points[] == _expected_beam_points(beam; flen)
        @test length(h.points[]) == 2 * length(BMO.rays(beam))
    end

    @testset "Beam: update_render! after moving a component" begin
        sys, lens, mir, beam = _lens_mirror_fixture()
        fig = Figure()
        ax = LScene(fig[1, 1])
        flen = 1.0
        h = Ext.live_render!(ax, beam; flen)
        n0 = length(ax.scene.plots)

        # Move the mirror sideways (still in the system, just to a new pose)
        translate3d!(mir, [5e-3, 0, 0])
        solve_system!(sys, beam)
        Ext.update_render!(h)

        @test length(ax.scene.plots) == n0 # no new plots created
        @test h.points[] == _expected_beam_points(beam; flen)
    end

    @testset "Beam: topology change keeps a single plot" begin
        sys, lens, mir, beam = _lens_mirror_fixture()
        fig = Figure()
        ax = LScene(fig[1, 1])
        flen = 1.0
        h = Ext.live_render!(ax, beam; flen)
        n0 = length(ax.scene.plots)
        n_segments_before = length(BMO.rays(beam))

        # Move the mirror completely out of the beam path
        translate3d!(mir, [10.0, 0, 0])
        solve_system!(sys, beam)
        n_segments_after = length(BMO.rays(beam))
        @test n_segments_after != n_segments_before

        Ext.update_render!(h)
        @test length(ax.scene.plots) == n0
        @test length(h.points[]) == 2 * n_segments_after
        @test h.points[] == _expected_beam_points(beam; flen)
    end

    @testset "Ray: single-segment handle" begin
        fig = Figure()
        ax = LScene(fig[1, 1])
        ray = Ray([0.0, 0, 0], [1.0, 0, 0])
        h = Ext.live_render!(ax, ray; flen = 0.5)
        @test length(h.points[]) == 2
        @test h.points[][1] == Point3f(0, 0, 0)
        @test h.points[][2] == Point3f(0.5, 0, 0)
        Ext.update_render!(h)
        @test length(h.points[]) == 2
    end

    @testset "Beam group with render_every" begin
        src = CollimatedSource([0.0, 0, 0], [0.0, 1.0, 0.0], 10e-3; num_rings = 4, num_rays = 80)
        fig = Figure()
        ax = LScene(fig[1, 1])
        render_every = 5
        h = Ext.live_render!(ax, src; render_every, flen = 1.0)
        expected_beams = 1:render_every:length(BMO.beams(src))
        @test length(h.points[]) == 2 * length(expected_beams)

        # rebuild after re-solving (no system, but exercises update_render! on a group)
        Ext.update_render!(h)
        @test length(h.points[]) == 2 * length(expected_beams)
    end

    @testset "Beam group of non-Beams errors via the linesegments path" begin
        b1 = AstigmaticGaussianBeamlet([0, 0, 0], [0, 1, 0], 1000e-9, 1e-3)
        bg = AstigmaticBeamGroup([b1], [0, 0, 0], [0, 1, 0])
        # `_collect_segments!` (linesegments path) still rejects non-`Beam` members;
        # `live_render!(axis, ::AstigmaticBeamGroup)` itself dispatches to the dedicated
        # mesh method instead, see the "AstigmaticBeamGroup" testset below.
        @test_throws ArgumentError Ext._collect_segments!(Point3f[], bg; flen = 1.0)
    end

    @testset "remove_render! deletes the plot" begin
        _, _, _, beam = _lens_mirror_fixture()
        fig = Figure()
        ax = LScene(fig[1, 1])
        n0 = length(ax.scene.plots)
        h = Ext.live_render!(ax, beam)
        @test length(ax.scene.plots) == n0 + 1
        Ext.remove_render!(h)
        @test length(ax.scene.plots) == n0
    end

    @testset "show is a compact one-liner" begin
        _, _, _, beam = _lens_mirror_fixture()
        fig = Figure()
        ax = LScene(fig[1, 1])
        h = Ext.live_render!(ax, beam)
        str = sprint(show, h)
        @test !occursin('\n', str)
        @test occursin("BeamRenderHandle", str)
    end
end

# Michelson interferometer, see TestMichelson.jl
function _michelson_fixture()
    l_0 = 0.1
    m1 = SquarePlanoMirror2D(BMO.inch)
    m2 = SquarePlanoMirror2D(BMO.inch)
    bs = ThinBeamsplitter(BMO.inch, reflectance = 0.5)
    pd = Detector(BMO.inch / 5)
    translate3d!(m1, [l_0, 0, 0])
    translate3d!(m2, [0, l_0, 0])
    translate3d!(pd, [-l_0, 0, 0])
    zrotate3d!(bs, deg2rad(45))
    zrotate3d!(m1, deg2rad(90))
    zrotate3d!(pd, deg2rad(90))
    system = System([m1, m2, bs, pd])
    λ = 635e-9
    gauss = GaussianBeamlet([0, -l_0, 0], [0, 1.0, 0], λ, 1e-4, P0 = 5e-3)
    solve_system!(system, gauss)
    return system, m1, m2, bs, gauss
end

function _count_gaussian_segments(g)
    n = length(BMO.rays(g.chief))
    for c in g.children
        n += _count_gaussian_segments(c)
    end
    return n
end

@testset "Live Gaussian beamlet rendering" begin
    Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)

    @testset "Mesh vertex/face count matches #segments × resolution" begin
        _, _, _, _, gauss = _michelson_fixture()
        fig = Figure()
        ax = LScene(fig[1, 1])
        r_res, z_res = 12, 16
        h = Ext.live_render!(ax, gauss; r_res, z_res)
        @test h isa Ext.GaussianRenderHandle

        nseg = _count_gaussian_segments(gauss)
        m = h.mesh_obs[]
        @test length(coordinates(m)) == nseg * r_res * z_res
        @test length(faces(m)) == nseg * 2 * (r_res - 1) * (z_res - 1)
    end

    @testset "Vertex radius at segment start matches gauss_parameters waist" begin
        _, _, _, _, gauss = _michelson_fixture()
        fig = Figure()
        ax = LScene(fig[1, 1])
        r_res, z_res = 12, 16
        h = Ext.live_render!(ax, gauss; r_res, z_res)

        m = h.mesh_obs[]
        verts = coordinates(m)
        # First ring of the first segment (child = gauss itself, l = 0)
        ray1 = BMO.rays(gauss.chief)[1]
        w0 = BMO.gauss_parameters(gauss, [0.0])[1][1]
        ring = verts[1:r_res]
        p0 = position(ray1)
        radii = [norm(Vector(v) .- Vector(p0)) for v in ring]
        @test all(r -> isapprox(r, w0; atol = 1e-6), radii)
    end

    @testset "update_render! rebuilds mesh; plot count constant across updates" begin
        system, m1, m2, bs, gauss = _michelson_fixture()
        fig = Figure()
        ax = LScene(fig[1, 1])
        h = Ext.live_render!(ax, gauss; r_res = 10, z_res = 12)
        n0 = length(ax.scene.plots)

        for dl in (1e-3, -2e-3, 5e-4)
            translate3d!(m2, [0, dl, 0])
            solve_system!(system, gauss)
            Ext.update_render!(h)
            @test length(ax.scene.plots) == n0
        end
        nseg = _count_gaussian_segments(gauss)
        @test length(coordinates(h.mesh_obs[])) == nseg * 10 * 12
    end

    @testset "show_beams overlays 3 additional linesegments plots" begin
        _, _, _, _, gauss = _michelson_fixture()
        fig = Figure()
        ax = LScene(fig[1, 1])
        n0 = length(ax.scene.plots)
        h = Ext.live_render!(ax, gauss; show_beams = true, r_res = 8, z_res = 10)
        @test length(ax.scene.plots) == n0 + 4 # mesh + chief + divergence + waist
        @test length(h.beam_obs) == 3
        Ext.update_render!(h)
        @test length(ax.scene.plots) == n0 + 4
        Ext.remove_render!(h)
        @test length(ax.scene.plots) == n0
    end

    @testset "remove_render!" begin
        _, _, _, _, gauss = _michelson_fixture()
        fig = Figure()
        ax = LScene(fig[1, 1])
        n0 = length(ax.scene.plots)
        h = Ext.live_render!(ax, gauss; r_res = 8, z_res = 10)
        @test length(ax.scene.plots) == n0 + 1
        Ext.remove_render!(h)
        @test length(ax.scene.plots) == n0
    end

    @testset "show is a compact one-liner" begin
        _, _, _, _, gauss = _michelson_fixture()
        fig = Figure()
        ax = LScene(fig[1, 1])
        h = Ext.live_render!(ax, gauss; r_res = 8, z_res = 10)
        str = sprint(show, h)
        @test !occursin('\n', str)
        @test occursin("GaussianRenderHandle", str)
    end

end

# AGB through a thin lens, see test/Rendering/TestRenderPolarization.jl
function _agb_lens_fixture()
    agb = AstigmaticGaussianBeamlet([0, 0, 0], [0, 1, 0], 1000e-9, 5e-3; support = [0, 0, 1])
    l = ThinLens(50e-3, 50e-3, BMO.inch, 1.5)
    translate3d!(l, [0, 50e-3, 0])
    sys = System([l])
    solve_system!(sys, agb; check_invariant = false)
    return sys, l, agb
end

function _count_astigmatic_segments(agb)
    n = length(BMO.rays(agb.c))
    for c in agb.children
        n += _count_astigmatic_segments(c)
    end
    return n
end

@testset "Live AstigmaticGaussianBeamlet rendering" begin
    Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)

    @testset "Mesh vertex/face count matches #segments × resolution" begin
        _, _, agb = _agb_lens_fixture()
        fig = Figure()
        ax = LScene(fig[1, 1])
        r_res, z_res = 12, 16
        h = Ext.live_render!(ax, agb; r_res, z_res)
        @test h isa Ext.GaussianRenderHandle

        nseg = _count_astigmatic_segments(agb)
        m = h.mesh_obs[]
        @test length(coordinates(m)) == nseg * r_res * z_res
        @test length(faces(m)) == nseg * 2 * (r_res - 1) * (z_res - 1)
    end

    @testset "First ring radii match waist_parameters" begin
        _, _, agb = _agb_lens_fixture()
        fig = Figure()
        ax = LScene(fig[1, 1])
        r_res, z_res = 16, 20
        h = Ext.live_render!(ax, agb; r_res, z_res)

        m = h.mesh_obs[]
        verts = coordinates(m)
        ring = verts[1:r_res]
        (p0, b, c) = BMO.waist_parameters(agb, 0.0)
        expected = [Point3f(BMO.ellipse(v, p0, b, c)) for v in LinRange(0, 2π, r_res)]
        @test all(isapprox.(ring, expected; atol = 1e-6))
    end

    @testset "update_render! after moving a component and re-solving" begin
        sys, l, agb = _agb_lens_fixture()
        fig = Figure()
        ax = LScene(fig[1, 1])
        h = Ext.live_render!(ax, agb; r_res = 10, z_res = 12)
        n0 = length(ax.scene.plots)

        for dl in (5e-3, -8e-3)
            translate3d!(l, [0, dl, 0])
            solve_system!(sys, agb; check_invariant = false)
            Ext.update_render!(h)
            @test length(ax.scene.plots) == n0
        end
        nseg = _count_astigmatic_segments(agb)
        @test length(coordinates(h.mesh_obs[])) == nseg * 10 * 12
    end

    @testset "update_render! after translate3d! on the beamlet itself" begin
        sys, l, agb = _agb_lens_fixture()
        fig = Figure()
        ax = LScene(fig[1, 1])
        h = Ext.live_render!(ax, agb; r_res = 10, z_res = 12)
        n0 = length(ax.scene.plots)

        # Moving the beamlet source resets it (empty!); re-solve before updating.
        translate3d!(agb, [2e-3, 0, 0])
        solve_system!(sys, agb; check_invariant = false)
        Ext.update_render!(h)

        @test length(ax.scene.plots) == n0
        nseg = _count_astigmatic_segments(agb)
        @test length(coordinates(h.mesh_obs[])) == nseg * 10 * 12
    end

    @testset "show_beams overlays 3 additional linesegments plots" begin
        _, _, agb = _agb_lens_fixture()
        fig = Figure()
        ax = LScene(fig[1, 1])
        n0 = length(ax.scene.plots)
        h = Ext.live_render!(ax, agb; show_beams = true, r_res = 8, z_res = 10)
        @test length(ax.scene.plots) == n0 + 4 # mesh + chief + divergence + waist
        @test length(h.beam_obs) == 3
        Ext.update_render!(h)
        @test length(ax.scene.plots) == n0 + 4
        Ext.remove_render!(h)
        @test length(ax.scene.plots) == n0
    end

    @testset "remove_render!" begin
        _, _, agb = _agb_lens_fixture()
        fig = Figure()
        ax = LScene(fig[1, 1])
        n0 = length(ax.scene.plots)
        h = Ext.live_render!(ax, agb; r_res = 8, z_res = 10)
        @test length(ax.scene.plots) == n0 + 1
        Ext.remove_render!(h)
        @test length(ax.scene.plots) == n0
    end
end

@testset "Live AstigmaticBeamGroup rendering" begin
    Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)

    @testset "render_every merges every n-th beamlet into one mesh" begin
        bg = CollimatedGaussianBeamletSource([0, 0, 0], [0, 1, 0], 10e-3, 1000e-9, 2e-3; n_grid = 4)
        fig = Figure()
        ax = LScene(fig[1, 1])
        n0 = length(ax.scene.plots)
        render_every, r_res, z_res = 3, 6, 5
        h = Ext.live_render!(ax, bg; render_every, r_res, z_res)
        @test h isa Ext.AstigmaticGroupRenderHandle
        @test length(ax.scene.plots) == n0 + 1

        n_rendered = length(1:render_every:length(BMO.beams(bg)))
        m = h.mesh_obs[]
        @test length(coordinates(m)) == n_rendered * r_res * z_res

        Ext.update_render!(h)
        @test length(ax.scene.plots) == n0 + 1
        @test length(coordinates(h.mesh_obs[])) == n_rendered * r_res * z_res
    end

    @testset "remove_render!" begin
        bg = CollimatedGaussianBeamletSource([0, 0, 0], [0, 1, 0], 10e-3, 1000e-9, 2e-3; n_grid = 3)
        fig = Figure()
        ax = LScene(fig[1, 1])
        n0 = length(ax.scene.plots)
        h = Ext.live_render!(ax, bg; render_every = 2, r_res = 6, z_res = 5)
        @test length(ax.scene.plots) == n0 + 1
        Ext.remove_render!(h)
        @test length(ax.scene.plots) == n0
    end

    @testset "show is a compact one-liner" begin
        bg = CollimatedGaussianBeamletSource([0, 0, 0], [0, 1, 0], 10e-3, 1000e-9, 2e-3; n_grid = 3)
        fig = Figure()
        ax = LScene(fig[1, 1])
        h = Ext.live_render!(ax, bg; render_every = 2, r_res = 6, z_res = 5)
        str = sprint(show, h)
        @test !occursin('\n', str)
        @test occursin("AstigmaticGroupRenderHandle", str)
    end
end

end # module
