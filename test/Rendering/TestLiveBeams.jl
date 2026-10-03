module TestLiveBeams

using BeamletOptics
using Makie
using Test
using LinearAlgebra
using GeometryBasics: coordinates, faces, Point3f

const BMO = BeamletOptics

const Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)

# Segment geometry of a ray
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

"""The points or the mesh that the plot `p` draws."""
_data(p) = p[1][]

"""The segment points of the `linesegments` plot of the ray handle `h`."""
_points(h) = _data(only(filter(p -> p isa Makie.LineSegments, BMO.render_plots(h))))

"""The envelope mesh of the Gaussian handle `h`."""
_mesh(h) = _data(only(filter(p -> p isa Makie.Mesh, BMO.render_plots(h))))

"""Returns the new plots of `render!(ax, x; kwargs...)` and of `live_render!(ax, x; kwargs...)`."""
function _static_and_live_plots(x; kwargs...)
    ax = LScene(Figure()[1, 1])
    n0 = length(ax.scene.plots)
    render!(ax, x; kwargs...)
    static = ax.scene.plots[(n0 + 1):end]
    ax = LScene(Figure()[1, 1])
    h = live_render!(ax, x; kwargs...)
    return static, BMO.render_plots(h), h
end

"""Checks that the `plots` are of the plot `types`, e.g. `Makie.Mesh`."""
_isa_all(plots, types) = length(plots) == length(types) && all(isa.(plots, types))

_same_data(a::AbstractVector, b::AbstractVector) = isequal(a, b)
_same_data(a, b) = coordinates(a) == coordinates(b) && faces(a) == faces(b)

"""Tests that `render!` and `live_render!` of `x` create the same plot types with the same geometry."""
function _test_same_plots(x; kwargs...)
    static, live, h = _static_and_live_plots(x; kwargs...)
    @test map(typeof, static) == map(typeof, live)
    @test all(_same_data(_data(s), _data(l)) for (s, l) in zip(static, live))
    return static, live, h
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

@testset "Live beam rendering" begin
    @test !isnothing(Ext)

    @testset "Beam: segment endpoints match the ray geometry" begin
        sys, lens, mir, beam = _lens_mirror_fixture()
        ax = LScene(Figure()[1, 1])
        flen = 1.0
        h = live_render!(ax, beam; flen)
        @test h isa BMO.AbstractBeamRenderHandle
        @test _points(h) == _expected_beam_points(beam; flen)
        @test length(_points(h)) == 2 * length(BMO.rays(beam))
    end

    @testset "Beam: update_render! after moving a component" begin
        sys, lens, mir, beam = _lens_mirror_fixture()
        ax = LScene(Figure()[1, 1])
        flen = 1.0
        h = live_render!(ax, beam; flen)
        n0 = length(ax.scene.plots)

        # Move the mirror sideways (still in the system, just to a new pose)
        translate3d!(mir, [5e-3, 0, 0])
        solve_system!(sys, beam)
        update_render!(h)

        @test length(ax.scene.plots) == n0 # no new plots created
        @test _points(h) == _expected_beam_points(beam; flen)
    end

    @testset "Beam: topology change keeps a single plot" begin
        sys, lens, mir, beam = _lens_mirror_fixture()
        ax = LScene(Figure()[1, 1])
        flen = 1.0
        h = live_render!(ax, beam; flen)
        n0 = length(ax.scene.plots)
        n_segments_before = length(BMO.rays(beam))

        # Move the mirror completely out of the beam path
        translate3d!(mir, [10.0, 0, 0])
        solve_system!(sys, beam)
        n_segments_after = length(BMO.rays(beam))
        @test n_segments_after != n_segments_before

        update_render!(h)
        @test length(ax.scene.plots) == n0
        @test length(_points(h)) == 2 * n_segments_after
        @test _points(h) == _expected_beam_points(beam; flen)
    end

    @testset "Ray: single-segment handle" begin
        ax = LScene(Figure()[1, 1])
        ray = Ray([0.0, 0, 0], [1.0, 0, 0])
        h = live_render!(ax, ray; flen = 0.5)
        @test length(_points(h)) == 2
        @test _points(h)[1] == Point3f(0, 0, 0)
        @test _points(h)[2] == Point3f(0.5, 0, 0)
        update_render!(h)
        @test length(_points(h)) == 2
    end

    @testset "Beam group with render_every" begin
        src = CollimatedSource([0.0, 0, 0], [0.0, 1.0, 0.0], 10e-3; num_rings = 4, num_rays = 80)
        ax = LScene(Figure()[1, 1])
        render_every = 5
        h = live_render!(ax, src; render_every, flen = 1.0)
        expected_beams = 1:render_every:length(BMO.beams(src))
        @test length(_points(h)) == 2 * length(expected_beams)

        # rebuild after re-solving (no system, but exercises update_render! on a group)
        update_render!(h)
        @test length(_points(h)) == 2 * length(expected_beams)
    end

    @testset "Beam group of non-Beams errors via the linesegments path" begin
        b1 = AstigmaticGaussianBeamlet([0, 0, 0], [0, 1, 0], 1000e-9, 1e-3)
        bg = AstigmaticBeamGroup([b1], [0, 0, 0], [0, 1, 0])
        # the linesegments geometry still rejects non-`Beam` members; `render!` and `live_render!`
        # of an `AstigmaticBeamGroup` dispatch to the mesh method instead, see below
        @test_throws ArgumentError Ext._ray_segments(bg; flen = 1.0)
    end

    @testset "remove_render! deletes the plot" begin
        _, _, _, beam = _lens_mirror_fixture()
        ax = LScene(Figure()[1, 1])
        n0 = length(ax.scene.plots)
        h = live_render!(ax, beam)
        @test length(ax.scene.plots) == n0 + 1
        remove_render!(h)
        @test length(ax.scene.plots) == n0
    end

    @testset "show is a compact one-liner" begin
        _, _, _, beam = _lens_mirror_fixture()
        ax = LScene(Figure()[1, 1])
        h = live_render!(ax, beam)
        str = sprint(show, h)
        @test !occursin('\n', str)
        @test occursin("BeamRenderHandle", str)
    end

    @testset "show_pos and polarization overlays follow update_render!" begin
        E0 = [0, 0, 1.0]
        mir = RoundPlanoMirror(25.4e-3, 5e-3)
        translate3d!(mir, [0, 60e-3, 0])
        zrotate3d!(mir, deg2rad(45))
        sys = System([mir])
        beam = Beam(PolarizedRay([0.0, 0, 0], [0.0, 1.0, 0], 1000e-9, E0))
        solve_system!(sys, beam)
        ax = LScene(Figure()[1, 1])
        n0 = length(ax.scene.plots)
        h = live_render!(ax, beam; show_pos = true, show_polarization = true, flen = 0.05)
        @test _isa_all(BMO.render_plots(h), [Makie.LineSegments, Makie.Scatter, Makie.Lines])
        @test length(ax.scene.plots) == n0 + 3
        _, pos, curve = BMO.render_plots(h)
        @test _data(pos) == Ext._ray_ends(beam)
        @test isequal(_data(curve), Ext._polarization_points(beam; flen = 0.05))

        translate3d!(mir, [0, 10e-3, 0])
        solve_system!(sys, beam)
        update_render!(h)
        @test length(ax.scene.plots) == n0 + 3
        @test _data(pos) == Ext._ray_ends(beam)
        @test isequal(_data(curve), Ext._polarization_points(beam; flen = 0.05))
        @test _points(h) == _expected_beam_points(beam; flen = 0.05)
        remove_render!(h)
        @test length(ax.scene.plots) == n0
    end
end

@testset "Live Gaussian beamlet rendering" begin
    @testset "Mesh vertex/face count matches #segments × resolution" begin
        _, _, _, _, gauss = _michelson_fixture()
        ax = LScene(Figure()[1, 1])
        r_res, z_res = 12, 16
        h = live_render!(ax, gauss; r_res, z_res)
        @test h isa BMO.AbstractBeamRenderHandle

        nseg = _count_gaussian_segments(gauss)
        m = _mesh(h)
        @test length(coordinates(m)) == nseg * r_res * z_res
        @test length(faces(m)) == nseg * 2 * r_res * (z_res - 1)
    end

    @testset "Vertex radius at segment start matches gauss_parameters waist" begin
        _, _, _, _, gauss = _michelson_fixture()
        ax = LScene(Figure()[1, 1])
        r_res, z_res = 12, 16
        h = live_render!(ax, gauss; r_res, z_res)

        verts = coordinates(_mesh(h))
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
        ax = LScene(Figure()[1, 1])
        h = live_render!(ax, gauss; r_res = 10, z_res = 12)
        n0 = length(ax.scene.plots)

        for dl in (1e-3, -2e-3, 5e-4)
            translate3d!(m2, [0, dl, 0])
            solve_system!(system, gauss)
            update_render!(h)
            @test length(ax.scene.plots) == n0
        end
        nseg = _count_gaussian_segments(gauss)
        @test length(coordinates(_mesh(h))) == nseg * 10 * 12
    end

    @testset "show_beams overlays 3 additional linesegments plots" begin
        _, _, _, _, gauss = _michelson_fixture()
        ax = LScene(Figure()[1, 1])
        n0 = length(ax.scene.plots)
        h = live_render!(ax, gauss; show_beams = true, r_res = 8, z_res = 10)
        @test length(ax.scene.plots) == n0 + 4 # mesh + chief + divergence + waist
        @test count(p -> p isa Makie.LineSegments, BMO.render_plots(h)) == 3
        update_render!(h)
        @test length(ax.scene.plots) == n0 + 4
        remove_render!(h)
        @test length(ax.scene.plots) == n0
    end

    @testset "remove_render!" begin
        _, _, _, _, gauss = _michelson_fixture()
        ax = LScene(Figure()[1, 1])
        n0 = length(ax.scene.plots)
        h = live_render!(ax, gauss; r_res = 8, z_res = 10)
        @test length(ax.scene.plots) == n0 + 1
        remove_render!(h)
        @test length(ax.scene.plots) == n0
    end

    @testset "show is a compact one-liner" begin
        _, _, _, _, gauss = _michelson_fixture()
        ax = LScene(Figure()[1, 1])
        h = live_render!(ax, gauss; r_res = 8, z_res = 10)
        str = sprint(show, h)
        @test !occursin('\n', str)
        @test occursin("BeamRenderHandle", str)
    end
end

@testset "Live AstigmaticGaussianBeamlet rendering" begin
    @testset "Mesh vertex/face count matches #segments × resolution" begin
        _, _, agb = _agb_lens_fixture()
        ax = LScene(Figure()[1, 1])
        r_res, z_res = 12, 16
        h = live_render!(ax, agb; r_res, z_res)
        @test h isa BMO.AbstractBeamRenderHandle

        nseg = _count_astigmatic_segments(agb)
        m = _mesh(h)
        @test length(coordinates(m)) == nseg * r_res * z_res
        @test length(faces(m)) == nseg * 2 * r_res * (z_res - 1)
    end

    @testset "First ring radii match waist_parameters" begin
        _, _, agb = _agb_lens_fixture()
        ax = LScene(Figure()[1, 1])
        r_res, z_res = 16, 20
        h = live_render!(ax, agb; r_res, z_res)

        ring = coordinates(_mesh(h))[1:r_res]
        (p0, b, c) = BMO.waist_parameters(agb, 0.0)
        # the ring runs counterclockwise about the beam, without a duplicate vertex at 2π
        d = BMO.direction(first(BMO.rays(agb.c)))
        b = dot(cross(b, c), d) < 0 ? -b : b
        expected = [Point3f(BMO.ellipse(v, p0, b, c)) for v in 2π .* (0:(r_res - 1)) ./ r_res]
        @test all(isapprox.(ring, expected; atol = 1e-6))
    end

    @testset "Closed rings, all faces point outwards" begin
        _, _, agb = _agb_lens_fixture()
        _, _, _, _, gauss = _michelson_fixture()
        # the chief ray segments in the order of the mesh, i.e. depth-first
        segments(b, chief) = vcat(BMO.rays(chief(b)), (segments(c, chief) for c in b.children)...)
        r_res, z_res = 12, 8
        nf = 2 * r_res * (z_res - 1)
        for (beam, chief) in ((agb, b -> b.c), (gauss, b -> b.chief))
            m = _mesh(live_render!(LScene(Figure()[1, 1]), beam; r_res, z_res))
            verts, fs = coordinates(m), faces(m)
            for (k, ray) in enumerate(segments(beam, chief))
                d = direction(ray)
                # every face normal has a positive component away from the segment axis
                outward = map(fs[((k - 1) * nf + 1):(k * nf)]) do f
                    p1, p2, p3 = (Vector(verts[i]) for i in f)
                    v = (p1 + p2 + p3) / 3 - position(ray)
                    dot(cross(p2 - p1, p3 - p1), v - dot(v, d) * d) > 0
                end
                @test all(outward)
            end
        end
    end

    @testset "update_render! after moving a component and re-solving" begin
        sys, l, agb = _agb_lens_fixture()
        ax = LScene(Figure()[1, 1])
        h = live_render!(ax, agb; r_res = 10, z_res = 12)
        n0 = length(ax.scene.plots)

        for dl in (5e-3, -8e-3)
            translate3d!(l, [0, dl, 0])
            solve_system!(sys, agb; check_invariant = false)
            update_render!(h)
            @test length(ax.scene.plots) == n0
        end
        nseg = _count_astigmatic_segments(agb)
        @test length(coordinates(_mesh(h))) == nseg * 10 * 12
    end

    @testset "update_render! after translate3d! on the beamlet itself" begin
        sys, l, agb = _agb_lens_fixture()
        ax = LScene(Figure()[1, 1])
        h = live_render!(ax, agb; r_res = 10, z_res = 12)
        n0 = length(ax.scene.plots)

        # Moving the beamlet source resets it (empty!); re-solve before updating.
        translate3d!(agb, [2e-3, 0, 0])
        solve_system!(sys, agb; check_invariant = false)
        update_render!(h)

        @test length(ax.scene.plots) == n0
        nseg = _count_astigmatic_segments(agb)
        @test length(coordinates(_mesh(h))) == nseg * 10 * 12
    end

    @testset "show_beams overlays 3 additional linesegments plots" begin
        _, _, agb = _agb_lens_fixture()
        ax = LScene(Figure()[1, 1])
        n0 = length(ax.scene.plots)
        h = live_render!(ax, agb; show_beams = true, r_res = 8, z_res = 10)
        @test length(ax.scene.plots) == n0 + 4 # mesh + chief + divergence + waist
        @test count(p -> p isa Makie.LineSegments, BMO.render_plots(h)) == 3
        update_render!(h)
        @test length(ax.scene.plots) == n0 + 4
        remove_render!(h)
        @test length(ax.scene.plots) == n0
    end

    @testset "show_waist and polarization follow update_render!" begin
        sys, l, agb = _agb_lens_fixture()
        ax = LScene(Figure()[1, 1])
        h = live_render!(ax, agb; show_waist = true, show_polarization = true, r_res = 8, z_res = 10)
        _, waist, curve = BMO.render_plots(h)
        @test waist isa Makie.Scatter && curve isa Makie.Lines
        translate3d!(l, [0, 5e-3, 0])
        solve_system!(sys, agb; check_invariant = false)
        update_render!(h)
        @test _data(waist) == coordinates(_mesh(h))
        @test isequal(_data(curve), Ext._polarization_points(agb; flen = 0.1))
        remove_render!(h)
    end

    @testset "remove_render!" begin
        _, _, agb = _agb_lens_fixture()
        ax = LScene(Figure()[1, 1])
        n0 = length(ax.scene.plots)
        h = live_render!(ax, agb; r_res = 8, z_res = 10)
        @test length(ax.scene.plots) == n0 + 1
        remove_render!(h)
        @test length(ax.scene.plots) == n0
    end
end

@testset "Live AstigmaticBeamGroup rendering" begin
    @testset "render_every merges every n-th beamlet into one mesh" begin
        bg = CollimatedGaussianBeamletSource([0, 0, 0], [0, 1, 0], 10e-3, 1000e-9, 2e-3; n_grid = 4)
        ax = LScene(Figure()[1, 1])
        n0 = length(ax.scene.plots)
        render_every, r_res, z_res = 3, 6, 5
        h = live_render!(ax, bg; render_every, r_res, z_res)
        @test h isa BMO.AbstractBeamRenderHandle
        @test length(ax.scene.plots) == n0 + 1

        n_rendered = length(1:render_every:length(BMO.beams(bg)))
        @test length(coordinates(_mesh(h))) == n_rendered * r_res * z_res

        update_render!(h)
        @test length(ax.scene.plots) == n0 + 1
        @test length(coordinates(_mesh(h))) == n_rendered * r_res * z_res
    end

    @testset "remove_render!" begin
        bg = CollimatedGaussianBeamletSource([0, 0, 0], [0, 1, 0], 10e-3, 1000e-9, 2e-3; n_grid = 3)
        ax = LScene(Figure()[1, 1])
        n0 = length(ax.scene.plots)
        h = live_render!(ax, bg; render_every = 2, r_res = 6, z_res = 5)
        @test length(ax.scene.plots) == n0 + 1
        remove_render!(h)
        @test length(ax.scene.plots) == n0
    end

    @testset "show is a compact one-liner" begin
        bg = CollimatedGaussianBeamletSource([0, 0, 0], [0, 1, 0], 10e-3, 1000e-9, 2e-3; n_grid = 3)
        ax = LScene(Figure()[1, 1])
        h = live_render!(ax, bg; render_every = 2, r_res = 6, z_res = 5)
        str = sprint(show, h)
        @test !occursin('\n', str)
        @test occursin("BeamRenderHandle", str)
    end
end

@testset "render! and live_render! share one drawing path" begin
    _, _, _, beam = _lens_mirror_fixture()
    _, _, _, _, gauss = _michelson_fixture()
    _, _, agb = _agb_lens_fixture()
    src = CollimatedSource([0.0, 0, 0], [0.0, 1.0, 0.0], 10e-3; num_rings = 4, num_rays = 80)
    bg = CollimatedGaussianBeamletSource([0, 0, 0], [0, 1, 0], 10e-3, 1000e-9, 2e-3; n_grid = 3)
    pbeam = Beam(PolarizedRay([0.0, 0, 0], [1.0, 0, 0], 1000e-9, [0, 0, 1.0]))
    pray = PolarizedRay([0.0, 0, 0], [1.0, 0, 0], 1000e-9, [0, 0, 1.0])

    @testset "$name" for (name, x, kw, types) in (
            ("ray", first(BMO.rays(beam)), (; flen = 0.2, show_pos = true), [Makie.LineSegments, Makie.Scatter]),
            ("Beam", beam, (; flen = 0.2, show_pos = true), [Makie.LineSegments, Makie.Scatter]),
            ("beam group of rays", src, (; render_every = 3, flen = 0.2), [Makie.LineSegments]),
            ("polarized ray", pray, (; flen = 0.1, show_polarization = true), [Makie.LineSegments, Makie.Lines]),
            ("polarized Beam", pbeam, (; flen = 0.1, show_polarization = true, pol_λ = 0.01),
                [Makie.LineSegments, Makie.Lines]),
            ("GaussianBeamlet", gauss, (; r_res = 8, z_res = 10, show_beams = true, show_pos = true),
                [Makie.Mesh, Makie.LineSegments, Makie.Scatter, Makie.LineSegments, Makie.Scatter,
                    Makie.LineSegments, Makie.Scatter]),
            ("AstigmaticGaussianBeamlet", agb, (; r_res = 8, z_res = 10, show_waist = true, show_beams = true,
                show_polarization = true), [Makie.Mesh, Makie.Scatter, Makie.LineSegments,
                    Makie.LineSegments, Makie.LineSegments, Makie.Lines]),
            ("AstigmaticBeamGroup", bg, (; r_res = 6, z_res = 5, render_every = 2), [Makie.Mesh]),
        )
        static, live, h = _test_same_plots(x; kw...)
        @test _isa_all(live, types)
        @test BMO.rendered(h) === x
    end

    @testset "default resolutions of render! and live_render!" begin
        # render! draws finer than live_render! for GaussianBeamlets and beam groups, see the docstrings
        nseg = _count_gaussian_segments(gauss)
        static, live, _ = _static_and_live_plots(gauss)
        @test length(coordinates(_data(static[1]))) == nseg * 50 * 100
        @test length(coordinates(_data(live[1]))) == nseg * 24 * 40
        static, live, _ = _static_and_live_plots(bg)
        n_rendered = length(1:5:length(BMO.beams(bg)))
        @test length(coordinates(_data(static[1]))) == n_rendered * 64 * 100
        @test length(coordinates(_data(live[1]))) == n_rendered * 10 * 8
    end

    @testset "kwargs of render!" begin
        static, _, _ = _static_and_live_plots(beam; color = :red, linewidth = 3.0, alpha = 0.3)
        @test only(static).linewidth[] == 3.0
        @test only(static).alpha[] == 0.3
        @test_throws ArgumentError render!(LScene(Figure()[1, 1]), beam; show_polarization = true)
        @test_throws ArgumentError live_render!(LScene(Figure()[1, 1]), src; show_polarization = true)
    end
end

@testset "Render handle protocol of beam handles" begin
    _, _, _, beam = _lens_mirror_fixture()
    _, _, _, _, gauss = _michelson_fixture()
    _, _, agb = _agb_lens_fixture()
    src = CollimatedSource([0.0, 0, 0], [0.0, 1.0, 0.0], 10e-3; num_rings = 4, num_rays = 80)
    bg = CollimatedGaussianBeamletSource([0, 0, 0], [0, 1, 0], 10e-3, 1000e-9, 2e-3; n_grid = 3)
    ax = LScene(Figure()[1, 1])
    for (x, kw, settings) in (
            (beam, (; flen = 0.3), (; flen = 0.3, render_every = 1)),
            (src, (; render_every = 3), (; flen = 1.0, render_every = 3)),
            (gauss, (; r_res = 8, z_res = 10), (; flen = 0.1, render_every = 1, r_res = 8, z_res = 10)),
            (agb, (;), (; flen = 0.1, render_every = 1, r_res = 64, z_res = 100)),
            (bg, (;), (; flen = 0.1, render_every = 5, r_res = 10, z_res = 8)),
        )
        h = live_render!(ax, x; kw...)
        @test h isa BMO.AbstractBeamRenderHandle
        @test BMO.rendered(h) === x
        @test BMO.render_settings(h) == settings
        @test !isempty(BMO.render_plots(h))
        @test all(p -> any(q -> q === p, ax.scene.plots), BMO.render_plots(h))
        remove_render!(h)
        @test isempty(BMO.render_plots(h))
    end
end

@testset "render_settings! of beam handles" begin
    _, _, _, beam = _lens_mirror_fixture()
    _, _, _, _, gauss = _michelson_fixture()
    _, _, agb = _agb_lens_fixture()
    src = CollimatedSource([0.0, 0, 0], [0.0, 1.0, 0.0], 10e-3; num_rings = 4, num_rays = 80)
    bg = CollimatedGaussianBeamletSource([0, 0, 0], [0, 1, 0], 10e-3, 1000e-9, 2e-3; n_grid = 3)
    ax = LScene(Figure()[1, 1])

    @testset "rays, beams and groups" begin
        h = live_render!(ax, beam; flen = 0.3)
        plots = copy(BMO.render_plots(h))
        @test BMO.render_settings!(h; flen = 0.7) === h
        @test BMO.render_settings(h) == (; flen = 0.7, render_every = 1)
        @test _points(h) == _expected_beam_points(beam; flen = 0.7)
        @test BMO.render_plots(h) == plots
        # same as rendering with the new setting
        @test _points(h) == _points(live_render!(ax, beam; flen = 0.7))
        # render_every does not apply to a single beam
        BMO.render_settings!(h; render_every = 3)
        @test BMO.render_settings(h).render_every == 1

        h = live_render!(ax, src; render_every = 5, show_pos = true)
        n = length(_points(h))
        BMO.render_settings!(h; render_every = 1, flen = 2.0)
        @test BMO.render_settings(h) == (; flen = 2.0, render_every = 1)
        @test length(_points(h)) > n
        @test _points(h) == _points(live_render!(ax, src; render_every = 1, flen = 2.0))
    end

    @testset "Gaussian beamlets" begin
        h = live_render!(ax, gauss; r_res = 8, z_res = 10, show_beams = true)
        BMO.render_settings!(h; flen = 0.4, r_res = 12, z_res = 6)
        @test BMO.render_settings(h) == (; flen = 0.4, render_every = 1, r_res = 12, z_res = 6)
        ref = live_render!(ax, gauss; r_res = 12, z_res = 6, flen = 0.4, show_beams = true)
        @test _same_data(_mesh(h), _mesh(ref))
        # the generating rays follow flen
        gen = filter(p -> p isa Makie.LineSegments, BMO.render_plots(h))
        genref = filter(p -> p isa Makie.LineSegments, BMO.render_plots(ref))
        @test length(gen) == 3
        @test all(_data(a) == _data(b) for (a, b) in zip(gen, genref))

        h = live_render!(ax, agb)
        BMO.render_settings!(h; flen = 0.25, r_res = 16)
        @test _same_data(_mesh(h), _mesh(live_render!(ax, agb; flen = 0.25, r_res = 16)))

        h = live_render!(ax, bg; show_waist = true)
        BMO.render_settings!(h; render_every = 1)
        @test BMO.render_settings(h).render_every == 1
        @test _same_data(_mesh(h), _mesh(live_render!(ax, bg; render_every = 1, show_waist = true)))
    end

    @testset "errors" begin
        h = live_render!(ax, beam; flen = 0.3)
        @test_throws ArgumentError BMO.render_settings!(h; r_res = 10)
        @test_throws ArgumentError BMO.render_settings!(h; color = :red)
        @test_throws ArgumentError BMO.render_settings!(h; flen = -1.0)
        @test_throws ArgumentError BMO.render_settings!(h; flen = Inf)
        @test_throws ArgumentError BMO.render_settings!(h; render_every = 0)
        @test_throws ArgumentError BMO.render_settings!(h; flen = 0.5, r_res = 10)
        @test BMO.render_settings(h) == (; flen = 0.3, render_every = 1)
        @test BMO.render_settings!(h) === h
    end
end

end # module
