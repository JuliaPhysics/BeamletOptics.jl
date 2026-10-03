module TestPathSegments

using BeamletOptics
using Test
using LinearAlgebra: norm

const BMO = BeamletOptics

# Michelson-like layout: beam along +y, 45 deg beamsplitter at the origin, mirrors at 0.1 m (+x) and
# 0.2 m (+y)
function michelson()
    bs = ThinBeamsplitter(BMO.inch, reflectance = 0.5)
    zrotate3d!(bs, deg2rad(45))
    m1 = SquarePlanoMirror2D(BMO.inch)
    translate3d!(m1, [0.1, 0, 0])
    zrotate3d!(m1, deg2rad(90))
    m2 = SquarePlanoMirror2D(BMO.inch)
    translate3d!(m2, [0, 0.2, 0])
    return System([bs, m1, m2])
end

@testset "path_segments" begin
    system = michelson()
    beam = Beam([0, -0.1, 0], [0, 1.0, 0])
    solve_system!(system, beam)
    segs = path_segments(beam; flen = 0.05)

    @testset "tree structure and lengths" begin
        @test first(segs).parent == 0
        @test first(segs).depth == 1
        @test first(segs).s_start == 0
        @test first(segs).start ≈ [0, -0.1, 0]
        @test any(s -> s.depth == 3, segs)
        for (i, s) in enumerate(segs)
            @test s.s_stop - s.s_start ≈ norm(s.stop - s.start)
            @test s.opl_stop - s.opl_start ≈ norm(s.stop - s.start) # vacuum
            @test s.λ == 1e-6
            if s.parent == 0
                @test s.s_start == 0
            else
                # a segment starts where its parent ends
                @test s.start ≈ segs[s.parent].stop
                @test s.s_start ≈ segs[s.parent].s_stop
                @test s.opl_start ≈ segs[s.parent].opl_stop
            end
        end
        # lengths add up to the beam methods of every branch
        # (the final segments are in the order of the leaf beams)
        leaves = collect(BMO.Leaves(beam))
        for (leaf, i) in zip(leaves, findall(s -> s.final, segs))
            @test segs[i].s_start ≈ length(leaf)
            @test segs[i].opl_start ≈ BMO.optical_path_length(leaf)
        end
    end

    @testset "final rays" begin
        finals = filter(s -> s.final, segs)
        @test length(finals) == length(collect(BMO.Leaves(beam)))
        @test all(s -> s.s_stop - s.s_start ≈ 0.05, finals)
        @test all(s -> norm(s.stop - s.start) ≈ 0.05, finals)
        @test path_segments(beam; flen = 1) != segs
    end

    @testset "untraced beam and refractive index" begin
        b = Beam([0, 0, 0], [1.0, 0, 0])
        s = only(path_segments(b))
        @test s.final && s.s_stop == 1.0
        lens = SphericalLens(0.1, -0.1, 5e-3, BMO.inch, 1.5)
        sys = System([lens])
        b2 = Beam([0, -0.1, 0], [0, 1.0, 0])
        solve_system!(sys, b2)
        segs2 = path_segments(b2; flen = 0.0)
        @test maximum(s -> s.opl_stop - s.opl_start - (s.s_stop - s.s_start), segs2) > 1e-3
        @test last(filter(s -> !s.final, segs2)).opl_stop ≈ BMO.optical_path_length(b2)
    end

    @testset "beamlets" begin
        g = GaussianBeamlet([0, -0.1, 0], [0, 1.0, 0], 1e-6, 1e-3)
        solve_system!(system, g)
        gs = path_segments(g; flen = 0.05)
        cs = path_segments(g.chief; flen = 0.05)
        # the chief beam alone is the first branch, the children hang on the beamlet
        @test cs == gs[1:length(cs)]
        @test length(gs) > 3
        @test maximum(s -> s.depth, gs) == 3

        a = AstigmaticGaussianBeamlet([0, -0.1, 0], [0, 1.0, 0], 1e-6, 1e-3)
        solve_system!(system, a)
        as = path_segments(a; flen = 0.05)
        @test path_segments(a.c; flen = 0.05) == as[1:length(a.c.rays)]
        @test length(as) == length(gs)
        @test eltype(as) == BMO.PathSegment{Float64}
    end

    @testset "beam groups" begin
        src = UniformDiscSource([0, -0.1, 0], [0, 1.0, 0], 2e-3, 1e-6; num_rays = 7)
        solve_system!(system, src)
        gsegs = path_segments(src; flen = 0.05)
        @test Set(s.beam for s in gsegs) == Set(1:length(src))
        for (i, b) in enumerate(BMO.beams(src))
            idx = findall(s -> s.beam == i, gsegs)
            ref = path_segments(b; flen = 0.05)
            @test length(idx) == length(ref)
            for (j, k) in enumerate(idx)
                # parent indices of a group refer to its whole vector
                @test gsegs[k].parent == (ref[j].parent == 0 ? 0 : ref[j].parent + first(idx) - 1)
                @test gsegs[k].stop == ref[j].stop
                @test gsegs[k].s_stop == ref[j].s_stop
            end
        end
    end
end

end # module
