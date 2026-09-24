module TestBugFixes

using BeamletOptics
using LinearAlgebra
using Test
using Logging

const BMO = BeamletOptics

const mm = 1e-3

@testset "Issue#11" begin
    # https://github.com/JuliaPhysics/BeamletOptics.jl/issues/11
    # A very narrow point source produced a scattered spot diagram behind a lens system that had
    # been moved and rotated before (suspected rounding errors of `_world_to_sdf`).
    @testset "Narrow point source through a double Gauss lens" begin
        # Based on https://www.pencilofrays.com/double-gauss-sonnar-comparison/
        l1 = SphericalLens(48.88mm, 182.96mm, 8.89mm, 52.3mm, λ -> 1.62286)
        l23 = SphericalDoubletLens(36.92mm, Inf, 23.06mm, 15.11mm, 2.31mm,
            45.11mm, λ -> 1.58565, λ -> 1.67764)
        l45 = SphericalDoubletLens(-23.91mm, Inf, -36.92mm, 1.92mm, 7.77mm,
            40.01mm, λ -> 1.57046, λ -> 1.64128)
        l6 = SphericalLens(1063.24mm, -48.88mm, 6.73mm, 45.11mm, λ -> 1.62286)
        l_23 = thickness(l1) + 0.38mm
        l_45 = l_23 + thickness(l23) + 9.14mm + 13.36mm
        l_6 = l_45 + thickness(l45) + 0.38mm
        translate3d!(l23, [0, l_23, 0])
        translate3d!(l45, [0, l_45, 0])
        translate3d!(l6, [0, l_6, 0])
        detector = Detector(5mm)
        test_setup = ObjectGroup([ObjectGroup([l1, l23, l45, l6]), detector])
        # Move, rotate and reset the setup before tracing, as in the issue
        translate3d!(test_setup, [0.05, 0.05, 0.05])
        xrotate3d!(test_setup, deg2rad(60))
        zrotate3d!(test_setup, deg2rad(45))
        reset_rotation3d!(test_setup)
        reset_translation3d!(test_setup)
        translate_to3d!(detector, [0, 0.147, 0])
        source = PointSource([0, -0.5, 0], [0, 1, 0], 5e-5, 486.0e-9, num_rays = 1000,
            num_rings = 10)
        solve_system!(System([test_setup]), source)
        @test all(hit -> norm(hit) < 2e-7, spot_diagram(detector))
    end

    @testset "Narrow point source through a tilted concave asphere" begin
        # The numeric normal of aspheres, introduced for this issue, kicked near-axis rays out of
        # their meridional plane at concave aspheres, see Issue#80. A point source on the axis of
        # a rotationally symmetric lens: every ray must stay in its meridional plane.
        n = 1.458
        lens = Lens(CircularFlatSurface(30mm),
            EvenAsphericalSurface((n - 1) * 50mm, 30mm, -n^2, [0.0]), 4mm, x -> n)
        xrotate3d!(lens, deg2rad(0.5))
        zrotate3d!(lens, deg2rad(0.2))
        axis = BMO.orientation(lens)[:, 2]
        source = PointSource(Vector(position(lens) .- 0.1 .* axis), Vector(axis), 5e-5, 486.0e-9;
            num_rays = 1000, num_rings = 10)
        solve_system!(System(lens), source)
        # out-of-plane component of each outgoing ray, the on-axis ray has no meridional plane
        kicks = map(BMO.beams(source)) do beam
            m = cross(axis, BMO.direction(first(beam.rays)))
            norm(m) < 1e-15 ? 0.0 : abs(dot(BMO.direction(last(beam.rays)), normalize(m)))
        end
        @test length(kicks) == 1000
        @test maximum(kicks) < 1e-12
    end
end

@testset "Issue#14" begin
    pd_res = 1000
    pd_size = 10mm
    pd = Detector(pd_size)
    BMO.zrotate3d!(pd, deg2rad(90))
    BMO.translate3d!(pd, [0.46, 0, 0])
    # Setup beam
    y_0 = 0.2
    beam = BMO.GaussianBeamlet([0, y_0, 0], [0.46, -y_0, 0], 532e-9, 2.5mm, P0 = 10e-3)
    # Solve system
    system = BMO.System([pd])
    empty!(pd)
    BMO.solve_system!(system, beam)

    pd_pwr = optical_power(
        pd;
        # restore pre-v0.11 behavior, e.g. no autolims
        n = pd_res,
        x_min = -pd_size / 2,
        x_max = pd_size / 2,
        z_min = -pd_size / 2,
        z_max = pd_size / 2
    )

    @test pd_pwr≈10e-3 atol=1e-5
end

@testset "Issue#22 and Issue#23" begin
    # https://github.com/JuliaPhysics/BeamletOptics.jl/issues/22
    # https://github.com/JuliaPhysics/BeamletOptics.jl/issues/23

    mutable struct TestSubstrate{T, S <: BMO.AbstractShape{T}, N} <:
                   BMO.AbstractRefractiveOptic{T, N}
        const shape::S
        n::N
    end

    BMO.refractive_index(ts::TestSubstrate, ::Real) = ts.n
    set_index(ts::TestSubstrate, new) = (ts.n = new)
    get_index(ts::TestSubstrate) = ts.n

    "Shifts the phase of the beamlet by a specific amount in [rad]."
    function shift_phase(gb::BMO.GaussianBeamlet, phase::Real)
        BMO.electric_field!(gb, BMO.electric_field(gb) * exp(im * phase))
        return nothing
    end

    ref_signal(ϕ, A) = (cos(ϕ) + 1) / 2 * A

    # setup system for tests below
    splitter = CubeBeamsplitter(10mm, n -> 1)
    substrate_length = 10mm
    substrate = TestSubstrate(BMO.CylinderSDF(5mm, substrate_length / 2), 1.5)

    pd_size = 10mm
    pd_res = 250
    detector = Detector(pd_size)

    translate3d!(substrate, [0, -25mm, 0])
    translate3d!(detector, [0, 40mm, 0])

    system = System([substrate, splitter, detector])

    start_offset = 50mm

    @testset "Testing electric_field calculation - non-imaging ref. index change" begin
        indices = (1, 10, 100, 1000)
        for index in indices
            set_index(substrate, index)
            phi = LinRange(0, 2pi, 30)
            int = zeros(length(phi))
            for (i, p) in enumerate(phi)
                gb_prb = GaussianBeamlet([0, -start_offset, 0], [0, 1, 0], 1e-6, 0.5mm)
                gb_ref = GaussianBeamlet([start_offset, 0, 0], [-1, 0, 0], 1e-6, 0.5mm)
                empty!(detector)
                shift_phase(gb_ref, p)
                solve_system!(system, gb_prb)
                solve_system!(system, gb_ref)
                int[i] = BMO.optical_power(detector)
            end
            @test isapprox(BMO.visibility(int), 1, atol = 1e-2)
        end
    end

    @testset "Testing electric_field calculation - ref. index based phase shift" begin
        λ = 1e-6
        gb_prb = GaussianBeamlet([0, -start_offset, 0], [0, 1, 0], λ, 0.5mm)
        gb_ref = GaussianBeamlet([start_offset, 0, 0], [-1, 0, 0], λ, 0.5mm)

        n_lambdas = substrate_length / BMO.wavelength(gb_prb)

        n_factors = LinRange(0, 1, 50)
        pwr = zeros(length(n_factors))
        # Increase ref. index of substrate until one additional λ of OPL has been introduced
        for (i, nf) in enumerate(n_factors)
            set_index(substrate, 1 + 1 / n_lambdas * nf)
            empty!(detector)
            solve_system!(system, gb_prb)
            solve_system!(system, gb_ref)
            pd_pwr = optical_power(
                detector;
                # restore pre-v0.11 behavior, e.g. no autolims
                n = pd_res,
                x_min = -pd_size / 2,
                x_max = pd_size / 2,
                z_min = -pd_size / 2,
                z_max = pd_size / 2
            )
            @test isapprox(pd_pwr, ref_signal(2pi * nf, 2e-3), atol = 1e-8)
        end
        # Test if opl difference is indeed one λ
        delta = BMO.optical_path_length(gb_prb)
        delta -= BMO.optical_path_length(gb_ref)
        delta /= λ
        @test delta ≈ 1
    end

    @testset "Testing electric_field mutation during retracing" begin
        gb_prb = GaussianBeamlet([0, -start_offset, 0], [0, 1, 0], 1e-6, 0.5mm)
        gb_ref = GaussianBeamlet([start_offset, 0, 0], [-1, 0, 0], 1e-6, 0.5mm)
        phis = LinRange(0, 2pi, 50)
        pwr = zeros(length(phis))
        # Vary starting phase by 0...2pi via retracing
        for (i, phi) in enumerate(phis)
            empty!(detector)
            solve_system!(system, gb_prb)
            solve_system!(system, gb_ref)
            pd_pwr = optical_power(
                detector;
                # restore pre-v0.11 behavior, e.g. no autolims
                n = pd_res,
                x_min = -pd_size / 2,
                x_max = pd_size / 2,
                z_min = -pd_size / 2,
                z_max = pd_size / 2
            )
            @test isapprox(pd_pwr, ref_signal(phi, 2e-3), atol = 1e-8)
            shift_phase(gb_prb, step(phis))
        end
    end
end

@testset "Issue#51" begin
    # https://github.com/JuliaPhysics/BeamletOptics.jl/issues/51
    mirror = BeamletOptics.SphericalMirror(0.1, 0.01, 0.2)
    system = StaticSystem([mirror])
    beam = Beam([0, -0.19, 0.07], [0.0, 1, 0])
    solve_system!(system, beam)
    r1 = BMO.rays(beam)[1]
    r2 = BMO.rays(beam)[2]
    r3 = BMO.rays(beam)[3]
    @testset "Testing SDF sphere marching surface bug regression" begin
        @test length(BMO.rays(beam)) == 3
        @test BMO.object(BMO.intersection(r1)) === mirror
        @test BMO.object(BMO.intersection(r2)) === mirror
        @test isnothing(BMO.intersection(r3))
        @test dot(BMO.direction(r1), BMO.direction(r3)) < 0
    end
end

@testset "Meniscus lens on-axis normal" begin
    # Issue: AD normal of MeniscusLensSDF flipped on the optical axis, on-axis rays missed the lens.
    for (r1, r2) in ((69.21mm, 433.84mm), (-433.84mm, -69.21mm))
        lens = SphericalLens(r1, r2, 9.33mm, 70mm, λ -> 1.671)
        @test BMO.shape(lens) isa BMO.MeniscusLensSDF
        # rotate to test transformation of the normal
        zrotate3d!(lens, deg2rad(30))
        dir = orientation(lens)[:, 2]
        beam = Beam(Ray(position(lens) - 0.05 * dir, dir))
        solve_system!(System([lens]), beam)
        @test length(BMO.rays(beam)) == 3
        @test BMO.refractive_index.(BMO.rays(beam)) == [1, 1.671, 1]
        @test abs(dot(BMO.direction(last(BMO.rays(beam))), dir)) ≈ 1
    end
end

@testset "Issue#80" begin
    # https://github.com/JuliaPhysics/BeamletOptics.jl/issues/80
    # Near-axis rays through hyperbolic plano-concave and plano-convex lenses.
    # A plano-conic lens with k = -n² images a collimated beam that enters through the plane
    # face exactly onto a point: every outgoing ray passes through the (real or virtual)
    # focus at `f` from the conic vertex. Near the vertex the concave SDF is a sliver thinner
    # than the finite-difference step of `numeric_gradient`, which used to give wrong normals
    # for |h| ≲ 10 µm.
    n = 1.458
    f = 50mm
    ct = 4mm
    d = 30mm
    concave = Lens(CircularFlatSurface(d),
        EvenAsphericalSurface((n - 1) * f, d, -n^2, [0.0]), ct, x -> n)
    convex = Lens(CircularFlatSurface(d),
        EvenAsphericalSurface(-(n - 1) * f, d, -n^2, [0.0]), 2ct, x -> n)
    # The conic vertex lies on the back face, `thickness` along the lens axis
    for (lens, thickness, f_signed) in ((concave, ct, -f), (convex, 2ct, f)),
        tilt in (0.0, deg2rad(0.5))

        xrotate3d!(lens, tilt)
        system = System(lens)
        e_x, axis = BMO.orientation(lens)[:, 1], BMO.orientation(lens)[:, 2]
        focus = position(lens) .+ (thickness + f_signed) .* axis
        # Below h ≈ 0.1 µm the sag (h²/2R ≈ 1e-13 m) is smaller than the ray-marching
        # tolerance, the plane face is hit instead; the angle error stays below h / f.
        for h in (1e-6, 1e-5, 1e-4, 1e-3, 5mm)
            beam = Beam(Ray(position(lens) .- 0.1 .* axis .+ h .* e_x, Vector(axis)))
            solve_system!(system, beam)
            p, out = position(last(beam.rays)), BMO.direction(last(beam.rays))
            # lateral miss of the focus by the (extended) outgoing ray, relative to h
            t = dot(focus .- p, axis) / dot(out, axis)
            miss = dot(p .+ t .* out .- focus, e_x)
            @test abs(miss) < 1e-6 * h
        end
        xrotate3d!(lens, -tilt)
    end
end

@testset "PlateBeamsplitter BoundsError" begin
    # Issue: Newborn child beams only contain one ray. Indexing rays(beam)[id]
    # where id > 1 throws a BoundsError in interact3d for AstigmaticGaussianBeamlet.

    # 1 lens followed by 1 plate beamsplitter
    lens = ThinLens(Inf, 50mm, 25mm, 1.5)
    BMO.translate3d!(lens, [0, 50mm, 0])

    pbs = RoundPlateBeamsplitter(25mm, 5mm, λ -> 1.5)
    BMO.translate3d!(pbs, [0, 100mm, 0])

    system = System([lens, pbs])

    # Astigmatic beamlet starting at origin
    agb = AstigmaticGaussianBeamlet([0, 0, 0], [0, 1, 0], 1e-6, 1mm)

    # This should not throw BoundsError
    @test_nowarn solve_system!(system, agb)

    # Verify we have children (splitting happened)
    @test length(BMO.children(agb)) == 2
end

@testset "Retrace Tail Trimming Consistency" begin
    # Issue: When retracing reaches the old tail (i == n_c), a new segment
    # can be pushed and then rejected by check_optical_invariant.
    # Stale n_c caused inconsistent trimming between chief and aux beams.

    # We need a system where we can trigger an invariant violation during retrace.
    # We'll use a very small invariant threshold to make it easy to trigger.
    strict_threshold = 1e-15

    # A simple plane surface
    surf = RoundPlanoMirror(50mm, 5mm)
    BMO.translate3d!(surf, [0, 50mm, 0])
    system = System([surf])

    # Initial trace: 1 segment (2 rays)
    agb = AstigmaticGaussianBeamlet([0, 0, 0], [0, 1, 0], 1e-6, 1mm)

    with_logger(NullLogger()) do
        # Use strict threshold for this call
        solve_system!(system, agb; threshold = strict_threshold)

        n_initial = length(BMO.rays(agb.c))
        @test n_initial == 2

        # Now move the surface further away such that retracing will
        # reach the old tail and try to push a new segment.
        BMO.translate3d!(surf, [0, 10mm, 0])

        # Also tilt it or change parameters to trigger invariant violation if possible.
        BMO.retrace_system!(system, agb; threshold = strict_threshold)
    end

    # Verification of consistency
    lengths = map(b -> length(BMO.rays(b)), BMO._component_beams(agb))
    @test all(==(lengths[1]), lengths)
end

end # MODULE
