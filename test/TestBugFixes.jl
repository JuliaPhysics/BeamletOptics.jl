module TestBugFixes

using BeamletOptics
using LinearAlgebra
using Test
using Logging

const BMO = BeamletOptics

const mm = 1e-3

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
    mirror = BeamletOptics.ConcaveSphericalMirror(0.1, 0.01, 0.2)
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
    old_threshold = BMO.INVARIANT_THRESHOLD[]
    BMO.INVARIANT_THRESHOLD[] = 1e-15 # Extremely strict

    try
        # A simple plane surface
        surf = RoundPlanoMirror(50mm, 5mm)
        BMO.translate3d!(surf, [0, 50mm, 0])
        system = System([surf])

        # Initial trace: 1 segment (2 rays)
        agb = AstigmaticGaussianBeamlet([0, 0, 0], [0, 1, 0], 1e-6, 1mm)
        
        with_logger(NullLogger()) do
            solve_system!(system, agb)

            n_initial = length(BMO.rays(agb.c))
            @test n_initial == 2

            # Now move the surface further away such that retracing will
            # reach the old tail and try to push a new segment.
            BMO.translate3d!(surf, [0, 10mm, 0])

            # Also tilt it or change parameters to trigger invariant violation if possible.
            BMO.retrace_system!(system, agb)
        end

        # Verification of consistency
        lengths = map(b -> length(BMO.rays(b)), BMO._component_beams(agb))
        @test all(==(lengths[1]), lengths)

    finally
        BMO.INVARIANT_THRESHOLD[] = old_threshold
    end
end

end # MODULE
