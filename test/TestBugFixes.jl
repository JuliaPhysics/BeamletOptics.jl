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

    @testset "Testing electric_field mutation across repeated solves" begin
        gb_prb = GaussianBeamlet([0, -start_offset, 0], [0, 1, 0], 1e-6, 0.5mm)
        gb_ref = GaussianBeamlet([start_offset, 0, 0], [-1, 0, 0], 1e-6, 0.5mm)
        phis = LinRange(0, 2pi, 50)
        pwr = zeros(length(phis))
        # Vary starting phase by 0...2pi, re-solving the system each time
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
        @test !isnothing(BMO.intersection(r1))
        @test !isnothing(BMO.intersection(r2))
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

@testset "Coincident refractive-reflective boundary" begin
    # Test for coincident boundary bug where a mirror placed exactly at the lens/glass back boundary
    # is correctly processed and reflects rays instead of letting them leak/refract straight through.
    
    λ = 354.84e-9
    B1, B2, B3 = 0.6961663, 0.4079426, 0.8974794
    C1, C2, C3 = 0.0684043^2, 0.1162414^2, 9.896161^2
    SE = SellmeierEquation(B1, B2, B3, C1, C2, C3)
    
    d_glass = 17.5mm
    L_bs = 30mm
    
    # Lens element (glass block)
    glass_block = Lens(RectangularFlatSurface(L_bs), RectangularFlatSurface(L_bs), d_glass, SE)
    translate3d!(glass_block, [0.0, L_bs/2, 0.0])
    
    # Mirror placed EXACTLY at the back face of the glass block (no gap, gap = 0.0)
    mirror = SquarePlanoMirror(L_bs, 5mm)
    translate3d!(mirror, [0.0, L_bs/2 + d_glass, 0.0])
    
    # Detector placed to capture the returning reflected ray
    detector = Detector(16mm)
    translate3d!(detector, [0.0, -10mm, 0.0])
    
    # Beam pointing along +y, starting at y = 10 mm
    beam = GaussianBeamlet([0.0, 10mm, 0.0], [0.0, 1.0, 0.0], λ, 0.5mm)
    
    system = System([glass_block, mirror, detector])
    
    empty!(detector)
    solve_system!(system, beam)
    
    @test detector.hits !== nothing
    @test length(detector.hits) == 1
end

@testset "Auxiliary beamlet TIR before splitting" begin
    struct Splitting5050Coating end
    BMO.coating_behavior(::Splitting5050Coating, ray) = Splitting()
    BMO.get_jones_matrix(::Splitting5050Coating, θi, λ, n1, n2, is_reflected; from_front=true) = BMO.SPBasis(1/sqrt(2), 0, 0, 1/sqrt(2))

    mesh = BMO.QuadraticFlatMesh(100mm)
    coating = Coating(mesh, Splitting5050Coating(), normal_filter=[0.0, -1.0, 0.0])
    
    # Ray incident from n_incident = 1.5 to n_transmitted = 1.0 (glass to air)
    # Critical angle = asin(1/1.5) ≈ 41.8103 deg
    # Chief ray angle = 41.7 deg (refracts, near critical angle)
    θ = deg2rad(41.7)
    dir_c = [sin(θ), cos(θ), 0.0]
    
    # Waist size 0.001mm (1 um) gives divergence angle large enough for divergence ray to exceed critical angle
    w0 = 0.001mm
    λ = 1000e-9
    
    gauss = GaussianBeamlet([0.0, -10mm, 0.0], dir_c, λ, w0)
    
    n_inc = 1.5
    n_trans = 1.0
    normal = [0.0, -1.0, 0.0]
    
    tir_triggered = false
    cb = () -> (tir_triggered = true; return nothing)
    
    BMO._propagate_splitting_gaussian_beamlet(
        System([coating]), coating, Splitting5050Coating(), gauss, 1, n_inc, n_trans, normal, true, cb
    )
    @test tir_triggered

    agb = AstigmaticGaussianBeamlet([0.0, -10mm, 0.0], dir_c, λ, w0)
    tir_triggered_agb = false
    cb_agb = () -> (tir_triggered_agb = true; return nothing)
    
    BMO._propagate_splitting_astigmatic_beamlet(
        System([coating]), coating, Splitting5050Coating(), agb, 1, n_inc, n_trans, normal, true, cb_agb
    )
    @test tir_triggered_agb
end
end # MODULE
