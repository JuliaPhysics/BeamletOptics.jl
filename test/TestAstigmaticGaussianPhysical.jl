module TestAstigmaticGaussianPhysical

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics
const mm = 1e-3
const nm = 1e-9

@testset "AGB Physical Validation" begin

    # Use support=[1,0,0] so that w0_x → X-axis, w0_y → Z-axis
    @testset "Free-space field — elliptical profile (wx ≠ wy)" begin
        λ = 633e-9
        w0_x = 1e-4
        w0_y = 2e-4
        beam = AstigmaticGaussianBeamlet(
            [0.0, 0, 0], [0, 1, 0], λ, w0_x, w0_y;
            E0 = [0, 0, 1], support = [1, 0, 0])

        z_eval = 0.05
        zRx = π * w0_x^2 / λ
        zRy = π * w0_y^2 / λ
        wx_z = w0_x * sqrt(1 + (z_eval / zRx)^2)
        wy_z = w0_y * sqrt(1 + (z_eval / zRy)^2)

        I0 = BMO.intensity(beam, [0, 0, 0], z_eval)

        # X-axis profile → w0_x direction (s1=[1,0,0])
        rs_x = LinRange(-3wx_z, 3wx_z, 51)
        Ix = [BMO.intensity(beam, [r, 0, 0], z_eval) for r in rs_x]
        Ix_ana = [I0 * exp(-2 * r^2 / wx_z^2) for r in rs_x]
        @test all(isapprox.(Ix, Ix_ana; rtol = 1e-3))

        # Z-axis profile → w0_y direction (s2=[0,0,-1])
        rs_z = LinRange(-3wy_z, 3wy_z, 51)
        Iz = [BMO.intensity(beam, [0, 0, r], z_eval) for r in rs_z]
        Iz_ana = [I0 * exp(-2 * r^2 / wy_z^2) for r in rs_z]
        @test all(isapprox.(Iz, Iz_ana; rtol = 1e-3))
    end

    @testset "Energy conservation during propagation" begin
        λ = 633e-9
        w0 = 1e-4
        beam = AstigmaticGaussianBeamlet(
            [0.0, 0, 0], [0, 1, 0], λ, w0;
            E0 = [0, 0, 1], support = [1, 0, 0])

        # I_peak(z) * |area(z)| should be constant (energy conservation)
        function peak_times_area(beam, z)
            I_peak = BMO.intensity(beam, [0, 0, 0], z)
            p0, i = BMO.point_on_beam(beam, z)
            h1, _, h2, _, _ = BMO.parabasal_ray_parameters(beam, p0, i)
            dir = BMO.direction(BMO.rays(beam.c)[i])
            area = abs(BMO._pseudo_cross2d(h1, h2, dir))
            return I_peak * area
        end

        P0 = peak_times_area(beam, 0.0)
        P1 = peak_times_area(beam, 0.1)
        P2 = peak_times_area(beam, 0.5)
        @test isapprox(P0, P1; rtol = 1e-6)
        @test isapprox(P0, P2; rtol = 1e-6)
    end

    @testset "Optical power invariant" begin
        λ = 633e-9
        w0_x = 1e-4
        w0_y = 3e-4
        beam = AstigmaticGaussianBeamlet(
            [0.0, 0, 0], [0, 1, 0], λ, w0_x, w0_y;
            E0 = [0, 0, 1], support = [1, 0, 0])
        P = BMO.optical_power(beam)
        @test P > 0
        @test isfinite(P)
    end

    @testset "Peak intensity scaling with beam area" begin
        λ = 633e-9
        w0_x = 1e-4
        w0_y = 2e-4
        beam = AstigmaticGaussianBeamlet(
            [0.0, 0, 0], [0, 1, 0], λ, w0_x, w0_y;
            E0 = [0, 0, 1], support = [1, 0, 0])

        I0 = BMO.intensity(beam, [0, 0, 0], 0.0)
        p0_0, i0 = BMO.point_on_beam(beam, 0.0)
        h1_0, _, h2_0, _, _ = BMO.parabasal_ray_parameters(beam, p0_0, i0)
        dir0 = BMO.direction(BMO.rays(beam.c)[i0])
        area0 = abs(BMO._pseudo_cross2d(h1_0, h2_0, dir0))

        for z in [0.02, 0.05, 0.1, 0.2]
            Iz = BMO.intensity(beam, [0, 0, 0], z)
            p0_z, iz = BMO.point_on_beam(beam, z)
            h1_z, _, h2_z, _, _ = BMO.parabasal_ray_parameters(beam, p0_z, iz)
            dir_z = BMO.direction(BMO.rays(beam.c)[iz])
            area_z = abs(BMO._pseudo_cross2d(h1_z, h2_z, dir_z))
            expected_ratio = area0 / area_z
            @test isapprox(Iz / I0, expected_ratio; rtol = 1e-4)
        end
    end

    @testset "Astigmatic Gouy phase" begin
        λ = 633e-9
        w0_x = 1e-4
        w0_y = 2e-4
        beam = AstigmaticGaussianBeamlet(
            [0.0, 0, 0], [0, 1, 0], λ, w0_x, w0_y;
            E0 = [0, 0, 1], support = [1, 0, 0])

        zRx = π * w0_x^2 / λ
        zRy = π * w0_y^2 / λ
        k = 2π / λ

        for z_eval in [0.01, 0.03, 0.05]
            E = BMO.electric_field(beam, [0, 0, 0], z_eval)
            actual_phase = angle(E * exp(-im * k * z_eval))
            # AGB now matches standard convention → negative Gouy phase
            expected_gouy = -(0.5 * atan(z_eval / zRx) + 0.5 * atan(z_eval / zRy))
            @test isapprox(actual_phase, expected_gouy; atol = 1e-3)
        end
    end

    @testset "Transverse phase curvature" begin
        λ = 633e-9
        w0_x = 1e-4
        w0_y = 2e-4
        beam = AstigmaticGaussianBeamlet(
            [0.0, 0, 0], [0, 1, 0], λ, w0_x, w0_y;
            E0 = [0, 0, 1], support = [1, 0, 0])

        z_eval = 0.1
        k = 2π / λ
        E_center = BMO.electric_field(beam, [0, 0, 0], z_eval)

        # X direction → w0_x → zRx
        x_off = 1e-5
        E_x = BMO.electric_field(beam, [x_off, 0, 0], z_eval)
        Δϕ_x = angle(E_x / E_center)
        zRx = π * w0_x^2 / λ
        Rx = z_eval * (1 + (zRx / z_eval)^2)
        expected_Δϕ_x = k * x_off^2 / (2 * Rx)
        @test isapprox(Δϕ_x, expected_Δϕ_x; rtol = 5e-2)

        # Z direction → w0_y → zRy
        z_off = 1e-5
        E_z = BMO.electric_field(beam, [0, 0, z_off], z_eval)
        Δϕ_z = angle(E_z / E_center)
        zRy = π * w0_y^2 / λ
        Ry = z_eval * (1 + (zRy / z_eval)^2)
        expected_Δϕ_z = k * z_off^2 / (2 * Ry)
        @test isapprox(Δϕ_z, expected_Δϕ_z; rtol = 5e-2)
    end

    # ── Polarization ──

    @testset "Polarization preserved in free-space" begin
        λ = 633e-9
        w0 = 1e-4
        beam = AstigmaticGaussianBeamlet(
            [0.0, 0, 0], [0, 1, 0], λ, w0;
            E0 = [0, 0, 1], support = [1, 0, 0])
        chief = first(BMO.rays(beam.c))
        E = BMO.polarization(chief)
        @test isapprox(normalize(real.(E)), [0, 0, 1]; atol = 1e-10)
    end

    @testset "Polarization rotated by mirror reflection" begin
        λ = 633e-9
        w0 = 1e-4
        m = BMO.SquarePlanoMirror(0.05, 0.01)
        BMO.translate3d!(m, [0, 0.05, 0])
        system = BMO.System([m])
        beam = AstigmaticGaussianBeamlet(
            [0.0, 0, 0], [0, 1, 0], λ, w0;
            E0 = [0, 0, 1], support = [1, 0, 0])
        solve_system!(system, beam)
        reflected_chief = last(BMO.rays(beam.c))
        @test isapprox(BMO.direction(reflected_chief)[2], -1.0; atol = 0.01)
        E = BMO.polarization(reflected_chief)
        d = BMO.direction(reflected_chief)
        @test abs(dot(real.(E), d)) < 1e-10
    end

    @testset "Auto-E0 orthogonality for tilted beams" begin
        for dir in [[1, 1, 0], [0, 1, 1], [1, 1, 1], [1, 0, 0]]
            beam = AstigmaticGaussianBeamlet(
                [0.0, 0, 0], Float64.(dir), 633e-9, 1e-4)
            chief = first(BMO.rays(beam.c))
            E = BMO.polarization(chief)
            d = BMO.direction(chief)
            @test abs(dot(real.(E), d)) < 1e-10
            # Ensure magnitude corresponds to the default P0 = 1e-3
            @test isapprox(BMO.optical_power(beam), 1e-3; atol = 1e-10)
        end
    end

    @testset "Detector integrated power" begin
        λ = 633e-9
        w0 = 0.5mm
        pd = BMO.Detector(10mm)
        BMO.translate3d!(pd, [0, 0.05, 0])
        system = BMO.System([pd])
        beam = AstigmaticGaussianBeamlet(
            [0.0, 0, 0], [0, 1, 0], λ, w0;
            E0 = [0, 0, 1], support = [1, 0, 0])
        solve_system!(system, beam)

        xs, zs, E = BMO.electric_field(pd; n = 50)
        @test size(E) == (50, 50)
        @test any(!iszero, E)

        I = abs2.(E)
        max_idx = argmax(I)
        mid = size(I, 1) ÷ 2 + 1
        @test abs(max_idx[1] - mid) < 5
        @test abs(max_idx[2] - mid) < 5
    end

    @testset "Detector power self-consistency" begin
        # Detector power at two different distances should scale with beam area
        λ = 633e-9
        w0 = 0.5mm
        pd_size = 10mm

        pd1 = BMO.Detector(pd_size)
        BMO.translate3d!(pd1, [0, 0.05, 0])
        s1 = BMO.System([pd1])
        b1 = AstigmaticGaussianBeamlet(
            [0.0, 0, 0], [0, 1, 0], λ, w0;
            E0 = [0, 0, 1], support = [1, 0, 0])
        solve_system!(s1, b1)
        P1 = BMO.optical_power(pd1;
            n = 80, x_min = -pd_size / 2, x_max = pd_size / 2,
            z_min = -pd_size / 2, z_max = pd_size / 2)

        pd2 = BMO.Detector(pd_size)
        BMO.translate3d!(pd2, [0, 0.10, 0])
        s2 = BMO.System([pd2])
        b2 = AstigmaticGaussianBeamlet(
            [0.0, 0, 0], [0, 1, 0], λ, w0;
            E0 = [0, 0, 1], support = [1, 0, 0])
        solve_system!(s2, b2)
        P2 = BMO.optical_power(pd2;
            n = 80, x_min = -pd_size / 2, x_max = pd_size / 2,
            z_min = -pd_size / 2, z_max = pd_size / 2)

        # Power should be conserved (both detectors capture full beam)
        @test isapprox(P1, P2; rtol = 0.05)
    end

    @testset "Beamsplitter power conservation" begin
        λ = 633e-9
        w0 = 0.5mm
        bs = BMO.ThinBeamsplitter(10mm; reflectance = 0.5)
        BMO.zrotate3d!(bs, deg2rad(45))

        pd_t = BMO.Detector(10mm)
        pd_r = BMO.Detector(10mm)
        BMO.translate3d!(pd_t, [0, 0.1, 0])
        BMO.zrotate3d!(pd_t, deg2rad(180))
        BMO.translate3d!(pd_r, [0.1, 0, 0])
        BMO.zrotate3d!(pd_r, deg2rad(90))

        # Reference: single beam straight to detector (no BS)
        pd_ref = BMO.Detector(10mm)
        BMO.translate3d!(pd_ref, [0, 0.1, 0])
        BMO.zrotate3d!(pd_ref, deg2rad(180))

        pd_size = 10mm
        kw = (n = 80, x_min = -pd_size / 2, x_max = pd_size / 2,
            z_min = -pd_size / 2, z_max = pd_size / 2)

        # Measure reference power (no BS)
        sys_ref = BMO.System([pd_ref])
        b_ref = AstigmaticGaussianBeamlet(
            [0.0, -0.1, 0], [0, 1, 0], λ, w0;
            E0 = [0, 0, 1], support = [1, 0, 0])
        solve_system!(sys_ref, b_ref)
        P_ref = BMO.optical_power(pd_ref; kw...)

        # Measure split power
        system = BMO.System([bs, pd_t, pd_r])
        beam = AstigmaticGaussianBeamlet(
            [0.0, -0.1, 0], [0, 1, 0], λ, w0;
            E0 = [0, 0, 1], support = [1, 0, 0])
        solve_system!(system, beam)
        @test length(beam.children) == 2

        P_t = BMO.optical_power(pd_t; kw...)
        P_r = BMO.optical_power(pd_r; kw...)

        # P_t + P_r ≈ P_ref (lossless beamsplitter)
        @test isapprox(P_t + P_r, P_ref; rtol = 0.15)
    end

    @testset "Michelson fringe contrast" begin
        λ = 633e-9
        w0 = 0.5mm
        l0 = 0.1

        m1 = BMO.SquarePlanoMirror2D(BMO.inch)
        m2 = BMO.SquarePlanoMirror2D(BMO.inch)
        bs = BMO.ThinBeamsplitter(BMO.inch; reflectance = 0.5)
        pd = BMO.Detector(10mm)

        BMO.translate3d!(m1, [l0, 0, 0])
        BMO.translate3d!(m2, [0, l0, 0])
        BMO.translate3d!(pd, [-l0, 0, 0])
        BMO.zrotate3d!(bs, deg2rad(45))
        BMO.zrotate3d!(m1, deg2rad(90))
        BMO.zrotate3d!(pd, deg2rad(90))

        system = BMO.System([m1, m2, bs, pd])

        pd_size = 10mm
        Npts = 9
        offsets = LinRange(0, λ, Npts)
        powers = zeros(Npts)

        for (i, δ) in enumerate(offsets)
            BMO.translate_to3d!(m2, [0, l0 + δ, 0])
            beam = AstigmaticGaussianBeamlet(
                [0.0, -l0, 0], [0, 1, 0], λ, w0;
                E0 = [0, 0, 1], support = [1, 0, 0])
            BMO.empty!(pd)
            solve_system!(system, beam)
            powers[i] = BMO.optical_power(pd;
                n = 60, x_min = -pd_size / 2, x_max = pd_size / 2,
                z_min = -pd_size / 2, z_max = pd_size / 2)
        end

        # Signal periodic with period λ
        @test isapprox(powers[1], powers[end]; rtol = 0.1)

        # Fringe contrast > 50%
        P_max = maximum(powers)
        P_min = minimum(powers)
        contrast = (P_max - P_min) / (P_max + P_min)
        @test contrast > 0.5
    end

    @testset "Field symmetry" begin
        λ = 633e-9
        w0 = 1e-4
        beam = AstigmaticGaussianBeamlet(
            [0.0, 0, 0], [0, 1, 0], λ, w0;
            E0 = [0, 0, 1], support = [1, 0, 0])

        z_eval = 0.05
        for r in [1e-5, 3e-5, 5e-5]
            Ep_x = BMO.electric_field(beam, [r, 0, 0], z_eval)
            Em_x = BMO.electric_field(beam, [-r, 0, 0], z_eval)
            @test isapprox(Ep_x, Em_x; rtol = 1e-10)

            Ep_z = BMO.electric_field(beam, [0, 0, r], z_eval)
            Em_z = BMO.electric_field(beam, [0, 0, -r], z_eval)
            @test isapprox(Ep_z, Em_z; rtol = 1e-10)
        end
    end

    @testset "Stigmatic AGB matches standard GaussianBeamlet" begin
        λ = 1000e-9
        w0 = 1e-4
        # Both use same support direction
        agb = AstigmaticGaussianBeamlet(
            [0.0, 0, 0], [0, 1, 0], λ, w0;
            E0 = [0, 0, 1], support = [1, 0, 0])
        gb = GaussianBeamlet(
            [0.0, 0, 0], [0, 1, 0], λ, w0;
            support = [1, 0, 0])

        z_eval = 0.08
        rs = LinRange(-3e-3, 3e-3, 51)

        E_agb = [BMO.electric_field(agb, [r, 0, 0], z_eval) for r in rs]
        E_gb = [BMO.electric_field(gb, abs(r), z_eval) for r in rs]

        # Compare intensity profiles (normalized to peak)
        I_agb = abs2.(E_agb)
        I_gb = abs2.(E_gb)
        I_agb_norm = I_agb / maximum(I_agb)
        I_gb_norm = I_gb / maximum(I_gb)
        @test all(isapprox.(I_agb_norm, I_gb_norm; rtol = 1e-3))
    end

    @testset "Invariant Violation Test" begin
        λ = 532nm
        w0 = 0.1mm

        # Create a highly non-paraxial setup:
        # A flat plate tilted at an extreme angle (80 degrees)
        surf1 = BMO.CircularFlatSurface(10.0mm)
        surf2 = BMO.CircularFlatSurface(10.0mm)

        # n=2.0 (high index)
        lens = BMO.Lens(surf1, surf2, 1.0mm, λ -> 2.0)

        # Tilt the lens by 80 degrees
        BMO.zrotate3d!(lens, deg2rad(80))
        BMO.translate3d!(lens, [0, 1.0mm, 0])

        # Set threshold very low via Ref
        old_threshold = BeamletOptics.INVARIANT_THRESHOLD[]
        BeamletOptics.INVARIANT_THRESHOLD[] = 1e-12

        system = BMO.System([lens])

        # Create a beamlet pointing at the lens
        beam = AstigmaticGaussianBeamlet(
            [0.0, 0, 0], [0, 1, 0], λ, 0.05mm; support = [1, 0, 0])

        # Reset beam for the actual test
        beam = AstigmaticGaussianBeamlet(
            [0.4mm, 0, 0], [0, 1, 0], λ, 0.05mm; support = [1, 0, 0])

        @test_logs (:warn, r"Lagrange invariant violation") match_mode=:any BMO.solve_system!(
            system, beam)

        # Verify that the beam stopped tracing once the invariant failed
        # (solve_system calls trace_system! which has a break on invariant failure)
        @test length(BMO.rays(beam.c)) < 10 # It should stop early

        # Reset threshold
        BMO.INVARIANT_THRESHOLD[] = old_threshold
    end
end # outer testset

end # MODULE
