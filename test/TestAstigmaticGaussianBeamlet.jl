module TestAstigmaticGaussianBeamlet

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics

@testset "AstigmaticGaussianBeamlet" begin
    λ0 = 1000e-9
    w0 = 1e-4

    @testset "Construction" begin
        beam = AstigmaticGaussianBeamlet([0.0, 0, 0], [0, 1, 0], λ0, w0; E0=[0,0,1], support=[0,0,1])
        @test isa(beam, AstigmaticGaussianBeamlet{Float64})
        @test BMO.wavelength(beam) == λ0
        @test BMO.direction(beam) ≈ [0, 1, 0]
        @test length(BMO.rays(beam.c)) == 1
        @test length(BMO.rays(beam.wxp)) == 1
        @test length(BMO.rays(beam.dyp)) == 1
        @test isnothing(beam.parent)
        @test isempty(beam.children)
    end

    @testset "Component beam consistency" begin
        beam = AstigmaticGaussianBeamlet([0.0, 0, 0], [0, 1, 0], λ0, w0; E0=[0,0,1], support=[0,0,1])
        all_beams = BMO._component_beams(beam)
        @test length(all_beams) == 9
        aux_beams = BMO._aux_beams(beam)
        @test length(aux_beams) == 8
        for b in all_beams
            @test length(BMO.rays(b)) == 1
        end
    end

    @testset "Free-space field matches analytic Gaussian" begin
        beam = AstigmaticGaussianBeamlet([0.0, 0, 0], [0, 1, 0], λ0, w0; E0=[0,0,1], support=[0,0,1])
        dir = normalize([0, 1, 0])
        ex = normalize(cross(dir, [0,0,1]))
        z_eval = 0.08

        rs = LinRange(-5e-3, 5e-3, 201)
        E_num = [BMO.parabasal_field(beam, ex * r, z_eval) for r in rs]
        E_ana = BMO.electric_field.(rs, z_eval, 1.0, w0, λ0)

        mid = (length(rs) + 1) ÷ 2
        phase = angle(E_num[mid] / E_ana[mid])
        E_ana_aligned = E_ana .* exp(im * phase)

        @test all(isapprox.(E_num, E_ana_aligned; atol=1e-5, rtol=1e-3))
    end

    @testset "Waist parameters" begin
        beam = AstigmaticGaussianBeamlet([0.0, 0, 0], [0, 1, 0], λ0, w0; E0=[0,0,1], support=[0,0,1])
        p0, w1, w2 = BMO.waist_parameters(beam, 0.0)
        @test p0 ≈ [0, 0, 0]
        # At the waist, the beam size should be ~w0
        @test norm(w1) ≈ w0 atol=1e-6
        @test norm(w2) ≈ w0 atol=1e-6
    end

    @testset "Rayleigh range" begin
        beam = AstigmaticGaussianBeamlet([0.0, 0, 0], [0, 1, 0], λ0, w0; E0=[0,0,1], support=[0,0,1])
        rx, ry = BMO.rayleigh_range(beam)
        expected_zr = π * w0^2 / λ0
        @test rx ≈ expected_zr atol=1e-6
        @test ry ≈ expected_zr atol=1e-6
        
        rx2, ry2 = BMO.rayleigh_range(beam; M2=2.0)
        @test rx2 ≈ expected_zr / 2.0 atol=1e-6
        @test ry2 ≈ expected_zr / 2.0 atol=1e-6
    end

    @testset "Batch and In-place Field Evaluation" begin
        beam = AstigmaticGaussianBeamlet([0.0, 0, 0], [0, 1, 0], λ0, w0; E0=[0,0,1], support=[0,0,1])
        dir = normalize([0, 1, 0])
        ex = normalize(cross(dir, [0,0,1]))
        z_eval = 0.05
        rs = [ex * r for r in LinRange(-2e-4, 2e-4, 11)]

        E_single = [BMO.parabasal_field(beam, r, z_eval) for r in rs]
        E_batch = BMO.parabasal_field(beam, rs, z_eval)
        @test E_batch ≈ E_single

        E_out = zeros(ComplexF64, size(rs))
        BMO.electric_field!(E_out, beam, rs, z_eval)
        @test E_out ≈ E_single
    end

    @testset "Waist parameters with offsets" begin
        z0_x = 0.05
        z0_y = 0.10
        beam = AstigmaticGaussianBeamlet(
            [0.0, 0, 0], [0, 1, 0], λ0, w0, 2*w0;
            z0_x = z0_x, z0_y = z0_y,
            E0=[0,0,1], support=[1,0,0]
        )
        
        # X waist should be at z0_x
        _, wx_waist, _ = BMO.waist_parameters(beam, z0_x)
        @test norm(wx_waist) ≈ w0 atol=1e-6
        
        # Y waist should be at z0_y
        _, _, wy_waist = BMO.waist_parameters(beam, z0_y)
        @test norm(wy_waist) ≈ 2*w0 atol=1e-6
        
        # Verify invariants hold
        @test BMO.check_optical_invariant(beam, 1; threshold=1e-12)
    end

    @testset "Lens interaction" begin
        s1 = BMO.CylinderSDF(BMO.inch / 2, BMO.inch)
        l1 = BMO.Lens(s1, n -> 1.5)
        BMO.translate3d!(l1, [0, 0.1, 0.0])
        BMO.zrotate3d!(l1, deg2rad(90))

        system = BMO.System([l1])
        beam = AstigmaticGaussianBeamlet([0.0, 0, 0], [0, 1, 0], λ0, w0; E0=[0,0,1], support=[0,0,1])
        solve_system!(system, beam)

        # Should have 3 segments: initial → enter lens → exit lens → free
        n_rays = length(BMO.rays(beam.c))
        @test n_rays == 3
        for b in BMO._aux_beams(beam)
            @test length(BMO.rays(b)) == n_rays
        end
    end

    @testset "Mirror interaction" begin
        m = BMO.SquarePlanoMirror(0.05, 0.01)
        BMO.translate3d!(m, [0, 0.05, 0])
        system = BMO.System([m])
        beam = AstigmaticGaussianBeamlet([0.0, 0, 0], [0, 1, 0], λ0, w0; E0=[0,0,1], support=[0,0,1])
        solve_system!(system, beam)

        # Should have 2 segments: initial → reflected
        n_rays = length(BMO.rays(beam.c))
        @test n_rays == 2
        for b in BMO._aux_beams(beam)
            @test length(BMO.rays(b)) == n_rays
        end
        # Check reflected direction is approximately [0, -1, 0]
        last_chief = last(BMO.rays(beam.c))
        @test isapprox(BMO.direction(last_chief)[2], -1.0, atol=0.01)
    end

    @testset "ThinBeamsplitter interaction" begin
        bs = BMO.ThinBeamsplitter(0.05)
        BMO.translate3d!(bs, [0, 0.03, 0])
        BMO.zrotate3d!(bs, deg2rad(45))

        system = BMO.System([bs])
        beam = AstigmaticGaussianBeamlet([0.0, 0, 0], [0, 1, 0], λ0, w0; E0=[0,0,1], support=[0,0,1])
        solve_system!(system, beam)

        # Should create 2 children (transmitted + reflected)
        @test length(beam.children) == 2
        t_child = beam.children[1]
        r_child = beam.children[2]
        @test isa(t_child, AstigmaticGaussianBeamlet)
        @test isa(r_child, AstigmaticGaussianBeamlet)
        # Both children should have consistent component beams
        for child in beam.children
            n = length(BMO.rays(child.c))
            for b in BMO._aux_beams(child)
                @test length(BMO.rays(b)) == n
            end
        end
    end

    @testset "Two-lens astigmatic system" begin
        s1 = BMO.CylinderSDF(BMO.inch / 2, BMO.inch)
        s2 = BMO.CylinderSDF(BMO.inch / 2, BMO.inch)
        l1 = BMO.Lens(s1, n -> 1.5)
        l2 = BMO.Lens(s2, n -> 1.5)
        BMO.translate3d!(l1, [0, 0.1, 0.0])
        BMO.translate3d!(l2, [0, 0.15, 0.0])
        BMO.zrotate3d!(l1, deg2rad(90))
        BMO.xrotate3d!(l2, deg2rad(60))

        system = BMO.System([l1, l2])
        beam = AstigmaticGaussianBeamlet([0.0, 0, 0], [0, 1, 0], λ0, w0; E0=[0,0,1], support=[0,0,1])
        solve_system!(system, beam)

        # 2 lenses × 2 surfaces + 1 initial = 5 segments
        n_rays = length(BMO.rays(beam.c))
        @test n_rays == 5
        for b in BMO._aux_beams(beam)
            @test length(BMO.rays(b)) == n_rays
        end
    end

    @testset "Re-solving after a kinematic change" begin
        s1 = BMO.CylinderSDF(BMO.inch / 2, BMO.inch)
        l1 = BMO.Lens(s1, n -> 1.5)
        BMO.translate3d!(l1, [0, 0.1, 0.0])
        BMO.zrotate3d!(l1, deg2rad(90))

        system = BMO.System([l1])
        beam = AstigmaticGaussianBeamlet([0.0, 0, 0], [0, 1, 0], λ0, w0; E0=[0,0,1], support=[0,0,1])
        solve_system!(system, beam)

        n_initial = length(BMO.rays(beam.c))

        # Move lens slightly and solve again
        BMO.translate3d!(l1, [0, 0.005, 0])
        solve_system!(system, beam)

        # Should still have same number of segments
        @test length(BMO.rays(beam.c)) == n_initial
        for b in BMO._aux_beams(beam)
            @test length(BMO.rays(b)) == n_initial
        end
    end

    @testset "Detector interaction" begin
        pd = BMO.Detector(0.05)
        BMO.translate3d!(pd, [0, 0.1, 0])
        system = BMO.System([pd])
        beam = AstigmaticGaussianBeamlet([0.0, 0, 0], [0, 1, 0], λ0, w0)
        solve_system!(system, beam)
        
        @test !isnothing(BMO.hits(pd))
        @test length(BMO.hits(pd)) == 1
        @test isa(BMO.hits(pd)[1], BMO.AstigmaticGaussianBeamletHit)
        
        # Test field calculation
        xs, zs, E = BMO.electric_field(pd; n=20)
        @test size(E) == (20, 20)
        @test any(!iszero, E)
    end
end

end # MODULE
