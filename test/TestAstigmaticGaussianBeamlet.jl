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

    @testset "Support vector" begin
        pos = [0.0, 0, 0]
        dir = [0.0, 1, 0]
        wxp_start(beam) = position(first(BMO.rays(beam.wxp)))
        wyp_start(beam) = position(first(BMO.rays(beam.wyp)))
        ref = AstigmaticGaussianBeamlet(pos, dir, λ0, w0, 2w0; support = [1.0, 0, 0])
        @test wxp_start(ref) ≈ [w0, 0, 0]
        @test wyp_start(ref) ≈ [0, 0, -2w0]
        # the length of an orthogonal support vector does not matter
        for scale in (1e-9, 5, 1e9)
            beam = AstigmaticGaussianBeamlet(pos, dir, λ0, w0, 2w0; support = scale .* [1.0, 0, 0])
            @test wxp_start(beam) ≈ wxp_start(ref)
            @test wyp_start(beam) ≈ wyp_start(ref)
        end
        # nor does the length of the direction
        @test wxp_start(AstigmaticGaussianBeamlet(pos, 3 .* dir, λ0, w0, 2w0; support = [2.0, 0, 0])) ≈ wxp_start(ref)
        # deviations within the threshold pass, e.g. a support vector from a rotation
        @test AstigmaticGaussianBeamlet(pos, dir, λ0, w0; support = [1.0, 1e-12, 0]) isa AstigmaticGaussianBeamlet
        @test AstigmaticGaussianBeamlet(pos, dir, λ0, w0; support = 1e9 .* [1.0, 1e-12, 0]) isa AstigmaticGaussianBeamlet
        # not orthogonal
        @test_throws ArgumentError AstigmaticGaussianBeamlet(pos, dir, λ0, w0; support = [1.0, 1, 0])
        @test_throws ArgumentError AstigmaticGaussianBeamlet(pos, dir, λ0, w0; support = [1.0, 1e-6, 0])
        @test_throws ArgumentError AstigmaticGaussianBeamlet(pos, dir, λ0, w0; support = 1e-9 .* [1.0, 1e-6, 0])
        # parallel and anti-parallel
        @test_throws ArgumentError AstigmaticGaussianBeamlet(pos, dir, λ0, w0; support = [0.0, 1, 0])
        @test_throws ArgumentError AstigmaticGaussianBeamlet(pos, dir, λ0, w0; support = [0.0, -2, 0])
        # no direction
        @test_throws ArgumentError AstigmaticGaussianBeamlet(pos, dir, λ0, w0; support = [0.0, 0, 0])
        # the symmetric constructor passes the support vector on
        @test_throws ArgumentError AstigmaticGaussianBeamlet(pos, dir, λ0, w0; support = [1, 1, 0])
        @test wxp_start(AstigmaticGaussianBeamlet(pos, dir, λ0, w0; support = [3, 0, 0])) ≈ [w0, 0, 0]
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

    @testset "Solving again" begin
        s1 = BMO.CylinderSDF(BMO.inch / 2, BMO.inch)
        l1 = BMO.Lens(s1, n -> 1.5)
        BMO.translate3d!(l1, [0, 0.1, 0.0])
        BMO.zrotate3d!(l1, deg2rad(90))

        system = BMO.System([l1])
        beam = AstigmaticGaussianBeamlet([0.0, 0, 0], [0, 1, 0], λ0, w0; E0=[0,0,1], support=[0,0,1])
        solve_system!(system, beam)

        n_initial = length(BMO.rays(beam.c))

        # Move lens slightly
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

@testset "Gouy phase through a focus" begin
    # The beamlet starts D before its waist, so between the start and points behind the
    # focus its Gouy phase changes by more than π/2. Both field paths must follow the
    # analytic -atan(z/zR) continuously (the principal branch of √(a_ref/a) jumped by π).
    # The parabasal rays are paraxial, which leaves O(θ²) ≈ 5e-4 rad at θ = λ/(π w0).
    λ, w0, D = 1e-6, 20e-6, 5e-3
    k, zR = 2π / λ, π * w0^2 / λ
    wrap(x) = mod2pi(x + π) - π
    gouy(s) = -atan((s - D) / zR) - atan(D / zR)          # relative to the start
    for s in (D - zR, D, D + 0.5zR, D + 2zR, D + 5zR)
        beam() = AstigmaticGaussianBeamlet([0.0, -D, 0], [0.0, 1, 0], λ, w0, w0;
            z0_x = D, z0_y = D, support = [1.0, 0, 0])
        E = BMO.parabasal_field(beam(), zeros(3), s)
        @test wrap(angle(E[argmax(abs.(E))]) - k * s - gouy(s)) ≈ 0 atol = 2e-3

        det = Detector(0.05)
        BMO.translate3d!(det, [0, s - D, 0])
        solve_system!(System([det]), beam())
        ψ = BMO.beamlet_hit_field(only(BMO.hits(det)), [0.0, s - D, 0])
        @test wrap(angle(ψ) - k * s - gouy(s)) ≈ 0 atol = 2e-3
    end
end

@testset "Beam radius within a medium" begin
    # Inside a glass window the transverse term of the field enters with the wavenumber
    # n k0: the beam keeps its radius (with k0 alone it was √n too wide) and spreads with
    # the reduced distance l/n.
    λ, w0, n, d = 1e-6, 100e-6, 1.5, 0.05
    zR = π * w0^2 / λ
    beam = AstigmaticGaussianBeamlet([0.0, -d, 0], [0.0, 1, 0], λ, w0; support = [1.0, 0, 0])
    solve_system!(System([SphericalLens(Inf, Inf, 20e-3, BMO.inch, λ -> n)]), beam)
    for l in (0.0, 5e-3, 15e-3)                         # depth in the glass
        w = w0 * sqrt(1 + ((d + l / n) / zR)^2)
        field(r) = BMO.parabasal_field(beam, [r, 0.0, 0], d + l + 1e-9)
        @test abs(field(w) / field(0.0)) ≈ exp(-1) rtol = 1e-4
    end
end

end # MODULE
