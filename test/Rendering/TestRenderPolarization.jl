module TestRenderPolarization

using BeamletOptics
using Makie
using Test
using LinearAlgebra

const BMO = BeamletOptics

@testset "Polarization rendering" begin
    Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)
    @test !isnothing(Ext)

    # Shared fixture: 5 mm collimated AGB through an f = 50 mm ThinLens -- the
    # case that exposes the unbounded focus gain (afb6e18).
    agb = AstigmaticGaussianBeamlet([0, 0, 0], [0, 1, 0], 1000e-9, 5e-3; support = [0, 0, 1])
    l = ThinLens(50e-3, 50e-3, BMO.inch, 1.5)
    translate3d!(l, [0, 50e-3, 0])
    solve_system!(System([l]), agb; check_invariant = false)
    flen = 0.12

    (_, b_ref, c_ref) = BMO.waist_parameters(agb, 0.0)
    r_ref = (norm(b_ref) + norm(c_ref)) / 2

    finite_pts(pts) = filter(p -> all(isfinite, p), pts)
    transverse(p) = hypot(p[1], p[3])

    @testset "Amplitude bound (regression)" begin
        gain_max = 10.0
        scale = 1.0
        pts = Ext._polarization_points(agb; flen, scale, gain_max)
        fpts = finite_pts(pts)
        @test !isempty(fpts)
        bound = scale * r_ref * gain_max * (1 + 1e-6)
        @test maximum(transverse, fpts) <= bound
    end

    @testset "Scene sanity" begin
        pts = Ext._polarization_points(agb; flen)
        fpts = finite_pts(pts)
        @test maximum(transverse, fpts) < length(agb) + flen
    end

    @testset "focus_exponent = 0 unchanged" begin
        scale = 1.0
        pts = Ext._polarization_points(agb; flen, scale, focus_exponent = 0.0)
        fpts = finite_pts(pts)
        @test maximum(transverse, fpts) ≈ scale * r_ref atol = 1e-9
    end

    @testset "gain_max is honoured" begin
        scale = 1.0
        pts1 = Ext._polarization_points(agb; flen, scale, gain_max = 1.0)
        pts10 = Ext._polarization_points(agb; flen, scale, gain_max = 10.0)
        pts100 = Ext._polarization_points(agb; flen, scale, gain_max = 100.0)

        amp1 = maximum(transverse, finite_pts(pts1))
        amp10 = maximum(transverse, finite_pts(pts10))
        amp100 = maximum(transverse, finite_pts(pts100))

        @test amp1 ≈ scale * r_ref atol = 1e-9
        @test amp100 > amp10
    end

    @testset "Transversality" begin
        # Ray travels along +x, so the ray's own x-coordinate legitimately grows
        # with arc length; what must stay zero is the field-modulation term's
        # x-component (E0 is already orthogonal to the +x direction). Reconstruct
        # the expected on-axis x-position (matching the uniform t-sampling in
        # `_field_segment!`) and check the residual instead of the raw x value.
        # Note: pts are `Point3f` (Float32), so the tolerance is set accordingly
        # rather than the 1e-12 that would apply to a Float64 computation.
        E0 = [0, 0, 1.0]
        ray = PolarizedRay([0.0, 0, 0], [1.0, 0, 0], 1000e-9, E0)
        beam = Beam(ray)
        flen_t = 0.1
        pts = Ext._polarization_points(beam; flen = flen_t)
        fpts = finite_pts(pts)
        @test !isempty(fpts)
        N = length(fpts)
        ts = range(0, flen_t; length = N)
        disp_x = [p[1] - t for (p, t) in zip(fpts, ts)]
        @test maximum(abs, disp_x) < 1e-5
    end

    @testset "Narrowed signature" begin
        @test_throws MethodError Ext._polarization_points(Beam(Ray([0.0, 0, 0], [1.0, 0, 0])); flen = 0.1)
    end

    @testset "User-facing guard still friendly" begin
        @test_throws ArgumentError Ext._check_polarized(Beam(Ray([0.0, 0, 0], [1.0, 0, 0])))
        @test_throws ArgumentError Ext._check_polarized(Ray([0.0, 0, 0], [1.0, 0, 0]))
    end

    @testset "pol_λ guard" begin
        local pts
        @test_logs (:warn,) match_mode = :any begin
            pts = Ext._polarization_points(agb; flen, λ_vis = 1e-9)
        end
        @test length(pts) < 1e5
    end
end

end # module
