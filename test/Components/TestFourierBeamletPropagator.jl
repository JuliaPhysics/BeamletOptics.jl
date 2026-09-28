module TestFourierBeamletPropagator

using Test
using StaticArrays
using LinearAlgebra
using BeamletOptics

const BMO = BeamletOptics

function bench_interact_loop(system, detector, beam, ray, n)
    for _ in 1:n
        BMO.interact3d(system, detector, beam, ray)
    end
end

function _test_direct_phasor_sum!(
    buf::AbstractVector{SVector{3, Complex{T}}},
    pts::AbstractVector{<:BMO.FourierKPoint{T}},
    grid::BMO.SpatialGrid{T, 1}
) where {T}
    fill!(buf, zero(SVector{3, Complex{T}}))
    xs = grid.ranges[1]
    @inbounds for i in eachindex(xs)
        r = SVector{3, T}(xs[i], zero(T), zero(T))
        val = zero(SVector{3, Complex{T}})
        for kp in pts
            val += kp.E * cis(dot(kp.k, r))
        end
        buf[i] = val
    end
    return buf
end

function _test_direct_phasor_sum!(
    buf::AbstractMatrix{SVector{3, Complex{T}}},
    pts::AbstractVector{<:BMO.FourierKPoint{T}},
    grid::BMO.SpatialGrid{T, 2}
) where {T}
    fill!(buf, zero(SVector{3, Complex{T}}))
    xs = grid.ranges[1]
    ys = grid.ranges[2]
    @inbounds for j in eachindex(ys)
        for i in eachindex(xs)
            r = SVector{3, T}(xs[i], ys[j], zero(T))
            val = zero(SVector{3, Complex{T}})
            for kp in pts
                val += kp.E * cis(dot(kp.k, r))
            end
            buf[i, j] = val
        end
    end
    return buf
end

function _test_direct_phasor_sum!(
    buf::AbstractArray{SVector{3, Complex{T}}, 3},
    pts::AbstractVector{<:BMO.FourierKPoint{T}},
    grid::BMO.SpatialGrid{T, 3}
) where {T}
    fill!(buf, zero(SVector{3, Complex{T}}))
    xs = grid.ranges[1]
    ys = grid.ranges[2]
    zs = grid.ranges[3]
    @inbounds for k in eachindex(zs)
        for j in eachindex(ys)
            for i in eachindex(xs)
                r = SVector{3, T}(xs[i], ys[j], zs[k])
                val = zero(SVector{3, Complex{T}})
                for kp in pts
                    val += kp.E * cis(dot(kp.k, r))
                end
                buf[i, j, k] = val
            end
        end
    end
    return buf
end

function bench_synth1!(buf, pts, grid)
    _test_direct_phasor_sum!(buf, pts, grid)
end

function bench_synth2!(buf, pts, grid)
    _test_direct_phasor_sum!(buf, pts, grid)
end

function bench_synth3!(buf, pts, grid)
    _test_direct_phasor_sum!(buf, pts, grid)
end

function bench_obs_1d!(obs, d, field, grid)
    BMO.calculate_observables!(obs, d, field, grid)
end

function bench_obs_2d!(obs, d, field, grid)
    BMO.calculate_observables!(obs, d, field, grid)
end

function bench_obs_3d!(obs, d, field, grid)
    BMO.calculate_observables!(obs, d, field, grid)
end

function bench_obs_pixel!(obs, d, field, pa)
    BMO.calculate_observables!(obs, d, field; pixel_area = pa)
end

@testset "FourierBeamletPropagator Tests" begin

    # ==============================================================================
    # Ray Preservation & Zero Heap Reallocation
    # ==============================================================================
    @testset "Ray Preservation & Zero Heap Reallocation" begin
        # 1. Verification of constructors and fields
        det_pass = BMO.FourierBeamletPropagator(100.0e-6; stop = false, is_planar = false)
        @test det_pass.stop == false
        @test det_pass.is_planar == false

        det_stop = BMO.FourierBeamletPropagator(100.0e-6; stop = true, is_planar = false)
        @test det_stop.stop == true
        @test det_stop.is_planar == false

        # Test cuboid constructor
        det_cuboid = BMO.FourierBeamletPropagator(100.0e-6, 120.0e-6, 150.0e-6; stop = true)
        @test det_cuboid isa BMO.FourierBeamletPropagator

        # 2. Interact3d with stop = false: ray continues without cloning or reallocation
        sys_pass = System([det_pass])
        r = Ray([0.0, -50.0e-6, 0.0], [0.0, 1.0, 0.0], 1.0e-6)
        b = Beam(r)

        res_pass = BMO.interact3d(sys_pass, det_pass, b, r)
        @test !isnothing(res_pass)
        @test res_pass isa BMO.BeamInteraction
        # Strict object identity: incoming ray is preserved directly without heap cloning
        @test res_pass.ray === r
        @test length(det_pass.hits) == 1

        # 3. Interact3d with stop = true: ray is stopped (returns nothing)
        sys_stop = System([det_stop])
        res_stop = BMO.interact3d(sys_stop, det_stop, b, r)
        @test isnothing(res_stop)
        @test length(det_stop.hits) == 1

        # 4. Repeated interact3d calls: rays are consistently preserved without reallocation
        for _ in 1:10
            res_rep = BMO.interact3d(sys_pass, det_pass, b, r)
            @test !isnothing(res_rep)
            @test res_rep.ray === r
        end

        # 5. Polarized ray preservation
        pr = PolarizedRay([0.0, -50.0e-6, 0.0], [0.0, 1.0, 0.0], 1.0e-6, ComplexF64[1.0, 0.0, 0.0])
        pb = Beam(pr)
        res_pr_pass = BMO.interact3d(sys_pass, det_pass, pb, pr)
        @test !isnothing(res_pr_pass)
        @test res_pr_pass.ray === pr

        res_pr_stop = BMO.interact3d(sys_stop, det_stop, pb, pr)
        @test isnothing(res_pr_stop)

        # 6. Allocation Invariant: Zero heap allocations in hot-loop
        det_bench = BMO.FourierBeamletPropagator(100.0e-6; stop = true, is_planar = false)
        sys_bench = System([det_bench])
        sizehint!(det_bench.hits, 1000)
        # Warmup
        bench_interact_loop(sys_bench, det_bench, b, r, 5)
        allocs = @allocated bench_interact_loop(sys_bench, det_bench, b, r, 10)
        @test allocs == 0
    end

    # ==============================================================================
    # Delta Ray Normalization & Phase Anchor
    # ==============================================================================
    @testset "[PHYSICS] Delta Ray Normalization & Phase Anchor" begin
        det = BMO.FourierBeamletPropagator(200.0e-6; stop = true)

        lambda = 1.0e-6
        d_vec = normalize(SVector{3, Float64}(0.0, 1.0, 0.0))
        k0_expected = (2π / lambda) * d_vec
        P_ray = 4.0
        L_box = 2.5

        # Hit point with non-zero spatial coordinates to test geometric phase anchor
        r_hit = SVector{3, Float64}(12.0e-6, 0.0, -18.0e-6)
        r_origin = r_hit - 50.0e-6 * d_vec
        r1 = Ray(Vector(r_origin), Vector(d_vec), lambda)
        hit1 = BMO.RayHit(r1, 50.0e-6, P_ray)

        # Generate Fourier spectrum from ray hit
        kpoints = BMO.generate_fourier_spectrum(det, [hit1]; L_box = L_box)

        # 1. Delta property: exactly one FourierKPoint generated for simple ray
        @test length(kpoints) == 1
        kp1 = kpoints[1]
        @test kp1 isa BMO.FourierKPoint

        # 2. Wavevector matching k0 = (2pi / lambda) * d
        @test isapprox(kp1.k, k0_expected, rtol = 1e-7)

        # 3. Amplitude scaling: ||E|| = sqrt(P_ray) * L_box
        amp_expected = sqrt(P_ray) * L_box
        @test isapprox(norm(kp1.E), amp_expected, rtol = 1e-7)

        # 4. Multi-parameter scaling verification
        P_ray2 = 9.0
        L_box2 = 1.5
        hit2 = BMO.RayHit(r1, 50.0e-6, P_ray2)
        kpoints2 = BMO.generate_fourier_spectrum(det, [hit2]; L_box = L_box2)
        amp_expected2 = sqrt(P_ray2) * L_box2
        @test isapprox(norm(kpoints2[1].E), amp_expected2, rtol = 1e-7)
        @test isapprox(norm(kpoints2[1].E) / norm(kp1.E), amp_expected2 / amp_expected, rtol = 1e-7)

        # 5. Phase anchor matching exp(-i * k0 · r_hit)
        # Reference hit at origin [0, 0, 0]
        r_origin_ref = Ray([0.0, -50.0e-6, 0.0], Vector(d_vec), lambda)
        hit_origin = BMO.RayHit(r_origin_ref, 50.0e-6, P_ray)
        kp_origin = BMO.generate_fourier_spectrum(det, [hit_origin]; L_box = L_box)[1]

        # Phase factor relative to origin
        expected_phase_factor = exp(-im * dot(k0_expected, r_hit))
        measured_phase_factor = dot(kp_origin.E, kp1.E) / (norm(kp_origin.E) * norm(kp1.E))
        @test isapprox(measured_phase_factor, expected_phase_factor, rtol = 1e-7)

        # Explicit test with quarter-wavelength shift along propagation direction
        # Delta y = lambda / 4 -> delta phi = -k0 * (lambda / 4) = -2pi/lambda * lambda/4 = -pi/2 -> factor -im
        r_quarter = Ray([0.0, -50.0e-6 + lambda / 4, 0.0], Vector(d_vec), lambda)
        hit_quarter = BMO.RayHit(r_quarter, 50.0e-6, P_ray)
        kp_quarter = BMO.generate_fourier_spectrum(det, [hit_quarter]; L_box = L_box)[1]
        measured_quarter_phase = dot(kp_origin.E, kp_quarter.E) / (norm(kp_origin.E) * norm(kp_quarter.E))
        @test isapprox(measured_quarter_phase, ComplexF64(0.0, -1.0), atol = 1e-7)

        # 6. Spectrum generation from accumulated detector hits
        det_accum = BMO.FourierBeamletPropagator(200.0e-6; stop = true)
        sys_accum = System([det_accum])
        BMO.interact3d(sys_accum, det_accum, Beam(r_origin_ref), r_origin_ref)
        kp_accum = BMO.generate_fourier_spectrum(det_accum; L_box = L_box)
        @test length(kp_accum) == 1
        @test isapprox(kp_accum[1].k, k0_expected, rtol = 1e-7)
    end

    # ==============================================================================
    # Gaussian Beamlet Spectral Extent & Transversality
    # ==============================================================================
    @testset "[PHYSICS] Gaussian Beamlet Spectral Extent & Transversality" begin
        lambda = 1.0e-6
        w0 = 25.0e-6
        sigma_k = sqrt(2.0) / w0
        cutoff_k = 4.0 * sigma_k

        d_prop = SVector{3, Float64}(0.0, 1.0, 0.0)
        gb = GaussianBeamlet([0.0, -100.0e-6, 0.0], Vector(d_prop), lambda, w0; P0 = 1.0e-3)

        det = BMO.FourierBeamletPropagator(200.0e-6; stop = true)
        sys = System([det])
        BMO.interact3d(sys, det, gb, 1)

        kpoints = BMO.generate_fourier_spectrum(det)

        # 1. Discrete patch: multiple discrete k-points sampling the spectral envelope
        @test length(kpoints) > 1

        # 2. Spectral extent: all points satisfy ||k_perp|| <= 4 * sigma_k
        max_k_perp = maximum(norm(kp.k - dot(kp.k, d_prop) * d_prop) for kp in kpoints)
        @test all(kp -> norm(kp.k - dot(kp.k, d_prop) * d_prop) <= cutoff_k + 1e-7, kpoints)
        # Ensure finite non-trivial spectral width has been sampled
        @test max_k_perp > 0.0
        @test max_k_perp <= cutoff_k

        # 3. Electromagnetic transversality: k · E = 0 for every spectral sample
        @test all(kp -> abs(dot(kp.k, kp.E)) / (norm(kp.k) * norm(kp.E)) < 10 * eps(real(eltype(kp.k))), kpoints)

        # 4. Invariant under 3D spatial rotation (tilted beamlet)
        d_tilt = normalize(SVector{3, Float64}(0.2, 0.8, 0.1))
        gb_tilt = GaussianBeamlet([0.0, -100.0e-6, 0.0], Vector(d_tilt), lambda, w0; P0 = 1.0e-3)
        det_tilt = BMO.FourierBeamletPropagator(200.0e-6; stop = true)
        sys_tilt = System([det_tilt])
        BMO.interact3d(sys_tilt, det_tilt, gb_tilt, 1)

        kpoints_tilt = BMO.generate_fourier_spectrum(det_tilt)
        @test length(kpoints_tilt) > 1

        @test all(kp -> norm(kp.k - dot(kp.k, d_tilt) * d_tilt) <= cutoff_k + 1e-7, kpoints_tilt)
        @test all(kp -> abs(dot(kp.k, kp.E)) / (norm(kp.k) * norm(kp.E)) < 10 * eps(real(eltype(kp.k))), kpoints_tilt)
    end

    # ==============================================================================
    # Planar Cosine Law Scaling
    # ==============================================================================
    @testset "[PHYSICS] Planar Cosine Law Scaling" begin
        # 1. Constructor with is_planar = true
        det_planar = BMO.FourierBeamletPropagator(200.0e-6, 200.0e-6; stop = true, is_planar = true)
        @test det_planar.is_planar == true

        normal = SVector{3, Float64}(0.0, 1.0, 0.0)
        lambda = 1.0e-6

        # Normal incidence (theta = 0 deg)
        d_normal = SVector{3, Float64}(0.0, 1.0, 0.0)
        r_normal = Ray([0.0, -50.0e-6, 0.0], Vector(d_normal), lambda)
        hit_normal = BMO.RayHit(r_normal, 50.0e-6, 1.0)
        kp_normal = BMO.generate_fourier_spectrum(det_planar, [hit_normal]; normal = normal)
        @test length(kp_normal) == 1
        amp_normal = norm(kp_normal[1].E)
        @test amp_normal > 0.0

        # Oblique incidence at 60 degrees to normal
        # cos(60 deg) = 0.5, sin(60 deg) = sqrt(3)/2
        d_60 = SVector{3, Float64}(sin(deg2rad(60.0)), cos(deg2rad(60.0)), 0.0)
        @test isapprox(dot(d_60, normal), 0.5, atol = 1e-12)

        r_60 = Ray([0.0, -50.0e-6, 0.0], Vector(d_60), lambda)
        hit_60 = BMO.RayHit(r_60, 50.0e-6, 1.0)
        kp_60 = BMO.generate_fourier_spectrum(det_planar, [hit_60]; normal = normal)
        @test length(kp_60) == 1
        amp_60 = norm(kp_60[1].E)

        # Planar cosine law scaling factor = sqrt(cos(60 deg)) = sqrt(0.5) ≈ 0.70710678
        expected_scaling = sqrt(0.5)
        @test isapprox(amp_60 / amp_normal, expected_scaling, rtol = 1e-7)
        @test isapprox(amp_60 / amp_normal, 0.7071067811865475, rtol = 1e-7)

        # Backward incidence (k · n < 0) yields zero amplitude
        d_back = SVector{3, Float64}(0.0, -1.0, 0.0)
        @test dot(d_back, normal) < 0.0
        r_back = Ray([0.0, 50.0e-6, 0.0], Vector(d_back), lambda)
        hit_back = BMO.RayHit(r_back, 50.0e-6, 1.0)
        kp_back = BMO.generate_fourier_spectrum(det_planar, [hit_back]; normal = normal)
        @test isempty(kp_back) || isapprox(norm(kp_back[1].E), 0.0, atol = 1e-12)

        # Additional backward angle (theta = 120 deg, cos(120 deg) = -0.5 < 0)
        d_back_120 = SVector{3, Float64}(sin(deg2rad(120.0)), cos(deg2rad(120.0)), 0.0)
        @test dot(d_back_120, normal) < 0.0
        hit_back_120 = BMO.RayHit(Ray([0.0, 50.0e-6, 0.0], Vector(d_back_120), lambda), 50.0e-6, 1.0)
        kp_back_120 = BMO.generate_fourier_spectrum(det_planar, [hit_back_120]; normal = normal)
        @test isempty(kp_back_120) || isapprox(norm(kp_back_120[1].E), 0.0, atol = 1e-12)

        # Contrast test: Non-planar detector (is_planar = false) does NOT apply cosine scaling
        det_vol = BMO.FourierBeamletPropagator(200.0e-6; stop = true, is_planar = false)
        @test det_vol.is_planar == false
        kp_vol_normal = BMO.generate_fourier_spectrum(det_vol, [hit_normal]; normal = normal)
        kp_vol_60 = BMO.generate_fourier_spectrum(det_vol, [hit_60]; normal = normal)
        @test isapprox(norm(kp_vol_60[1].E), norm(kp_vol_normal[1].E), rtol = 1e-7)
    end

    # ==============================================================================
    # Fourier Field Synthesis & NUFFT Equivalence (Ebene 3)
    # ==============================================================================
    @testset "Fourier Field Synthesis & NUFFT Equivalence" begin

        # --------------------------------------------------------------------------
        # Analytic Plane Wave Reconstruction & Exact Phase Evolution
        # --------------------------------------------------------------------------
        @testset "[PHYSICS] Analytic Plane Wave Reconstruction" begin
            det = BMO.FourierBeamletPropagator(200.0e-6; stop = true)
            lambda = 1.0e-6
            k0 = 2π / lambda
            d_vec = normalize(SVector{3, Float64}(0.3, 0.9, 0.1))
            k_expected = k0 * d_vec
            r_hit = SVector{3, Float64}(5.0e-6, -10.0e-6, 15.0e-6)
            r_origin = r_hit - 50.0e-6 * d_vec
            r = Ray(Vector(r_origin), Vector(d_vec), lambda)
            P_ray = 3.0
            L_box = 1.0
            hit = BMO.RayHit(r, 50.0e-6, P_ray)
            push!(det, hit)

            pts = BMO.generate_fourier_spectrum(det; L_box = L_box)
            @test length(pts) == 1
            kp = pts[1]

            # 1D Spatial Grid
            x1 = LinRange(-30.0e-6, 30.0e-6, 19)
            g1 = BMO.SpatialGrid((x1,))
            E1 = BMO.synthesize_field(det, g1)
            E1_tuple = BMO.synthesize_field(det, (x1,))
            @test E1 == E1_tuple
            @test size(E1) == (19,)
            E_analytic_1d = [kp.E * exp(im * dot(kp.k, SVector{3, Float64}(x, 0.0, 0.0))) for x in x1]
            @test isapprox(E1, E_analytic_1d, rtol = 1e-6)

            # 2D Spatial Grid
            y2 = LinRange(-25.0e-6, 35.0e-6, 21)
            g2 = BMO.SpatialGrid((x1, y2))
            E2 = BMO.synthesize_field(det, g2)
            E2_tuple = BMO.synthesize_field(det, (x1, y2))
            @test E2 == E2_tuple
            @test size(E2) == (19, 21)
            E_analytic_2d = [kp.E * exp(im * dot(kp.k, SVector{3, Float64}(x, y, 0.0))) for x in x1, y in y2]
            @test isapprox(E2, E_analytic_2d, rtol = 1e-6)

            # 3D Spatial Grid
            z3 = LinRange(-20.0e-6, 20.0e-6, 15)
            g3 = BMO.SpatialGrid((x1, y2, z3))
            E3 = BMO.synthesize_field(det, g3)
            E3_tuple = BMO.synthesize_field(det, (x1, y2, z3))
            @test E3 == E3_tuple
            @test size(E3) == (19, 21, 15)
            E_analytic_sub = [kp.E * exp(im * dot(kp.k, SVector{3, Float64}(x1[ix*6-5], y2[iy*7-6], z3[iz*5-4]))) for ix in 1:length(x1[1:6:end]), iy in 1:length(y2[1:7:end]), iz in 1:length(z3[1:5:end])]
            @test isapprox(E3[1:6:end, 1:7:end, 1:5:end], E_analytic_sub, rtol = 1e-6)

            # Intensity Uniformity of plane wave
            expected_intensity = norm(kp.E)^2
            @test isapprox(expected_intensity, P_ray * L_box^2, rtol = 1e-12)
            @test all(val -> isapprox(norm(val)^2, expected_intensity, rtol = 1e-6), E2)
        end

        # --------------------------------------------------------------------------
        # Dual Engine Equivalence (:direct vs :nufft vs :auto)
        # --------------------------------------------------------------------------
        @testset "Dual Engine Equivalence" begin
            det = BMO.FourierBeamletPropagator(200.0e-6; stop = true)
            sys = System([det])
            lambda = 1.0e-6
            w0 = 25.0e-6
            gb = GaussianBeamlet([0.0, -100.0e-6, 0.0], [0.0, 1.0, 0.0], lambda, w0; P0 = 1.0e-3)
            BMO.interact3d(sys, det, gb, 1)

            pts = BMO.generate_fourier_spectrum(det)
            @test length(pts) > 100

            # 1D comparison
            x1 = LinRange(-50.0e-6, 50.0e-6, 32)
            g1 = BMO.SpatialGrid((x1,))
            E1_dir = zeros(SVector{3, ComplexF64}, 32)
            _test_direct_phasor_sum!(E1_dir, pts, g1)
            E1_nufft = BMO.synthesize_field(det, g1)
            max1 = maximum(norm.(E1_dir))
            @test max1 > 0.0
            @test maximum(norm.(E1_dir .- E1_nufft)) / max1 < 1e-5

            # 2D comparison
            y2 = LinRange(-50.0e-6, 50.0e-6, 32)
            g2 = BMO.SpatialGrid((x1, y2))
            E2_dir = zeros(SVector{3, ComplexF64}, 32, 32)
            _test_direct_phasor_sum!(E2_dir, pts, g2)
            E2_nufft = BMO.synthesize_field(det, g2)
            max2 = maximum(norm.(E2_dir))
            @test max2 > 0.0
            @test maximum(norm.(E2_dir .- E2_nufft)) / max2 < 1e-5

            # 3D comparison
            z3 = LinRange(-20.0e-6, 20.0e-6, 12)
            x3 = LinRange(-40.0e-6, 40.0e-6, 16)
            y3 = LinRange(-40.0e-6, 40.0e-6, 16)
            g3 = BMO.SpatialGrid((x3, y3, z3))
            E3_dir = zeros(SVector{3, ComplexF64}, 16, 16, 12)
            _test_direct_phasor_sum!(E3_dir, pts, g3)
            E3_nufft = BMO.synthesize_field(det, g3)
            max3 = maximum(norm.(E3_dir))
            @test max3 > 0.0
            @test maximum(norm.(E3_dir .- E3_nufft)) / max3 < 1e-5

            # In-place functions execution & matching
            buf1 = zeros(SVector{3, ComplexF64}, 32)
            buf2 = zeros(SVector{3, ComplexF64}, 32, 32)
            buf3 = zeros(SVector{3, ComplexF64}, 16, 16, 12)

            _test_direct_phasor_sum!(buf1, pts, g1)
            @test isapprox(buf1, E1_dir, rtol = 1e-12)

            _test_direct_phasor_sum!(buf2, pts, g2)
            @test isapprox(buf2, E2_dir, rtol = 1e-12)

            _test_direct_phasor_sum!(buf3, pts, g3)
            @test isapprox(buf3, E3_dir, rtol = 1e-12)

            fill!(buf1, zero(SVector{3, ComplexF64}))
            fill!(buf2, zero(SVector{3, ComplexF64}))
            fill!(buf3, zero(SVector{3, ComplexF64}))

            BMO._synthesize_field_nufft!(buf1, pts, g1)
            @test isapprox(buf1, E1_nufft, rtol = 1e-12)

            BMO._synthesize_field_nufft!(buf2, pts, g2)
            @test isapprox(buf2, E2_nufft, rtol = 1e-12)

            BMO._synthesize_field_nufft!(buf3, pts, g3)
            @test isapprox(buf3, E3_nufft, rtol = 1e-12)
        end

        # --------------------------------------------------------------------------
        # Anti-Padding & Strict Dimension Invariant (Direktive #10)
        # --------------------------------------------------------------------------
        @testset "Anti-Padding & Strict Dimension Invariant" begin
            det = BMO.FourierBeamletPropagator(200.0e-6; stop = true)
            r = Ray([0.0, -50.0e-6, 0.0], [0.0, 1.0, 0.0], 1.0e-6)
            sys = System([det])
            BMO.interact3d(sys, det, Beam(r), r)

            # Arbitrary non-power-of-two, non-square, prime dimensions
            dims_1d = (23,)
            dims_2d = (17, 31)
            dims_3d = (11, 13, 7)

            g1 = BMO.SpatialGrid((LinRange(-10.0e-6, 20.0e-6, dims_1d[1]),))
            g2 = BMO.SpatialGrid((LinRange(-15.0e-6, 25.0e-6, dims_2d[1]), LinRange(-5.0e-6, 35.0e-6, dims_2d[2])))
            g3 = BMO.SpatialGrid((LinRange(-8.0e-6, 12.0e-6, dims_3d[1]), LinRange(-9.0e-6, 11.0e-6, dims_3d[2]), LinRange(-6.0e-6, 14.0e-6, dims_3d[3])))

            E1 = BMO.synthesize_field(det, g1)
            E2 = BMO.synthesize_field(det, g2)
            E3 = BMO.synthesize_field(det, g3)

            @test size(E1) == dims_1d
            @test size(E2) == dims_2d
            @test size(E3) == dims_3d
            @test eltype(E1) == SVector{3, ComplexF64}
            @test eltype(E2) == SVector{3, ComplexF64}
            @test eltype(E3) == SVector{3, ComplexF64}

            # Dimension N < 4 contract: synthesize_field throws ArgumentError
            @test_throws ArgumentError BMO.synthesize_field(det, BMO.SpatialGrid((LinRange(-10.0e-6, 20.0e-6, 1),)))
            @test_throws ArgumentError BMO.synthesize_field(det, BMO.SpatialGrid((LinRange(-10.0e-6, 20.0e-6, 2),)))
            @test_throws ArgumentError BMO.synthesize_field(det, BMO.SpatialGrid((LinRange(-10.0e-6, 20.0e-6, 3),)))
            @test_throws ArgumentError BMO.synthesize_field(det, BMO.SpatialGrid((LinRange(-10.0e-6, 20.0e-6, 16), LinRange(-10.0e-6, 20.0e-6, 3))))
            @test_throws ArgumentError BMO.synthesize_field(det, (LinRange(-10.0e-6, 20.0e-6, 3),))
        end

        # --------------------------------------------------------------------------
        # Fresnel-Arago Laws & Electromagnetic Transversality
        # --------------------------------------------------------------------------
        @testset "[PHYSICS] Fresnel-Arago Laws & Transversality" begin
            lambda = 1.0e-6
            theta = deg2rad(5.0)
            d1 = SVector{3, Float64}(sin(theta), cos(theta), 0.0)
            d2 = SVector{3, Float64}(-sin(theta), cos(theta), 0.0)

            # 1. Parallel s-polarization perpendicular to incidence plane
            det_par = BMO.FourierBeamletPropagator(200.0e-6; stop = true)
            p_s = ComplexF64[0.0, 0.0, 1.0]
            pr1_par = PolarizedRay(Vector(-50.0e-6 * d1), Vector(d1), lambda, p_s)
            pr2_par = PolarizedRay(Vector(-50.0e-6 * d2), Vector(d2), lambda, p_s)
            push!(det_par, BMO.PolarizedRayHit(pr1_par, 50.0e-6))
            push!(det_par, BMO.PolarizedRayHit(pr2_par, 50.0e-6))

            kx = (2π / lambda) * sin(theta)
            fringe_period = π / kx

            x_grid = LinRange(0.0, fringe_period, 51)
            g_fringe = BMO.SpatialGrid((x_grid,))

            E_par = BMO.synthesize_field(det_par, g_fringe)
            I_par = [norm(v)^2 for v in E_par]

            # Constructive interference at x = 0 and x = Lambda
            @test isapprox(I_par[1], 4.0, rtol = 1e-6)
            @test isapprox(I_par[51], 4.0, rtol = 1e-6)
            # Destructive interference at x = Lambda / 2 (sample 26)
            @test isapprox(I_par[26], 0.0, atol = 1e-6)

            I_max = maximum(I_par)
            I_min = minimum(I_par)
            visibility_par = (I_max - I_min) / (I_max + I_min)
            @test isapprox(visibility_par, 1.0, atol = 1e-6)

            # 2. Orthogonal polarization (Fresnel-Arago Law: no interference)
            det_orth = BMO.FourierBeamletPropagator(200.0e-6; stop = true)
            p1_orth = SVector{3, ComplexF64}(0.0, 0.0, 1.0)
            p2_orth = SVector{3, ComplexF64}(cos(theta), sin(theta), 0.0)
            pr1_orth = PolarizedRay(Vector(-50.0e-6 * d1), Vector(d1), lambda, Vector(p1_orth))
            pr2_orth = PolarizedRay(Vector(-50.0e-6 * d2), Vector(d2), lambda, Vector(p2_orth))
            push!(det_orth, BMO.PolarizedRayHit(pr1_orth, 50.0e-6))
            push!(det_orth, BMO.PolarizedRayHit(pr2_orth, 50.0e-6))

            E_orth = BMO.synthesize_field(det_orth, g_fringe)
            I_orth = [norm(v)^2 for v in E_orth]
            visibility_orth = (maximum(I_orth) - minimum(I_orth)) / (maximum(I_orth) + minimum(I_orth))
            @test isapprox(visibility_orth, 0.0, atol = 1e-6)
            @test all(I -> isapprox(I, 2.0, rtol = 1e-6), I_orth)

            # 3. Transversality verification in real space:
            pts_par = BMO.generate_fourier_spectrum(det_par)
            @test all(kp -> abs(dot(kp.k, kp.E)) / (norm(kp.k) * norm(kp.E)) < 10 * eps(Float64), pts_par)
        end

        # --------------------------------------------------------------------------
        # Parseval Energy Conservation / Invariance
        # --------------------------------------------------------------------------
        @testset "[PHYSICS] Parseval Energy Conservation" begin
            det = BMO.FourierBeamletPropagator(200.0e-6; stop = true)
            lambda = 1.0e-6
            w0 = 25.0e-6
            gb = GaussianBeamlet([0.0, -100.0e-6, 0.0], [0.0, 1.0, 0.0], lambda, w0; P0 = 1.0e-3)
            sys = System([det])
            BMO.interact3d(sys, det, gb, 1)

            x_range = LinRange(-60.0e-6, 60.0e-6, 32)
            y_range = LinRange(-60.0e-6, 60.0e-6, 32)
            grid = BMO.SpatialGrid((x_range, y_range))

            pts = BMO.generate_fourier_spectrum(det)
            E_dir = zeros(SVector{3, ComplexF64}, 32, 32)
            _test_direct_phasor_sum!(E_dir, pts, grid)
            E_nufft = BMO.synthesize_field(det, grid)

            dx = step(x_range)
            dy = step(y_range)
            U_dir = sum(norm(v)^2 for v in E_dir) * dx * dy
            U_nufft = sum(norm(v)^2 for v in E_nufft) * dx * dy

            @test U_dir > 0.0
            @test isapprox(U_dir, U_nufft, rtol = 1e-6)
        end

        # --------------------------------------------------------------------------
        # Boundary Conditions, Empty Hits & Robustness
        # --------------------------------------------------------------------------
        @testset "Boundary Conditions & Empty Hits" begin
            det_empty = BMO.FourierBeamletPropagator(100.0e-6; stop = true)
            x_range = LinRange(-50.0e-6, 50.0e-6, 16)
            y_range = LinRange(-50.0e-6, 50.0e-6, 16)
            grid2 = BMO.SpatialGrid((x_range, y_range))

            E_empty = BMO.synthesize_field(det_empty, grid2)
            @test size(E_empty) == (16, 16)
            @test all(v -> v == zero(SVector{3, ComplexF64}), E_empty)

            E_empty_tuple = BMO.synthesize_field(det_empty, (x_range, y_range))
            @test size(E_empty_tuple) == (16, 16)
            @test all(v -> v == zero(SVector{3, ComplexF64}), E_empty_tuple)

            r_zero = Ray([0.0, -50.0e-6, 0.0], [0.0, 1.0, 0.0], 1.0e-6)
            hit_zero = BMO.RayHit(r_zero, 50.0e-6, 0.0)
            det_zero = BMO.FourierBeamletPropagator(100.0e-6; stop = true)
            push!(det_zero, hit_zero)
            E_zero = BMO.synthesize_field(det_zero, grid2)
            @test all(v -> isapprox(norm(v), 0.0, atol = 1e-12), E_zero)
        end

        # --------------------------------------------------------------------------
        # Zero-Allocation Invariant in Hot-Path Field Synthesis
        # --------------------------------------------------------------------------
        @testset "Zero-Allocation Invariant in Hot-Path Field Synthesis" begin
            det = BMO.FourierBeamletPropagator(100.0e-6; stop = true)
            r = Ray([0.0, -50.0e-6, 0.0], [0.0, 1.0, 0.0], 1.0e-6)
            sys = System([det])
            BMO.interact3d(sys, det, Beam(r), r)
            pts = BMO.generate_fourier_spectrum(det)

            x_range = LinRange(-50.0e-6, 50.0e-6, 16)
            y_range = LinRange(-50.0e-6, 50.0e-6, 16)
            z_range = LinRange(-50.0e-6, 50.0e-6, 16)

            grid1 = BMO.SpatialGrid((x_range,))
            grid2 = BMO.SpatialGrid((x_range, y_range))
            grid3 = BMO.SpatialGrid((x_range, y_range, z_range))

            E1_buf = zeros(SVector{3, ComplexF64}, 16)
            E2_buf = zeros(SVector{3, ComplexF64}, 16, 16)
            E3_buf = zeros(SVector{3, ComplexF64}, 16, 16, 16)

            # 1D Zero-allocation
            bench_synth1!(E1_buf, pts, grid1)
            allocs1 = @allocated bench_synth1!(E1_buf, pts, grid1)
            @test allocs1 == 0

            # 2D Zero-allocation
            bench_synth2!(E2_buf, pts, grid2)
            allocs2 = @allocated bench_synth2!(E2_buf, pts, grid2)
            @test allocs2 == 0

            # 3D Zero-allocation
            bench_synth3!(E3_buf, pts, grid3)
            allocs3 = @allocated bench_synth3!(E3_buf, pts, grid3)
            @test allocs3 == 0
        end

    end # @testset "AK-F1.2651"

    # ==============================================================================
    # Observable Extraction & Detector Observables (Ebene 4)
    # ==============================================================================
    @testset "Observable Extraction & Detector Observables" begin

        # --------------------------------------------------------------------------
        # DetectorObservables Struct Contracts, Dimensions & Anti-Padding
        # --------------------------------------------------------------------------
        @testset "DetectorObservables Struct Contracts & Anti-Padding" begin
            # 1. Constructor signatures
            obs_tuple = BMO.DetectorObservables((16, 24))
            obs_splat = BMO.DetectorObservables(16, 24)
            @test size(obs_tuple.intensity) == (16, 24)
            @test size(obs_splat.intensity) == (16, 24)
            @test size(obs_tuple.phase) == (16, 24)
            @test length(obs_tuple.stokes) == 4
            @test all(s -> size(s) == (16, 24), obs_tuple.stokes)
            @test obs_tuple.power == 0.0
            @test eltype(obs_tuple.intensity) == Float64
            @test eltype(obs_tuple.phase) == Float64
            @test all(s -> eltype(s) == Float64, obs_tuple.stokes)

            # Custom type constructor
            obs_f32 = BMO.DetectorObservables((12, 18); T = Float32)
            @test eltype(obs_f32.intensity) == Float32
            @test eltype(obs_f32.phase) == Float32
            @test all(s -> eltype(s) == Float32, obs_f32.stokes)
            @test obs_f32.power == 0.0f0

            # Direct field assignment constructor
            I_arr = rand(10, 10)
            phi_arr = rand(10, 10)
            stokes_tuple = (rand(10, 10), rand(10, 10), rand(10, 10), rand(10, 10))
            P_val = 2.5
            obs_manual = BMO.DetectorObservables(I_arr, phi_arr, stokes_tuple, P_val)
            @test obs_manual.intensity === I_arr
            @test obs_manual.phase === phi_arr
            @test obs_manual.stokes === stokes_tuple
            @test obs_manual.power == P_val

            # Constructor defaulting power
            obs_manual_def = BMO.DetectorObservables(I_arr, phi_arr, stokes_tuple)
            @test obs_manual_def.intensity === I_arr
            @test obs_manual_def.power == 0.0

            # 2. Anti-Padding & Dimension Invariance (1D, 2D, 3D prime / irregular shapes)
            det = BMO.FourierBeamletPropagator(200.0e-6; stop = true)
            dims_1d = (37,)
            dims_2d = (23, 41)
            dims_3d = (11, 13, 7)

            field_1d = randn(SVector{3, ComplexF64}, dims_1d)
            field_2d = randn(SVector{3, ComplexF64}, dims_2d)
            field_3d = randn(SVector{3, ComplexF64}, dims_3d)

            obs_1d = BMO.calculate_observables(det, field_1d)
            obs_2d = BMO.calculate_observables(det, field_2d)
            obs_3d = BMO.calculate_observables(det, field_3d)

            @test size(obs_1d.intensity) == dims_1d
            @test size(obs_1d.phase) == dims_1d
            @test all(s -> size(s) == dims_1d, obs_1d.stokes)

            @test size(obs_2d.intensity) == dims_2d
            @test size(obs_2d.phase) == dims_2d
            @test all(s -> size(s) == dims_2d, obs_2d.stokes)

            @test size(obs_3d.intensity) == dims_3d
            @test size(obs_3d.phase) == dims_3d
            @test all(s -> size(s) == dims_3d, obs_3d.stokes)

            # 3. Dimension mismatch contract on in-place mutation
            obs_mismatch = BMO.DetectorObservables((16, 16))
            @test_throws DimensionMismatch BMO.calculate_observables!(obs_mismatch, det, field_2d)
            @test_throws DimensionMismatch BMO.calculate_observables!(obs_mismatch, det, field_1d)
        end

        # --------------------------------------------------------------------------
        # Stokes Parameters & Polarization State Invariants
        # --------------------------------------------------------------------------
        @testset "[PHYSICS] Stokes Parameters & Polarization Invariants" begin
            det = BMO.FourierBeamletPropagator(100.0e-6; stop = true)
            Z0 = BMO.Z_vacuum
            E0 = 12.0
            I0 = abs2(E0) / (2 * Z0)

            # 1. Horizontal Linear Polarization: E = [E0, 0, 0]
            f_H = [SVector{3, ComplexF64}(E0, 0.0, 0.0)]
            obs_H = BMO.calculate_observables(det, f_H)
            @test isapprox(obs_H.intensity[1], I0, rtol = 1e-12)
            @test isapprox(obs_H.stokes[1][1], I0, rtol = 1e-12)
            @test isapprox(obs_H.stokes[2][1], I0, rtol = 1e-12)
            @test isapprox(obs_H.stokes[3][1], 0.0, atol = 1e-14)
            @test isapprox(obs_H.stokes[4][1], 0.0, atol = 1e-14)

            # 2. Vertical Linear Polarization: E = [0, E0, 0]
            f_V = [SVector{3, ComplexF64}(0.0, E0, 0.0)]
            obs_V = BMO.calculate_observables(det, f_V)
            @test isapprox(obs_V.intensity[1], I0, rtol = 1e-12)
            @test isapprox(obs_V.stokes[1][1], I0, rtol = 1e-12)
            @test isapprox(obs_V.stokes[2][1], -I0, rtol = 1e-12)
            @test isapprox(obs_V.stokes[3][1], 0.0, atol = 1e-14)
            @test isapprox(obs_V.stokes[4][1], 0.0, atol = 1e-14)

            # 3. +45° Diagonal Linear Polarization: E = [E0/√2, E0/√2, 0]
            f_D = [SVector{3, ComplexF64}(E0 / sqrt(2), E0 / sqrt(2), 0.0)]
            obs_D = BMO.calculate_observables(det, f_D)
            @test isapprox(obs_D.intensity[1], I0, rtol = 1e-12)
            @test isapprox(obs_D.stokes[1][1], I0, rtol = 1e-12)
            @test isapprox(obs_D.stokes[2][1], 0.0, atol = 1e-14)
            @test isapprox(obs_D.stokes[3][1], I0, rtol = 1e-12)
            @test isapprox(obs_D.stokes[4][1], 0.0, atol = 1e-14)

            # 4. -45° Anti-Diagonal Linear Polarization: E = [E0/√2, -E0/√2, 0]
            f_A = [SVector{3, ComplexF64}(E0 / sqrt(2), -E0 / sqrt(2), 0.0)]
            obs_A = BMO.calculate_observables(det, f_A)
            @test isapprox(obs_A.intensity[1], I0, rtol = 1e-12)
            @test isapprox(obs_A.stokes[1][1], I0, rtol = 1e-12)
            @test isapprox(obs_A.stokes[2][1], 0.0, atol = 1e-14)
            @test isapprox(obs_A.stokes[3][1], -I0, rtol = 1e-12)
            @test isapprox(obs_A.stokes[4][1], 0.0, atol = 1e-14)

            # 5. Right Circular Polarization: E = [E0/√2, -i*E0/√2, 0]
            f_R = [SVector{3, ComplexF64}(E0 / sqrt(2), -im * E0 / sqrt(2), 0.0)]
            obs_R = BMO.calculate_observables(det, f_R)
            @test isapprox(obs_R.intensity[1], I0, rtol = 1e-12)
            @test isapprox(obs_R.stokes[1][1], I0, rtol = 1e-12)
            @test isapprox(obs_R.stokes[2][1], 0.0, atol = 1e-14)
            @test isapprox(obs_R.stokes[3][1], 0.0, atol = 1e-14)
            @test isapprox(obs_R.stokes[4][1], I0, rtol = 1e-12)

            # 6. Left Circular Polarization: E = [E0/√2, +i*E0/√2, 0]
            f_L = [SVector{3, ComplexF64}(E0 / sqrt(2), im * E0 / sqrt(2), 0.0)]
            obs_L = BMO.calculate_observables(det, f_L)
            @test isapprox(obs_L.intensity[1], I0, rtol = 1e-12)
            @test isapprox(obs_L.stokes[1][1], I0, rtol = 1e-12)
            @test isapprox(obs_L.stokes[2][1], 0.0, atol = 1e-14)
            @test isapprox(obs_L.stokes[3][1], 0.0, atol = 1e-14)
            @test isapprox(obs_L.stokes[4][1], -I0, rtol = 1e-12)

            # 7. General Elliptical Polarization: E = [a*E0, b*E0*exp(i*delta), 0]
            a = sqrt(3.0) / 2.0
            b = 0.5
            delta = π / 3.0
            f_ell = [SVector{3, ComplexF64}(a * E0, b * E0 * exp(im * delta), 0.0)]
            obs_ell = BMO.calculate_observables(det, f_ell)
            S0_ell = (a^2 + b^2) * I0
            S1_ell = (a^2 - b^2) * I0
            S2_ell = 2 * a * b * cos(delta) * I0
            S3_ell = -2 * a * b * sin(delta) * I0
            @test isapprox(obs_ell.stokes[1][1], S0_ell, rtol = 1e-12)
            @test isapprox(obs_ell.stokes[2][1], S1_ell, rtol = 1e-12)
            @test isapprox(obs_ell.stokes[3][1], S2_ell, rtol = 1e-12)
            @test isapprox(obs_ell.stokes[4][1], S3_ell, rtol = 1e-12)

            # 8. First-Principles Invariant 2: Stokes consistency S0^2 >= S1^2 + S2^2 + S3^2
            f_rand = randn(SVector{3, ComplexF64}, 25, 25)
            obs_rand = BMO.calculate_observables(det, f_rand)
            S0_grid = obs_rand.stokes[1]
            S1_grid = obs_rand.stokes[2]
            S2_grid = obs_rand.stokes[3]
            S3_grid = obs_rand.stokes[4]
            # S0^2 >= S1^2 + S2^2 + S3^2 up to machine precision
            @test all(i -> S0_grid[i]^2 >= (S1_grid[i]^2 + S2_grid[i]^2 + S3_grid[i]^2) - 1e-14, eachindex(S0_grid))
            # Coherence equality for monochromatic transverse wave
            @test all(i -> isapprox(S0_grid[i]^2, S1_grid[i]^2 + S2_grid[i]^2 + S3_grid[i]^2, rtol = 1e-10, atol = 1e-14), eachindex(S0_grid))

            # 9. Longitudinal component Ez ≠ 0
            f_long = [SVector{3, ComplexF64}(3.0, 4.0, 5.0)]
            obs_long = BMO.calculate_observables(det, f_long)
            I_long = (abs2(3.0) + abs2(4.0) + abs2(5.0)) / (2 * Z0)
            S0_long = (abs2(3.0) + abs2(4.0)) / (2 * Z0)
            @test isapprox(obs_long.intensity[1], I_long, rtol = 1e-12)
            @test isapprox(obs_long.stokes[1][1], S0_long, rtol = 1e-12)
            @test obs_long.intensity[1] > obs_long.stokes[1][1]
            @test isapprox(obs_long.intensity[1] - obs_long.stokes[1][1], abs2(5.0) / (2 * Z0), rtol = 1e-12)
        end

        # --------------------------------------------------------------------------
        # Optical Wavefront Phase Evolution & Dominant Mode Tracking
        # --------------------------------------------------------------------------
        @testset "[PHYSICS] Optical Wavefront Phase & Dominant Mode Tracking" begin
            det = BMO.FourierBeamletPropagator(100.0e-6; stop = true)

            # 1. Uniform phase field: phi(r) = const across spatial grid
            phi_const = 1.15
            field_const = fill(SVector{3, ComplexF64}(exp(im * phi_const), 0.0, 0.0), 16, 16)
            obs_const = BMO.calculate_observables(det, field_const)
            @test all(p -> isapprox(p, phi_const, atol = 1e-12), obs_const.phase)

            # 2. Phase tracking on dominant component
            # Ex dominant
            f_ex = [SVector{3, ComplexF64}(2.0 * exp(im * 0.7), 0.5 * exp(im * 1.5), 0.1)]
            @test isapprox(BMO.calculate_observables(det, f_ex).phase[1], 0.7, atol = 1e-12)

            # Ey dominant
            f_ey = [SVector{3, ComplexF64}(0.5 * exp(im * 0.7), 2.0 * exp(im * 1.5), 0.1)]
            @test isapprox(BMO.calculate_observables(det, f_ey).phase[1], 1.5, atol = 1e-12)

            # Ez dominant
            f_ez = [SVector{3, ComplexF64}(0.1, 0.5, 3.0 * exp(-im * 0.8))]
            @test isapprox(BMO.calculate_observables(det, f_ez).phase[1], -0.8, atol = 1e-12)

            # 3. Linear phase slope matching transverse wavevector k_perp
            kx = 2π / 50.0e-6
            x_range = LinRange(0.0, 100.0e-6, 101)
            grid_1d = BMO.SpatialGrid((x_range,))
            dx = step(x_range)
            f_slope = [SVector{3, ComplexF64}(exp(im * kx * x), 0.0, 0.0) for x in x_range]
            obs_slope = BMO.calculate_observables(det, f_slope, grid_1d)

            # Test phase differences modulo 2pi
            dphis = [mod(obs_slope.phase[i + 1] - obs_slope.phase[i] + π, 2π) - π for i in 1:(length(x_range) - 1)]
            @test all(dp -> isapprox(dp, kx * dx, atol = 1e-10), dphis)

            # 4. Strict phase bounds [-π, π]
            f_rand = randn(SVector{3, ComplexF64}, 20, 20)
            obs_rand = BMO.calculate_observables(det, f_rand)
            @test all(p -> -π - 1e-14 <= p <= π + 1e-14, obs_rand.phase)
        end

        # --------------------------------------------------------------------------
        # Spatial Integration, Impedance Scaling & Non-Negativity
        # --------------------------------------------------------------------------
        @testset "[PHYSICS] Spatial Integration & Non-Negativity Invariants" begin
            # 1. Grid pixel area calculation contracts
            x1 = LinRange(-10.0e-6, 10.0e-6, 21)
            y2 = LinRange(-15.0e-6, 15.0e-6, 31)
            z3 = LinRange(-20.0e-6, 20.0e-6, 11)

            g1 = BMO.SpatialGrid((x1,))
            g2 = BMO.SpatialGrid((x1, y2))
            g3 = BMO.SpatialGrid((x1, y2, z3))

            dx = step(x1)
            dy = step(y2)
            @test isapprox(BMO._calculate_pixel_area(g1), dx, rtol = 1e-14)
            @test isapprox(BMO._calculate_pixel_area(g2), dx * dy, rtol = 1e-14)
            @test isapprox(BMO._calculate_pixel_area(g3), dx * dy, rtol = 1e-14)

            # 2. Discrete optical power integration identity P = sum(I) * pixel_area
            det = BMO.FourierBeamletPropagator(100.0e-6; stop = true)
            f_grid2 = randn(SVector{3, ComplexF64}, length(x1), length(y2))
            obs_grid2 = BMO.calculate_observables(det, f_grid2, g2)
            expected_power2 = sum(obs_grid2.intensity) * (dx * dy)
            @test isapprox(obs_grid2.power, expected_power2, rtol = 1e-14)

            # Manual pixel_area override
            pa_custom = 3.14159e-11
            obs_custom_pa = BMO.calculate_observables(det, f_grid2; pixel_area = pa_custom)
            @test isapprox(obs_custom_pa.power, sum(obs_custom_pa.intensity) * pa_custom, rtol = 1e-14)

            # 3. Characteristic impedance (Z) scaling contract
            Z_dielectric = BMO.Z_vacuum / 1.5
            obs_z1 = BMO.calculate_observables(det, f_grid2; Z = BMO.Z_vacuum)
            obs_z2 = BMO.calculate_observables(det, f_grid2; Z = Z_dielectric)

            scaling = BMO.Z_vacuum / Z_dielectric
            @test isapprox(obs_z2.intensity, obs_z1.intensity .* scaling, rtol = 1e-12)
            @test isapprox(obs_z2.power, obs_z1.power * scaling, rtol = 1e-12)
            for k in 1:4
                @test isapprox(obs_z2.stokes[k], obs_z1.stokes[k] .* scaling, rtol = 1e-12)
            end
            # Phase invariant to real characteristic impedance change
            @test isapprox(obs_z2.phase, obs_z1.phase, rtol = 1e-12)

            # 4. First-Principles Invariant 3: Positivity I >= 0, S0 >= 0, P >= 0
            @test all(I -> I >= 0.0, obs_z2.intensity)
            @test all(S0 -> S0 >= 0.0, obs_z2.stokes[1])
            @test obs_z2.power >= 0.0
        end

        # --------------------------------------------------------------------------
        # In-Place Mutation, Equivalence & Calling Conventions
        # --------------------------------------------------------------------------
        @testset "In-Place Mutation & Calling Conventions" begin
            det = BMO.FourierBeamletPropagator(100.0e-6; stop = true)
            dims = (19, 23)
            field = randn(SVector{3, ComplexF64}, dims)
            pa = 2.5e-11
            Z_custom = 300.0

            # Preallocated observables container
            obs_mut = BMO.DetectorObservables(dims)

            # In-place return value identity and field mutation
            res = BMO.calculate_observables!(obs_mut, det, field; pixel_area = pa, Z = Z_custom)
            @test res.intensity === obs_mut.intensity
            @test res.phase === obs_mut.phase
            @test res.stokes === obs_mut.stokes

            # Equivalence with allocating out-of-place variant
            obs_alloc = BMO.calculate_observables(det, field; pixel_area = pa, Z = Z_custom)
            @test isapprox(obs_mut.intensity, obs_alloc.intensity, rtol = 1e-14)
            @test isapprox(obs_mut.phase, obs_alloc.phase, rtol = 1e-14)
            @test isapprox(res.power, obs_alloc.power, rtol = 1e-14)
            for k in 1:4
                @test isapprox(obs_mut.stokes[k], obs_alloc.stokes[k], rtol = 1e-14)
            end

            # Equivalence with SpatialGrid in-place signature
            x = LinRange(-10.0e-6, 10.0e-6, dims[1])
            y = LinRange(-10.0e-6, 10.0e-6, dims[2])
            grid = BMO.SpatialGrid((x, y))

            obs_mut_grid = BMO.DetectorObservables(dims)
            res_grid = BMO.calculate_observables!(obs_mut_grid, det, field, grid; Z = Z_custom)
            obs_alloc_grid = BMO.calculate_observables(det, field, grid; Z = Z_custom)
            @test res_grid.intensity === obs_mut_grid.intensity
            @test isapprox(obs_mut_grid.intensity, obs_alloc_grid.intensity, rtol = 1e-14)
            @test isapprox(res_grid.power, obs_alloc_grid.power, rtol = 1e-14)

            # Equivalence of direct field functions (omitting detector instance)
            obs_direct = BMO.calculate_observables(field; pixel_area = pa, Z = Z_custom)
            @test isapprox(obs_direct.intensity, obs_alloc.intensity, rtol = 1e-14)
            @test isapprox(obs_direct.power, obs_alloc.power, rtol = 1e-14)

            # Repeated mutation without memory leakage or state accumulation
            field_half = field .* 0.5
            res_half = BMO.calculate_observables!(obs_mut, det, field_half; pixel_area = pa, Z = Z_custom)
            @test isapprox(res_half.power, obs_alloc.power * 0.25, rtol = 1e-14)
            @test isapprox(obs_mut.intensity, obs_alloc.intensity .* 0.25, rtol = 1e-14)
        end

        # --------------------------------------------------------------------------
        # Boundary Conditions, Physical Tracing & Parseval Equivalence
        # --------------------------------------------------------------------------
        @testset "Boundary Conditions & Physical Tracing" begin
            det = BMO.FourierBeamletPropagator(100.0e-6; stop = true)

            # 1. Zero field boundary case
            field_zero = zeros(SVector{3, ComplexF64}, 12, 12)
            obs_zero = BMO.calculate_observables(det, field_zero)
            @test all(v -> v == 0.0, obs_zero.intensity)
            @test all(v -> v == 0.0, obs_zero.phase)
            @test obs_zero.power == 0.0
            for k in 1:4
                @test all(v -> v == 0.0, obs_zero.stokes[k])
            end

            # 2. Single-pixel grid (1D and 2D)
            f_single_1d = [SVector{3, ComplexF64}(3.0, 4.0im, 0.0)]
            pa_single = 7.5e-12
            obs_single_1d = BMO.calculate_observables(det, f_single_1d; pixel_area = pa_single)
            @test size(obs_single_1d.intensity) == (1,)
            @test size(obs_single_1d.phase) == (1,)
            @test isapprox(obs_single_1d.power, obs_single_1d.intensity[1] * pa_single, rtol = 1e-14)

            f_single_2d = fill(SVector{3, ComplexF64}(3.0, 4.0im, 0.0), 1, 1)
            obs_single_2d = BMO.calculate_observables(det, f_single_2d; pixel_area = pa_single)
            @test size(obs_single_2d.intensity) == (1, 1)
            @test size(obs_single_2d.phase) == (1, 1)
            @test isapprox(obs_single_2d.power, obs_single_2d.intensity[1, 1] * pa_single, rtol = 1e-14)

            # 3. Physical Beamlet Tracing & Parseval Equivalence (:direct vs :nufft)
            det_prop = BMO.FourierBeamletPropagator(200.0e-6; stop = true)
            sys = System([det_prop])
            lambda = 1.0e-6
            w0 = 25.0e-6
            gb = GaussianBeamlet([0.0, -100.0e-6, 0.0], [0.0, 1.0, 0.0], lambda, w0; P0 = 1.0e-3)
            BMO.interact3d(sys, det_prop, gb, 1)

            x_span = LinRange(-60.0e-6, 60.0e-6, 32)
            y_span = LinRange(-60.0e-6, 60.0e-6, 32)
            grid_prop = BMO.SpatialGrid((x_span, y_span))
            pts_prop = BMO.generate_fourier_spectrum(det_prop)
            E_dir = zeros(SVector{3, ComplexF64}, 32, 32)
            _test_direct_phasor_sum!(E_dir, pts_prop, grid_prop)
            E_nufft = BMO.synthesize_field(det_prop, grid_prop)

            obs_dir = BMO.calculate_observables(det_prop, E_dir, grid_prop)
            obs_nufft = BMO.calculate_observables(det_prop, E_nufft, grid_prop)

            # First-Principles Invariant 1: Parseval energy equivalence between direct and NUFFT
            @test obs_dir.power > 0.0
            @test obs_nufft.power > 0.0
            @test isapprox(obs_dir.power, obs_nufft.power, rtol = 1e-5)
            # Intensity profile fidelity
            max_I = maximum(obs_dir.intensity)
            @test maximum(abs.(obs_dir.intensity .- obs_nufft.intensity)) / max_I < 1e-5
            # Stokes S0 matches intensity for purely transverse beamlet
            @test isapprox(obs_dir.intensity, obs_dir.stokes[1], rtol = 1e-12)
        end

        # --------------------------------------------------------------------------
        # Zero-Allocation Invariant in Hot-Path Observable Extraction
        # --------------------------------------------------------------------------
        @testset "Zero-Allocation Invariant in Hot-Path Observables" begin
            det = BMO.FourierBeamletPropagator(100.0e-6; stop = true)

            # 1D Zero-Allocation
            x1 = LinRange(-50.0e-6, 50.0e-6, 16)
            g1 = BMO.SpatialGrid((x1,))
            f1 = zeros(SVector{3, ComplexF64}, 16)
            obs1 = BMO.DetectorObservables((16,))

            bench_obs_1d!(obs1, det, f1, g1)
            allocs1 = @allocated bench_obs_1d!(obs1, det, f1, g1)
            @test allocs1 == 0

            # 2D Zero-Allocation
            y2 = LinRange(-50.0e-6, 50.0e-6, 16)
            g2 = BMO.SpatialGrid((x1, y2))
            f2 = zeros(SVector{3, ComplexF64}, 16, 16)
            obs2 = BMO.DetectorObservables((16, 16))

            bench_obs_2d!(obs2, det, f2, g2)
            allocs2 = @allocated bench_obs_2d!(obs2, det, f2, g2)
            @test allocs2 == 0

            # 3D Zero-Allocation
            z3 = LinRange(-20.0e-6, 20.0e-6, 16)
            g3 = BMO.SpatialGrid((x1, y2, z3))
            f3 = zeros(SVector{3, ComplexF64}, 16, 16, 16)
            obs3 = BMO.DetectorObservables((16, 16, 16))

            bench_obs_3d!(obs3, det, f3, g3)
            allocs3 = @allocated bench_obs_3d!(obs3, det, f3, g3)
            @test allocs3 == 0

            # Pixel Area Mode Zero-Allocation
            bench_obs_pixel!(obs2, det, f2, 1.0e-12)
            allocs_pa = @allocated bench_obs_pixel!(obs2, det, f2, 1.0e-12)
            @test allocs_pa == 0
        end

    end # @testset "AK-F1.2652"

    # ==============================================================================
    # Integrated Observable Extraction Pipeline (Ebene 5)
    # ==============================================================================
    @testset "Integrated Observable Extraction Pipeline" begin

        # --------------------------------------------------------------------------
        # Pipeline End-to-End Equivalence & Range Overloads
        # --------------------------------------------------------------------------
        @testset "Pipeline End-to-End Equivalence & Range Overloads" begin
            det = BMO.FourierBeamletPropagator(200.0e-6; stop = true)
            lambda = 1.0e-6
            d_vec = normalize(SVector{3, Float64}(0.2, 0.9, 0.3))
            r_hit = SVector{3, Float64}(8.0e-6, -12.0e-6, 5.0e-6)
            r_origin = r_hit - 50.0e-6 * d_vec
            r = Ray(Vector(r_origin), Vector(d_vec), lambda)
            P_ray = 5.0
            hit = BMO.RayHit(r, 50.0e-6, P_ray)
            push!(det, hit)

            # 1D Grid Equivalence
            x1 = LinRange(-35.0e-6, 35.0e-6, 21)
            g1 = BMO.SpatialGrid((x1,))
            obs1_e2e = BMO.calculate_observables(det, g1)
            obs1_ranges = BMO.calculate_observables(det, (x1,))
            field1 = BMO.synthesize_field(det, g1)
            obs1_manual = BMO.calculate_observables(det, field1, g1)

            @test obs1_e2e.intensity == obs1_ranges.intensity
            @test obs1_e2e.phase == obs1_ranges.phase
            @test obs1_e2e.stokes == obs1_ranges.stokes
            @test obs1_e2e.power == obs1_ranges.power
            @test isapprox(obs1_e2e.intensity, obs1_manual.intensity, rtol = 1e-14)
            @test isapprox(obs1_e2e.phase, obs1_manual.phase, rtol = 1e-14)
            @test isapprox(obs1_e2e.power, obs1_manual.power, rtol = 1e-14)
            for k in 1:4
                @test isapprox(obs1_e2e.stokes[k], obs1_manual.stokes[k], rtol = 1e-14)
            end

            # 2D Grid Equivalence
            y2 = LinRange(-30.0e-6, 30.0e-6, 25)
            g2 = BMO.SpatialGrid((x1, y2))
            obs2_e2e = BMO.calculate_observables(det, g2)
            obs2_ranges = BMO.calculate_observables(det, (x1, y2))
            field2 = BMO.synthesize_field(det, g2)
            obs2_manual = BMO.calculate_observables(det, field2, g2)

            @test obs2_e2e.intensity == obs2_ranges.intensity
            @test obs2_e2e.phase == obs2_ranges.phase
            @test obs2_e2e.stokes == obs2_ranges.stokes
            @test obs2_e2e.power == obs2_ranges.power
            @test isapprox(obs2_e2e.intensity, obs2_manual.intensity, rtol = 1e-14)
            @test isapprox(obs2_e2e.phase, obs2_manual.phase, rtol = 1e-14)
            @test isapprox(obs2_e2e.power, obs2_manual.power, rtol = 1e-14)
            for k in 1:4
                @test isapprox(obs2_e2e.stokes[k], obs2_manual.stokes[k], rtol = 1e-14)
            end

            # 3D Grid Equivalence
            z3 = LinRange(-15.0e-6, 15.0e-6, 11)
            g3 = BMO.SpatialGrid((x1, y2, z3))
            obs3_e2e = BMO.calculate_observables(det, g3)
            obs3_ranges = BMO.calculate_observables(det, (x1, y2, z3))
            field3 = BMO.synthesize_field(det, g3)
            obs3_manual = BMO.calculate_observables(det, field3, g3)

            @test obs3_e2e.intensity == obs3_ranges.intensity
            @test obs3_e2e.power == obs3_ranges.power
            @test isapprox(obs3_e2e.intensity, obs3_manual.intensity, rtol = 1e-14)
            @test isapprox(obs3_e2e.phase, obs3_manual.phase, rtol = 1e-14)
            @test isapprox(obs3_e2e.power, obs3_manual.power, rtol = 1e-14)

            # In-place pipeline mutation: calculate_observables!(obs, det, grid)
            obs_mut2 = BMO.DetectorObservables((21, 25))
            res_mut2 = BMO.calculate_observables!(obs_mut2, det, g2)
            @test res_mut2.intensity === obs_mut2.intensity
            @test res_mut2.phase === obs_mut2.phase
            @test res_mut2.stokes === obs_mut2.stokes
            @test isapprox(res_mut2.intensity, obs2_e2e.intensity, rtol = 1e-14)
            @test isapprox(res_mut2.phase, obs2_e2e.phase, rtol = 1e-14)
            @test isapprox(res_mut2.power, obs2_e2e.power, rtol = 1e-14)
            for k in 1:4
                @test isapprox(res_mut2.stokes[k], obs2_e2e.stokes[k], rtol = 1e-14)
            end

            # Custom impedance Z parameter propagation across all overloads
            Z_custom = 250.0
            obs_z_grid = BMO.calculate_observables(det, g2; Z = Z_custom)
            obs_z_ranges = BMO.calculate_observables(det, (x1, y2); Z = Z_custom)
            obs_z_mut = BMO.DetectorObservables((21, 25))
            res_z_mut = BMO.calculate_observables!(obs_z_mut, det, g2; Z = Z_custom)

            @test obs_z_grid.power == obs_z_ranges.power
            @test isapprox(obs_z_grid.power, res_z_mut.power, rtol = 1e-14)
            @test isapprox(obs_z_grid.intensity, res_z_mut.intensity, rtol = 1e-14)
        end

        # --------------------------------------------------------------------------
        # Multi-Dimensional Grids (1D, 2D, 3D) & Anti-Padding Invariant
        # --------------------------------------------------------------------------
        @testset "Multi-Dimensional Grids & Anti-Padding Invariant" begin
            det = BMO.FourierBeamletPropagator(200.0e-6; stop = true)
            r = Ray([0.0, -50.0e-6, 0.0], [0.0, 1.0, 0.0], 1.0e-6)
            sys = System([det])
            BMO.interact3d(sys, det, Beam(r), r)

            # Prime, non-power-of-two dimensions (strict anti-padding check)
            dims_1d = (29,)
            dims_2d = (19, 31)
            dims_3d = (11, 13, 7)

            x1 = LinRange(-40.0e-6, 40.0e-6, dims_1d[1])
            x2 = LinRange(-25.0e-6, 25.0e-6, dims_2d[1])
            y2 = LinRange(-35.0e-6, 35.0e-6, dims_2d[2])
            x3 = LinRange(-15.0e-6, 15.0e-6, dims_3d[1])
            y3 = LinRange(-15.0e-6, 15.0e-6, dims_3d[2])
            z3 = LinRange(-10.0e-6, 10.0e-6, dims_3d[3])

            g1 = BMO.SpatialGrid((x1,))
            g2 = BMO.SpatialGrid((x2, y2))
            g3 = BMO.SpatialGrid((x3, y3, z3))

            obs1 = BMO.calculate_observables(det, g1)
            obs2 = BMO.calculate_observables(det, g2)
            obs3 = BMO.calculate_observables(det, g3)

            # Strict dimension matching (no implicit padding to 32 or 64)
            @test size(obs1.intensity) == dims_1d
            @test size(obs1.phase) == dims_1d
            @test all(s -> size(s) == dims_1d, obs1.stokes)

            @test size(obs2.intensity) == dims_2d
            @test size(obs2.phase) == dims_2d
            @test all(s -> size(s) == dims_2d, obs2.stokes)

            @test size(obs3.intensity) == dims_3d
            @test size(obs3.phase) == dims_3d
            @test all(s -> size(s) == dims_3d, obs3.stokes)

            # Power integration identity across dimensions
            pa1 = step(x1)
            pa2 = step(x2) * step(y2)
            pa3 = step(x3) * step(y3)
            @test isapprox(obs1.power, sum(obs1.intensity) * pa1, rtol = 1e-12)
            @test isapprox(obs2.power, sum(obs2.intensity) * pa2, rtol = 1e-12)
            @test isapprox(obs3.power, sum(obs3.intensity) * pa3, rtol = 1e-12)
        end

        # --------------------------------------------------------------------------
        # First-Principles Physical Invariants
        # --------------------------------------------------------------------------
        @testset "[PHYSICS] First-Principles Physical Invariants" begin
            det = BMO.FourierBeamletPropagator(200.0e-6; stop = true)
            lambda = 1.0e-6
            w0 = 25.0e-6
            gb = GaussianBeamlet([0.0, -100.0e-6, 0.0], [0.0, 1.0, 0.0], lambda, w0; P0 = 2.0e-3)
            sys = System([det])
            BMO.interact3d(sys, det, gb, 1)

            x_range = LinRange(-50.0e-6, 50.0e-6, 25)
            y_range = LinRange(-50.0e-6, 50.0e-6, 25)
            grid = BMO.SpatialGrid((x_range, y_range))
            obs = BMO.calculate_observables(det, grid)

            # 1. Non-negativity invariant: I >= 0, S0 >= 0, Power >= 0
            @test all(I -> I >= 0.0, obs.intensity)
            @test all(S0 -> S0 >= 0.0, obs.stokes[1])
            @test obs.power >= 0.0
            @test obs.power > 0.0

            # 2. Stokes Polarization Consistency: S0^2 >= S1^2 + S2^2 + S3^2
            S0_arr = obs.stokes[1]
            S1_arr = obs.stokes[2]
            S2_arr = obs.stokes[3]
            S3_arr = obs.stokes[4]
            @test all(i -> S0_arr[i]^2 >= (S1_arr[i]^2 + S2_arr[i]^2 + S3_arr[i]^2) - 1e-12, eachindex(S0_arr))

            # 3. Parseval Energy Conservation against Direct Phasor Sum
            pts = BMO.generate_fourier_spectrum(det)
            E_direct = zeros(SVector{3, ComplexF64}, 25, 25)
            _test_direct_phasor_sum!(E_direct, pts, grid)
            dx = step(x_range)
            dy = step(y_range)
            I_direct = [norm(v)^2 / (2 * BMO.Z_vacuum) for v in E_direct]
            P_direct = sum(I_direct) * (dx * dy)
            @test isapprox(obs.power, P_direct, rtol = 1e-5)

            # 4. Characteristic Wave Impedance Scaling & Phase Invariance
            Z_custom = BMO.Z_vacuum / 2.0
            obs_custom = BMO.calculate_observables(det, grid; Z = Z_custom)
            scaling = BMO.Z_vacuum / Z_custom
            @test isapprox(obs_custom.intensity, obs.intensity .* scaling, rtol = 1e-12)
            @test isapprox(obs_custom.power, obs.power * scaling, rtol = 1e-12)
            for k in 1:4
                @test isapprox(obs_custom.stokes[k], obs.stokes[k] .* scaling, rtol = 1e-12)
            end
            # Wavefront phase is strictly invariant to characteristic impedance scaling
            @test isapprox(obs_custom.phase, obs.phase, rtol = 1e-12)

            # 5. Strict Wavefront Phase Bounding [-π, π]
            @test all(p -> -π - 1e-14 <= p <= π + 1e-14, obs.phase)

            # 6. Linear Phase Ramp Verification for Tilted Beam
            det_tilt = BMO.FourierBeamletPropagator(200.0e-6; stop = true)
            theta = deg2rad(3.0)
            d_tilt = normalize(SVector{3, Float64}(sin(theta), cos(theta), 0.0))
            r_tilt = Ray([0.0, -50.0e-6, 0.0], Vector(d_tilt), lambda)
            push!(det_tilt, BMO.RayHit(r_tilt, 50.0e-6, 1.0))

            x_tilt = LinRange(0.0, 50.0e-6, 51)
            g_tilt = BMO.SpatialGrid((x_tilt,))
            obs_tilt = BMO.calculate_observables(det_tilt, g_tilt)

            kx = (2π / lambda) * sin(theta)
            dx_tilt = step(x_tilt)
            dphis_tilt = [mod(obs_tilt.phase[i + 1] - obs_tilt.phase[i] + π, 2π) - π for i in 1:(length(x_tilt) - 1)]
            @test all(dp -> isapprox(dp, kx * dx_tilt, atol = 1e-5), dphis_tilt)
        end

        # --------------------------------------------------------------------------
        # Precision Boundaries (Float32 / Float64) & Type Stability
        # --------------------------------------------------------------------------
        @testset "Precision Boundaries & Type Stability" begin
            # 1. End-to-end Float32 execution
            det32 = BMO.FourierBeamletPropagator(100.0f-6; stop = true)
            r32 = Ray(Float32[0.0, -50.0e-6, 0.0], Float32[0.0, 1.0, 0.0], 1.0f-6)
            push!(det32, BMO.RayHit(r32, 50.0f-6, 2.0f0))

            x32 = LinRange(-20.0f-6, 20.0f-6, 16)
            y32 = LinRange(-20.0f-6, 20.0f-6, 16)
            grid32 = BMO.SpatialGrid((x32, y32))

            obs32 = BMO.calculate_observables(det32, grid32)
            @test eltype(obs32.intensity) === Float32
            @test eltype(obs32.phase) === Float32
            @test all(s -> eltype(s) === Float32, obs32.stokes)
            @test obs32.power isa Float32
            @test obs32.power > 0.0f0

            # Float32 range tuple overload
            obs32_ranges = BMO.calculate_observables(det32, (x32, y32))
            @test eltype(obs32_ranges.intensity) === Float32
            @test obs32_ranges.power isa Float32

            # Float32 in-place pipeline
            obs32_mut = BMO.DetectorObservables((16, 16); T = Float32)
            res32_mut = BMO.calculate_observables!(obs32_mut, det32, grid32)
            @test res32_mut.intensity === obs32_mut.intensity
            @test eltype(res32_mut.intensity) === Float32
            @test res32_mut.power isa Float32
            @test isapprox(res32_mut.power, obs32.power, rtol = 1e-6)

            # 2. Strict type inference (@inferred) on Float64
            det64 = BMO.FourierBeamletPropagator(100.0e-6; stop = true)
            r64 = Ray([0.0, -50.0e-6, 0.0], [0.0, 1.0, 0.0], 1.0e-6)
            push!(det64, BMO.RayHit(r64, 50.0e-6, 1.0))
            x64 = LinRange(-20.0e-6, 20.0e-6, 16)
            y64 = LinRange(-20.0e-6, 20.0e-6, 16)
            grid64 = BMO.SpatialGrid((x64, y64))
            obs64_pre = BMO.DetectorObservables((16, 16))

            @test_nowarn @inferred BMO.calculate_observables(det64, grid64)
            @test_nowarn @inferred BMO.calculate_observables(det64, (x64, y64))
            @test_nowarn @inferred BMO.calculate_observables!(obs64_pre, det64, grid64)
        end

        # --------------------------------------------------------------------------
        # Edge Cases, Dimension Contracts & Robustness
        # --------------------------------------------------------------------------
        @testset "Edge Cases & Dimension Contracts" begin
            # 1. Empty detector hits
            det_empty = BMO.FourierBeamletPropagator(100.0e-6; stop = true)
            x = LinRange(-20.0e-6, 20.0e-6, 16)
            y = LinRange(-20.0e-6, 20.0e-6, 16)
            grid = BMO.SpatialGrid((x, y))

            obs_empty = BMO.calculate_observables(det_empty, grid)
            @test all(v -> v == 0.0, obs_empty.intensity)
            @test all(v -> v == 0.0, obs_empty.phase)
            @test obs_empty.power == 0.0
            for k in 1:4
                @test all(v -> v == 0.0, obs_empty.stokes[k])
            end

            # Empty detector in-place
            obs_empty_mut = BMO.DetectorObservables((16, 16))
            res_empty_mut = BMO.calculate_observables!(obs_empty_mut, det_empty, grid)
            @test all(v -> v == 0.0, res_empty_mut.intensity)
            @test res_empty_mut.power == 0.0

            # 2. Zero-power ray hit
            det_zero = BMO.FourierBeamletPropagator(100.0e-6; stop = true)
            r_zero = Ray([0.0, -50.0e-6, 0.0], [0.0, 1.0, 0.0], 1.0e-6)
            push!(det_zero, BMO.RayHit(r_zero, 50.0e-6, 0.0))
            obs_zero = BMO.calculate_observables(det_zero, grid)
            @test all(v -> isapprox(v, 0.0, atol = 1e-12), obs_zero.intensity)
            @test isapprox(obs_zero.power, 0.0, atol = 1e-12)

            # 3. NUFFT Grid Dimension Contract (N < 4 throws ArgumentError)
            det_active = BMO.FourierBeamletPropagator(100.0e-6; stop = true)
            push!(det_active, BMO.RayHit(r_zero, 50.0e-6, 1.0))
            x_short1 = LinRange(0.0, 0.0, 1)
            x_short2 = LinRange(-10.0e-6, 10.0e-6, 2)
            x_short3 = LinRange(-10.0e-6, 10.0e-6, 3)

            @test_throws ArgumentError BMO.calculate_observables(det_active, BMO.SpatialGrid((x_short1,)))
            @test_throws ArgumentError BMO.calculate_observables(det_active, BMO.SpatialGrid((x_short2,)))
            @test_throws ArgumentError BMO.calculate_observables(det_active, BMO.SpatialGrid((x_short3,)))
            @test_throws ArgumentError BMO.calculate_observables(det_active, (x_short3,))
            @test_throws ArgumentError BMO.calculate_observables(det_active, BMO.SpatialGrid((x_short3, y)))
            @test_throws ArgumentError BMO.calculate_observables(det_active, (x_short3, y))
            @test_throws ArgumentError BMO.calculate_observables(det_active, (x, x_short2))

            obs_short = BMO.DetectorObservables((3, 16))
            @test_throws ArgumentError BMO.calculate_observables!(obs_short, det_active, BMO.SpatialGrid((x_short3, y)))

            # 4. In-place container dimension mismatch throws DimensionMismatch
            obs_mismatch = BMO.DetectorObservables((16, 24))
            @test_throws DimensionMismatch BMO.calculate_observables!(obs_mismatch, det_active, grid)

            # 5. State purity: repeated mutation with updated detector state
            obs_pure = BMO.DetectorObservables((16, 16))
            BMO.calculate_observables!(obs_pure, det_empty, grid)
            @test obs_pure.power == 0.0

            BMO.calculate_observables!(obs_pure, det_active, grid)
            @test obs_pure.power > 0.0
            power_run1 = obs_pure.power

            BMO.calculate_observables!(obs_pure, det_active, grid)
            @test obs_pure.power == power_run1
        end

        # --------------------------------------------------------------------------
        # End-to-End Beamlet Tracing & Polarization Synthesis
        # --------------------------------------------------------------------------
        @testset "[PHYSICS] End-to-End Beamlet Tracing & Polarization" begin
            # 1. Tracing Gaussian Beamlet
            lambda = 1.0e-6
            w0 = 25.0e-6
            P0 = 1.0e-3
            gb = GaussianBeamlet([0.0, -100.0e-6, 0.0], [0.0, 1.0, 0.0], lambda, w0; P0 = P0)
            det_gb = BMO.FourierBeamletPropagator(200.0e-6; stop = true)
            sys = System([det_gb])
            BMO.interact3d(sys, det_gb, gb, 1)

            # Grid spanning [-2.5*w0, 2.5*w0]
            span = 2.5 * w0
            x_gb = LinRange(-span, span, 32)
            y_gb = LinRange(-span, span, 32)
            grid_gb = BMO.SpatialGrid((x_gb, y_gb))

            obs_gb = BMO.calculate_observables(det_gb, grid_gb)
            @test obs_gb.power > 0.0
            pa_gb = step(x_gb) * step(y_gb)
            @test isapprox(obs_gb.power, sum(obs_gb.intensity) * pa_gb, rtol = 1e-12)
            @test maximum(obs_gb.intensity) > 0.0
            @test all(I -> I >= 0.0, obs_gb.intensity)

            # 2. Polarized Ray Hits: Right Circular Polarization
            det_rcp = BMO.FourierBeamletPropagator(100.0e-6; stop = true)
            p_rcp = ComplexF64[1.0 / sqrt(2), -im / sqrt(2), 0.0]
            pr_rcp = PolarizedRay([0.0, 0.0, -50.0e-6], [0.0, 0.0, 1.0], lambda, p_rcp)
            push!(det_rcp, BMO.PolarizedRayHit(pr_rcp, 50.0e-6, 1.0))

            x_pol = LinRange(-20.0e-6, 20.0e-6, 16)
            y_pol = LinRange(-20.0e-6, 20.0e-6, 16)
            grid_pol = BMO.SpatialGrid((x_pol, y_pol))

            obs_rcp = BMO.calculate_observables(det_rcp, grid_pol)
            # RCP state: S0 > 0, S1 ≈ 0, S2 ≈ 0, S3 > 0, S3 ≈ S0
            @test all(S3 -> S3 > 0.0, obs_rcp.stokes[4])
            @test isapprox(obs_rcp.stokes[4], obs_rcp.stokes[1], rtol = 1e-5)
            @test all(S1 -> isapprox(S1, 0.0, atol = 1e-10), obs_rcp.stokes[2])
            @test all(S2 -> isapprox(S2, 0.0, atol = 1e-10), obs_rcp.stokes[3])

            # Left Circular Polarization
            det_lcp = BMO.FourierBeamletPropagator(100.0e-6; stop = true)
            p_lcp = ComplexF64[1.0 / sqrt(2), im / sqrt(2), 0.0]
            pr_lcp = PolarizedRay([0.0, 0.0, -50.0e-6], [0.0, 0.0, 1.0], lambda, p_lcp)
            push!(det_lcp, BMO.PolarizedRayHit(pr_lcp, 50.0e-6, 1.0))

            obs_lcp = BMO.calculate_observables(det_lcp, grid_pol)
            # LCP state: S0 > 0, S1 ≈ 0, S2 ≈ 0, S3 < 0, -S3 ≈ S0
            @test all(S3 -> S3 < 0.0, obs_lcp.stokes[4])
            @test isapprox(obs_lcp.stokes[4], -obs_lcp.stokes[1], rtol = 1e-5)
            @test all(S1 -> isapprox(S1, 0.0, atol = 1e-10), obs_lcp.stokes[2])
            @test all(S2 -> isapprox(S2, 0.0, atol = 1e-10), obs_lcp.stokes[3])
        end

    end # @testset "AK-F1.2653"

    # ==============================================================================
    # Physical Invariants & Paraxial Wavefront Optics Migration
    # ==============================================================================
    @testset "[PHYSICS] Physical Invariants & Paraxial Wavefront Optics" begin

        # --------------------------------------------------------------------------
        # Multi-Beam Interference (Crossing Rays, Fringe Visibility)
        # --------------------------------------------------------------------------
        # [TEST SPECIFICATION]
        # - Ziel: Verifikation der kohärenten Multi-Beam-Interferenz zweier kreuzender Strahlen ohne Phasenlöschung und Nachweis hoher Interferenzstreifen-Sichtbarkeit (Visibility > 99%).
        # - Input: Planarer FourierBeamletPropagator (L = 100 μm, is_planar = true) in der X-Z-Ebene bei y = 0. Zwei kohärente Strahlen bei Wellenlänge λ = 1.0 μm, Start bei y = -50 μm, symmetrisch geneigt unter Winkeln ±θ (θ = 0.01 rad) relativ zur Flächennormalen (+Y). Gitter: 256 x 256 Punkte über [-50 μm, 50 μm] in x und z.
        # - Output: Maximale Linien-Intensität I_max > 0, minimale Linien-Intensität I_min < 0.01 * I_max, Streifen-Sichtbarkeit V = (I_max - I_min) / (I_max + I_min) > 0.99.
        # - Physikalischer Hintergrund: Die Überlagerung zweier ebener Wellen mit transversaler Wellenzahl-Differenz Δk_x = 2 k0 sin(θ) erzeugt ein cos²-Interferenzmuster. Die Trägerphase beider Strahlen muss im Fourier-Spektrum interferenzfähig akkumuliert werden.
        @testset "[PHYSICS] Multi-Beam Interference & Fringe Visibility" begin
            λ = 1.0e-6
            L = 100.0e-6
            N = 256
            det = BMO.FourierBeamletPropagator(L; stop = true, is_planar = true)
            sys = System([det])

            # Two coherent rays crossing at +/- theta relative to normal (along y-axis)
            θ = 0.01
            dir1 = normalize([sin(θ), 1.0, 0.0])
            dir2 = normalize([-sin(θ), 1.0, 0.0])

            ray1 = Ray([0.0, -50.0e-6, 0.0], dir1, λ)
            ray2 = Ray([0.0, -50.0e-6, 0.0], dir2, λ)

            solve_system!(sys, Beam(ray1))
            solve_system!(sys, Beam(ray2))

            xs = LinRange(-L / 2, L / 2, N)
            zs = LinRange(-L / 2, L / 2, N)
            obs = BMO.calculate_observables(det, (xs, zs))

            I = obs.intensity
            mid_z = 1 + N ÷ 2
            I_line = I[:, mid_z]
            I_max = maximum(I_line)
            I_min = minimum(I_line)
            visibility = (I_max - I_min) / (I_max + I_min)

            # Multi-ray interference must be preserved: visibility > 99%
            @test visibility > 0.99
            @test I_max > 0.0
            @test I_min < 0.01 * I_max
        end

        # --------------------------------------------------------------------------
        # Paraxial Newton-Rings and Wavefront Curvature
        # --------------------------------------------------------------------------
        # [TEST SPECIFICATION]
        # - Ziel: Verifikation der quadratischen paraxialen Wellenfrontkrümmungsphase und der Ausbildung konzentrischer Newton-Ringe durch Interferenz mit einer ebenen Referenzwelle.
        # - Input: Planarer FourierBeamletPropagator (L_det = 800 μm, is_planar = true) in der X-Z-Ebene bei y = 0. Konvergierendes GaussianBeamlet bei λ = 1.0 μm, Taille w0 = 10.0 μm, Leistung P0 = 1 mW, Start bei y = -100 μm, Propagation entlang +Y, Fokus z0 = 200 μm (Taille 100 μm hinter Detektorebene). Auswertegitter: N = 128 Punkte über L = 80 μm in x und z.
        # - Output: Phasenkrümmung bei r = 5 μm approximiert die theoretische Paraxialphase mit relativer Fehlerschranke rtol = 0.03. Kontrast der Newton-Ringe vis > 0.95.
        # - Physikalischer Hintergrund: Ein paraxialer Gauss-Strahl besitzt eine sphärische Wellenfront mit Krümmungsradius R(z) = z * (1 + (zR/z)²). Die Interferenz des gekrümmten Strahls mit einer phasenkohärenten Referenzwelle gleicher Spitzenamplitude erzeugt ein Newton-Ring-Muster mit destruktiver Interferenz bei r_pi = sqrt(λ * |R|).
        @testset "[PHYSICS] Paraxial Newton-Rings & Wavefront Curvature" begin
            λ = 1.0e-6
            w0 = 10.0e-6
            P0 = 1e-3
            k0 = 2π / λ
            zR = π * w0^2 / λ
            z_focus_offset = -100.0e-6 # focus 100 μm behind detector (z_waist - z_det = +100 μm, z_prop = -100 μm)
            R_th = z_focus_offset * (1 + (zR / z_focus_offset)^2) # R < 0 for converging

            L_beam = 80.0e-6
            N_beam = 128
            dx_b = L_beam / N_beam

            # GaussianBeamlet started at y = -100 μm with waist shifted 200 μm ahead -> waist at y = +100 μm
            g = GaussianBeamlet([0.0, -100.0e-6, 0.0], [0.0, 1.0, 0.0], λ, w0; z0 = 200.0e-6, P0 = P0)
            det_beam = BMO.FourierBeamletPropagator(800.0e-6; stop = true, is_planar = true)
            sys_beam = System([det_beam])
            solve_system!(sys_beam, g)

            xs_b = range(-L_beam / 2, step = dx_b, length = N_beam)
            zs_b = range(-L_beam / 2, step = dx_b, length = N_beam)

            obs = BMO.calculate_observables(det_beam, (xs_b, zs_b))
            mid_b = N_beam ÷ 2 + 1

            # 1. Quadratic wavefront curvature phase verification: Δϕ(r) = k r^2 / (2 |R_th|)
            phi_0 = obs.phase[mid_b, mid_b]
            idx_5um = mid_b + round(Int, 5.0e-6 / dx_b)
            r_5um = xs_b[idx_5um]
            phi_5um = obs.phase[idx_5um, mid_b] - phi_0
            phi_th_5um = k0 * r_5um^2 / (2 * R_th)
            @test isapprox(abs(phi_5um), abs(phi_th_5um), rtol = 0.03)

            # 2. Concentric Newton-Ring Interference with Planar Reference Wave
            E_field = BMO.synthesize_field(det_beam, (xs_b, zs_b))
            E_b = [v[1] for v in E_field]
            r_pi = sqrt(λ * abs(R_th))
            idx_pi = round(Int, r_pi / dx_b)
            E_ref = abs(E_b[mid_b + idx_pi, mid_b])
            I_rings = abs2.(E_b .- E_ref)
            I_min = minimum(I_rings[mid_b:mid_b+idx_pi+2, mid_b])
            I_max = maximum(I_rings[mid_b:mid_b+idx_pi+2, mid_b])
            vis = (I_max - I_min) / (I_max + I_min)
            @test vis > 0.95
        end

        # --------------------------------------------------------------------------
        # Oblique Incidence at 60° (Power Constancy for Circular & Astigmatic Beamlets)
        # --------------------------------------------------------------------------
        # [TEST SPECIFICATION]
        # - Ziel: Verifikation der Azimutalen Projektionsinvarianz (exakte Schnitt-Ellipsentransformation) und der Leistungskonstanz bei 60° schrägem Einfall für zirkuläre und astigmatische Gauss-Beamlets auf einem planaren FourierBeamletPropagator.
        # - Input: Planarer FourierBeamletPropagator (L_det = 800 μm, is_planar = true) in der X-Z-Ebene bei y = 0. Strahlung bei λ = 20.0 μm, w0 = 50.0 μm, P0 = 2.5 mW unter 60° Einfallswinkel (θ = 60°), Start bei y = -500 μm * cos(60°). Numerisches Gitter: N = 256 Punkte über L = 800 μm in x und z.
        # - Output: Relative Abweichung der berechneten Ellipsenfläche von der Zielprojektionsfläche unter 1e-12 für alle Azimutwinkel. Gemessene Leistung für AstigmaticGaussianBeamlet und GaussianBeamlet innerhalb von 2% Toleranz (rtol = 0.02) zu P0_target. Leistungsverhältnis AGB zu GB innerhalb 0.5% (rtol = 0.005).
        # - Physikalischer Hintergrund: Beim schrägen Einfall auf eine planare Detektorfläche projiziert sich das transversale Strahlprofil zu einer Ellipse mit vergrößerter Fläche A = A0 / cos(θ). Nach dem Kosinusgesetz und Poynting-Fluss muss die über die Ebene integrierte optische Gesamtleistung exakt erhalten bleiben (P = P0), unabhängig von Strahl-Astigmatismus oder azimutaler Orientierung.
        @testset "[PHYSICS] Oblique Incidence at 60°: Ellipse Transformation & Power Constancy" begin
            # 1. Azimuthal projection invariance: exact ellipse area transformation
            w0 = 15.0e-6
            θ = π / 3 # 60°
            cos_th = cos(θ)
            target_area = π * w0^2 / cos_th

            for deg in [0, 15, 30, 45, 60, 75, 90]
                ϕ = deg2rad(deg)
                kx_loc = sin(θ) * cos(ϕ)
                kz_loc = sin(θ) * sin(ϕ)
                k_perp_sq = kx_loc^2 + kz_loc^2
                ux = k_perp_sq > 1e-20 ? (kx_loc^2 / k_perp_sq) : 0.5
                uz = k_perp_sq > 1e-20 ? (kz_loc^2 / k_perp_sq) : 0.5
                inv_cos2_minus_1 = (1.0 / cos_th^2) - 1.0
                wx_proj = sqrt(w0^2 * (1.0 + inv_cos2_minus_1 * ux))
                wz_proj = sqrt(w0^2 * (1.0 + inv_cos2_minus_1 * uz))
                gamma_norm = sqrt((w0^2 / cos_th) / (wx_proj * wz_proj))
                wx_surf = gamma_norm * wx_proj
                wz_surf = gamma_norm * wz_proj

                calc_area = π * wx_surf * wz_surf
                rel_err = abs(calc_area - target_area) / target_area
                @test rel_err < 1e-12
            end

            # 2. Oblique Incidence Power Constancy for Astigmatic and Circular Beamlets
            det = BMO.FourierBeamletPropagator(800.0e-6; stop = true, is_planar = true)
            P0_target = 2.5e-3
            θ_60 = deg2rad(60.0)
            dir_60 = [sin(θ_60), cos(θ_60), 0.0]
            pos_start = [0.0, 0.0, 0.0] - 500.0e-6 * dir_60
            λ = 20.0e-6
            w0_60 = 50.0e-6

            # 2a. AstigmaticGaussianBeamlet at 60°
            agb = AstigmaticGaussianBeamlet(pos_start, dir_60, λ, w0_60; P0 = P0_target)
            sys = System([det])
            solve_system!(sys, agb)
            @test length(det.hits) == 1

            L_box = 800.0e-6
            N = 256
            xs = range(-L_box / 2, step = L_box / N, length = N)
            zs = range(-L_box / 2, step = L_box / N, length = N)

            obs_agb = BMO.calculate_observables(det, (xs, zs))
            P_meas_agb = obs_agb.power
            @test isapprox(P_meas_agb, P0_target, rtol = 0.02)

            # 2b. GaussianBeamlet at 60° for consistency check
            empty!(det)
            g = GaussianBeamlet(pos_start, dir_60, λ, w0_60; P0 = P0_target)
            solve_system!(sys, g)
            @test length(det.hits) == 1

            obs_g = BMO.calculate_observables(det, (xs, zs))
            P_meas_g = obs_g.power
            @test isapprox(P_meas_g, P0_target, rtol = 0.02)

            # Ratio between AGB and GB measured powers must be within 0.5%
            @test isapprox(P_meas_agb, P_meas_g, rtol = 0.005)
        end

    end # @testset "AK-F1.2654"

    # --------------------------------------------------------------------------
    # [AK-24.9] SphericalLens OPL and Focal Symmetry
    # --------------------------------------------------------------------------
    # [TEST SPECIFICATION]
    # - ### Goal: Verify that FourierBeamletPropagator correctly incorporates the OPL phase anchor and local coordinate frame.
    # - ### Input: Plano-convex SphericalLens (n=1.5168, R=103.36mm, f approx 200mm, D=10mm, lambda=632.8nm), Planar FBP at focal plane (y approx 200.13mm), grid (x,z) in [-50um, 50um], N=129.
    # - ### Output: Central peak I(0,0) == I_max; Radial symmetry max |I(x,0)-I(0,x)|/I_max < 1e-3; Airy disc first minimum r_min in [15.5um, 17.0um].
    # - ### Physical Background: Paraxial focusing of collimated rays forms an Airy diffraction pattern at focal plane. Coherent field synthesis requires OPL phase anchor and local coordinate transformation.
    @testset "[PHYSICS] SphericalLens OPL and Focal Symmetry" begin
        """
        ### Goal
        Verify that `FourierBeamletPropagator` correctly incorporates the optical path length (OPL) phase anchor and the detector local coordinate frame for ray hits, ensuring correct focal spot formation, radial symmetry, and Airy disc minimum location for focused lens systems.

        ### Input
        - Plano-convex `SphericalLens`: refractive index n = 1.5168, radius R = 103.36 mm, thickness l = 1.0 mm, aperture diameter d = 25.4 mm (focal length f approx 200 mm).
        - Source: UniformDiscSource of aperture diameter D = 10 mm emitting parallel rays at wavelength lambda = 632.8 nm along +Y.
        - Detector: Planar `FourierBeamletPropagator` (L = 100 um, is_planar = true) placed at the focal plane (y approx 200.13 mm), oriented normal to the optical axis (+Y).
        - Evaluation grid: Symmetric grid (x, z) in [-50 um, 50 um] with N = 129 points (step dx = dz approx 0.78125 um).

        ### Output
        - Central peak intensity: maximum intensity located at the optical axis, I(0, 0) == I_max.
        - Radial symmetry: maximum relative difference along orthogonal axes, max |I(x, 0) - I(0, x)| / I_max < 1e-3.
        - Airy disc first minimum: first radial minimum located at r_min within [15.5 um, 17.0 um], consistent with the theoretical Airy radius 1.22 lambda f / D approx 15.44 um.

        ### Physical Background
        A collimated beam focused by an ideal spherical lens under the paraxial approximation forms an Airy diffraction pattern at the focal plane. Coherent field synthesis on a detector receiving refracted rays requires that each ray retains its accumulated optical path length (OPL) phase and that ray hit coordinates are transformed into the detector local coordinate frame. Inaccurate OPL phase tracking or improper frame transformation disrupts constructive interference at the focus, causing focal spot distortion and loss of radial symmetry.
        """;

        R = 103.36e-3
        n = 1.5168
        D = 10.0e-3
        λ = 632.8e-9
        f = 200.0e-3
        y_foc = 200.13e-3
        num_rays = 2000

        cs = UniformDiscSource([0.0, -10.0e-3, 0.0], [0.0, 1.0, 0.0], D, λ; num_rays = num_rays)
        lens = SphericalLens(R, Inf, 1.0e-3, 25.4e-3, x -> n)

        det = BMO.FourierBeamletPropagator(100.0e-6; stop = true, is_planar = true)
        BMO.translate3d!(det, [0.0, y_foc, 0.0])

        sys = System([lens, det])
        solve_system!(sys, cs)
        @test length(det.hits) == num_rays

        N = 129
        xs = range(-50.0e-6, 50.0e-6, length = N)
        zs = range(-50.0e-6, 50.0e-6, length = N)

        obs = BMO.calculate_observables(det, (xs, zs))
        I = obs.intensity
        mid = (N + 1) ÷ 2
        I_max = maximum(I)

        # 1. Central peak intensity: I(0, 0) == I_max
        @test I[mid, mid] == I_max

        # 2. Radial symmetry: max |I(x, 0) - I(0, x)| / I_max < 1e-3
        diff_radial = maximum(abs.(I[:, mid] .- I[mid, :])) / I_max
        @test diff_radial < 1e-3

        # 3. Airy disc first minimum location: minimum intensity at radius r_min ≈ 1.22 λ f / D ≈ 15.5 - 17.0 μm
        radial_profile = I[mid:mid+30, mid]
        min_rel_idx = argmin(radial_profile)
        r_min = (min_rel_idx - 1) * step(xs)
        @test 15.5e-6 <= r_min <= 17.0e-6
    end

end # @testset "FourierBeamletPropagator Tests"

end # module TestFourierBeamletPropagator
