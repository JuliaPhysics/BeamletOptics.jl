using Test
using BeamletOptics
using LinearAlgebra
const BMO = BeamletOptics

@testset "Bulk Absorption and Small-Signal Gain" begin
    λ = 1000e-9 # 1 μm

    @testset "AbsorbingMedium & GainMedium Functors" begin
        # AbsorbingMedium with constant α
        α0 = 100.0 # 100 1/m (1/cm)
        med_abs = AbsorbingMedium(1.5, α0)
        n_c = med_abs(λ)
        @test real_refractive_index(n_c) ≈ 1.5
        expected_κ = (α0 * λ) / (4π)
        @test extinction_coefficient(n_c) ≈ expected_κ
        @test absorption_coefficient(n_c, λ) ≈ α0

        # GainMedium with constant g
        g0 = 50.0 # 50 1/m
        med_gain = GainMedium(1.8, g0)
        n_gain = med_gain(λ)
        @test real_refractive_index(n_gain) ≈ 1.8
        expected_gain_κ = -(g0 * λ) / (4π)
        @test extinction_coefficient(n_gain) ≈ expected_gain_κ
        @test absorption_coefficient(n_gain, λ) ≈ -g0

        # bulk_attenuation_factor helper
        L = 0.01 # 10 mm
        att_p, att_f = bulk_attenuation_factor(n_c, λ, L)
        @test att_p ≈ exp(-α0 * L)
        @test att_f ≈ exp(-α0 * L / 2)

        att_gain_p, att_gain_f = bulk_attenuation_factor(n_gain, λ, L)
        @test att_gain_p ≈ exp(g0 * L)
        @test att_gain_f ≈ exp(g0 * L / 2)

        # Non-absorbing medium
        att_real_p, att_real_f = bulk_attenuation_factor(1.5, λ, L)
        @test att_real_p ≈ 1.0
        @test att_real_f ≈ 1.0
    end

    @testset "Ray Lambert-Beer Attenuation through Planar Slab" begin
        d = 10e-3 # 10 mm thickness
        α = 200.0 # 200 1/m
        med = AbsorbingMedium(1.5, α)

        # Flat plano-planar lens with AR coatings to isolate bulk attenuation
        slab = SphericalLens(Inf, Inf, d, 20e-3, med) |> with_coatings(
            front = SimpleARCoating(0.0),
            back = SimpleARCoating(0.0)
        )
        sys = System([slab])

        ray = Ray([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0], λ)
        beam = Beam(ray)
        solve_system!(sys, beam; retrace = false)

        # Ray 1: before slab, Ray 2: inside slab, Ray 3: after slab
        ray_inside = rays(beam)[2]
        ray_exit = rays(beam)[3]

        @test length(ray_inside) ≈ d
        @test ray_inside.weight ≈ 1.0 # Enters without loss due to AR coating
        # Exiting ray should be attenuated by Lambert-Beer factor exp(-α * d)
        @test ray_exit.weight ≈ exp(-α * d)
    end

    @testset "PolarizedRay Bulk Attenuation with Fresnel Interfaces" begin
        d = 5e-3 # 5 mm thickness
        α = 150.0 # 150 1/m
        med = AbsorbingMedium(1.5, α)

        # Uncoated plano-planar slab (Fresnel reflections at both interfaces)
        slab_uncoated = SphericalLens(Inf, Inf, d, 20e-3, med)
        sys = System([slab_uncoated])

        E0_init = [1.0, 0.0, 0.0]
        ray_pol = PolarizedRay([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0], λ, E0_init)
        beam_pol = Beam(ray_pol)
        solve_system!(sys, beam_pol; retrace = false)

        ray_inside = rays(beam_pol)[2]
        ray_exit = rays(beam_pol)[3]

        # Exact complex Fresnel transmission:
        n_c = med(λ)
        _, _, ts1, _ = BMO.fresnel_coefficients(0.0, n_c / 1.0)
        _, _, ts2, _ = BMO.fresnel_coefficients(0.0, 1.0 / n_c)
        expected_E = abs(ts1 * exp(-α * d / 2) * ts2)
        E_exit = abs(BMO.polarization(ray_exit)[1])
        @test E_exit ≈ expected_E
    end

    @testset "Laser Gain Medium (Small-Signal Gain)" begin
        d = 20e-3 # 20 mm crystal length
        g = 80.0  # 80 1/m gain
        crystal = SphericalLens(Inf, Inf, d, 10e-3, GainMedium(1.82, g)) |> with_coatings(
            front = SimpleARCoating(0.0),
            back = SimpleARCoating(0.0)
        )
        sys = System([crystal])

        ray = Ray([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0], λ)
        beam = Beam(ray)
        solve_system!(sys, beam; retrace = false)

        ray_exit = rays(beam)[3]
        @test ray_exit.weight ≈ exp(g * d)
    end

    @testset "Direct Complex Refractive Index Function" begin
        # Passing an anonymous function returning Complex directly: λ -> 1.6 + 0.0005im
        n_func = λ -> complex(1.6, 5e-4)
        α_calc = 4π * 5e-4 / λ
        d = 10e-3
        slab = SphericalLens(Inf, Inf, d, 20e-3, n_func) |> with_coatings(
            front = SimpleARCoating(0.0),
            back = SimpleARCoating(0.0)
        )
        sys = System([slab])

        ray = Ray([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0], λ)
        beam = Beam(ray)
        solve_system!(sys, beam; retrace = false)

        ray_exit = rays(beam)[3]
        @test ray_exit.weight ≈ exp(-α_calc * d)
    end

    @testset "Detector Hit Recording inside Absorbing Medium" begin
        d = 10e-3
        α = 100.0
        # Ray travelling in absorbing medium hitting a detector
        det = Detector(10e-3)
        translate3d!(det, [0.0, d, 0.0])

        sys = System([det])
        # Start ray with complex n
        n_c = complex(1.5, (α * λ) / (4π))
        ray = Ray([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], λ)
        ray.n = n_c
        beam = Beam(ray)
        solve_system!(sys, beam; retrace = false)

        hit = BMO.hits(det)[1]
        @test hit.ray.weight ≈ exp(-α * d)
    end
end
