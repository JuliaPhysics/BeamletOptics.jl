module TestGaussianBeamlet

using BeamletOptics
using Test

const BMO = BeamletOptics

const mm = 1e-3

@testset "Gaussian beamlet" begin
    @testset "Testing type definitions" begin
        @test isdefined(BMO, :GaussianBeamlet)
    end

    @testset "Testing analytical equations" begin
        λ = 500e-9
        w0 = 1mm
        M2 = 1
        zR = BMO.rayleigh_range(λ, w0, M2)
        # Test Rayleigh range and div. angle against Paschotta (https://www.rp-photonics.com/gaussian_beams.html)
        @test isapprox(zR, 6.28, atol = 1e-2)
        @test isapprox(BMO.beam_waist(zR, w0, zR), sqrt(2) * w0)
        @test isapprox(BMO.gouy_phase(zR, zR), -π / 4)
        @test isapprox(BMO.wavefront_curvature(zR, zR), 1 / (2 * zR))
        @test isapprox(BMO.divergence_angle(λ, w0, M2), 159e-6, atol = 1e-6)
    end

    @testset "Testing parameter correctness" begin
        # Gauss beam parameters
        y = -5:0.01:5       # m
        λ_1 = 500e-9        # m
        λ_2 = 1000e-9       # m
        P0 = 1              # W
        r = 0               # m
        w0_1 = 1mm         # m
        w0_2 = 2mm         # m
        M2_1 = 1mm         # m
        M2_2 = 2mm         # m
        E0_1 = BMO.electric_field(2 * P0 / (π * w0_1^2))
        E0_2 = BMO.electric_field(2 * P0 / (π * w0_2^2))
        gauss_1 = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0],
            λ_1,
            w0_1,
            M2 = M2_1,
            P0 = P0)
        gauss_2 = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0],
            λ_2,
            w0_2,
            M2 = M2_2,
            P0 = P0)
        # Calculate analytical values
        zr_1 = BMO.rayleigh_range(λ_1, w0_1, M2_1)
        zr_2 = BMO.rayleigh_range(λ_2, w0_2, M2_2)
        wa_1 = BMO.beam_waist.(y, w0_1, zr_1)
        wa_2 = BMO.beam_waist.(y, w0_2, zr_2)
        Ra_1 = BMO.wavefront_curvature.(y, zr_1)
        Ra_2 = BMO.wavefront_curvature.(y, zr_2)
        ψa_1 = BMO.gouy_phase.(y, zr_1)
        ψa_2 = BMO.gouy_phase.(y, zr_2)
        Ea_1 = BMO.electric_field.(r, y, E0_1, w0_1, λ_1, M2_1)
        Ea_2 = BMO.electric_field.(r, y, E0_2, w0_2, λ_2, M2_2)
        # Calculate numerical values
        wn_1, Rn_1, ψn_1, w0n_1 = BMO.gauss_parameters(gauss_1, y)
        wn_2, Rn_2, ψn_2, w0n_2 = BMO.gauss_parameters(gauss_2, y)
        En_1 = [BMO.electric_field(gauss_1, r, yi) for yi in y]
        En_2 = [BMO.electric_field(gauss_2, r, yi) for yi in y]
        # Compare beam diameter within 0.1 nm
        @test all(isapprox.(wa_1, wn_1, atol = 1e-10))
        @test all(isapprox.(wa_2, wn_2, atol = 1e-10))
        # Compare wavefront curvature
        @test all(isapprox.(Ra_1, Rn_1, atol = 5e-9))
        @test all(isapprox.(Ra_2, Rn_2, atol = 5e-9))
        # Compare Gouy phase within 0.1 μrad
        @test all(isapprox.(ψa_1, ψn_1, atol = 1e-7))
        @test all(isapprox.(ψa_2, ψn_2, atol = 1e-7))
        # Compare calculated waist radius
        @test all(isapprox.(w0_1, w0n_1))
        @test all(isapprox.(w0_2, w0n_2))
        # Compare calculate electric field at r
        @test all(isapprox.(Ea_1, En_1, atol = 1e-8))
        @test all(isapprox.(Ea_2, En_2, atol = 1e-7))
        # Compare beam power with original value
        @test P0 ≈ BMO.optical_power(gauss_1)
        @test P0 ≈ BMO.optical_power(gauss_2)
    end

    @testset "Testing propagation correctness" begin
        # Analytical result using complex q factor
        q_ana(q0::Complex, M::Matrix) = (M[1] * q0 + M[3]) / (M[2] * q0 + M[4])
        R_ana(q::Complex) = real(1 / q)
        w_ana(q::Complex, λ, n = 1) = sqrt(-λ / (π * n * imag(1 / q)))
        propagate_ABCD(d) = [1 d; 0 1]
        lensmaker_ABCD(f) = [1 0; -1/f 1]
        # Beam parameters
        λ = 1000e-9
        w0 = 1mm
        M2 = 1
        zr = BMO.rayleigh_range(λ, w0, M2)
        # Lens parameters
        R1 = 1
        R2 = 1
        lens_y_location = 0.1
        nl = 1.5
        f = BMO.lensmakers_eq(R1, -R2, nl)
        # Stuff
        dy = 0.001
        ys = 0:dy:1.5
        w_analytical = Vector{Float64}(undef, length(ys))
        R_analytical = Vector{Float64}(undef, length(ys))
        # Propagate using ABCD formalism
        q0 = 0 + zr * im
        for i in 1:length(ys)
            w_analytical[i] = w_ana(q0, λ)
            R_analytical[i] = R_ana(q0)
            # catch first lens
            if i * dy == lens_y_location
                q0 = q_ana(q0, lensmaker_ABCD(f))
                continue
            end
            q0 = q_ana(q0, propagate_ABCD(dy))
        end

        # Numerical result
        tl = BMO.ThinLensSDF(R1, R2, 0.025)
        lens = Lens(tl, x -> nl)
        system = System(lens)
        translate3d!(lens, [0, lens_y_location, 0])
        # Create and solve beam, calculate beam parameters
        gauss = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0],
            λ,
            w0,
            support = [1, 0, 0],
            M2 = 1)
        solve_system!(system, gauss)
        w_numerical, R_numerical, ψ_numerical, w0_numerical = BMO.gauss_parameters(
            gauss,
            ys)
        # Compare beam radius to within 1 μm
        @test all(isapprox.(w_analytical, w_numerical, atol = 1e-6))
        # Compare if radius of curvature agreement to within 1 cm above 95% and any NaN
        temp = isapprox.(R_analytical, R_numerical, atol = 1e-2)
        @test sum(temp) / length(temp) > 0.95
        @test !any(isnan.(R_numerical))
        # Compare if Gouy phase zero at waists
        w0, i = findmin(w_analytical)
        @test isapprox(ψ_numerical[1], 0, atol = 1mm)
        @test isapprox(ψ_numerical[i], 0, atol = 1mm)
        # Compare calculated waist after lens
        @test isapprox(w0_numerical[i], w0, atol = 1e-7)

        @testset "Testing isparaxial and istilted" begin
            # Before lens rotation
            @test BMO.istilted(system, gauss) == false
            @test BMO.isparaxial(system, gauss) == true
            # Tilt lens, test again with 30° threshold for paraxial approx.
            zrotate3d!(lens, deg2rad(45))
            solve_system!(system, gauss)
            @test BMO.istilted(system, gauss) == true
            @test BMO.isparaxial(system, gauss, deg2rad(30)) == false
        end
    end
end

end # MODULE