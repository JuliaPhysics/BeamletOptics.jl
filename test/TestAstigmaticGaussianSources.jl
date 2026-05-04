module TestAstigmaticGaussianSources

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics
const mm = 1e-3
const nm = 1e-9

@testset "Astigmatic Gaussian Beamlet Sources" begin
    λ = 532nm
    dir = [0.0, 1.0, 0.0]
    pos = [0.0, 0.0, 0.0]

    @testset "CollimatedGaussianBeamletSource" begin
        D = 2mm
        w0s = D / 10
        bg = CollimatedGaussianBeamletSource(pos, dir, D, λ, w0s; n_grid = 10)
        @test bg isa AstigmaticBeamGroup
        beams = BMO.beams(bg)

        @test length(beams) == 100
        # Check all point in the right direction
        @test all(normalize(BMO.direction(b)) ≈ dir for b in bg)
        # Check all have the correct sub-waist
        for b in beams
            _, w1, _ = BMO.waist_parameters(b, 0.0)
            @test norm(w1) ≈ w0s
        end
    end

    @testset "GaussianBeamletDecomposition" begin
        w0 = 1mm
        # A 21x21 grid will trim beamlets with amp < 1e-4, so length will be < 441
        bg = GaussianBeamletDecomposition(pos, dir, w0, λ; n_grid = 21)
        @test bg isa AstigmaticBeamGroup
        beams = BMO.beams(bg)

        @test length(beams) > 0
        @test length(beams) <= 441

        # Find the beamlet with the maximum amplitude
        center_beam = argmax([norm(BMO.polarization(b)) for b in bg])

        # The peak amplitude beamlet should be at the center (0, 0, 0)
        @test norm(BMO.position(bg[center_beam]) - pos)≈0.0 atol=1e-5
    end

    @testset "SphericalGaussianBeamletSource" begin
        θ = 0.1 # radians
        bg = SphericalGaussianBeamletSource(pos, dir, θ, λ; num_rings = 5, num_rays = 500)
        @test bg isa AstigmaticBeamGroup

        @test length(bg) == 500
        # Check they all originate from pos
        @test all(BMO.position(b) ≈ pos for b in bg)
    end

    @testset "WavefrontBeamletDecomposition" begin
        # Generate a synthetic wavefront (a tilted plane wave)
        nx, ny = 11, 11
        x = LinRange(-1mm, 1mm, nx)
        y = LinRange(-1mm, 1mm, ny)

        amp = ones(nx, ny)
        # Tilt in x by angle alpha
        alpha = 0.001
        k = 2π / λ
        phase = [k * xi * alpha for xi in x, yi in y]

        bg = WavefrontBeamletDecomposition(x, y, amp, phase, dir, λ)
        @test bg isa AstigmaticBeamGroup

        @test length(bg) == nx * ny

        # We check the direction of the beamlets
        # The true local propagation direction should be tilted by alpha
        for b in bg
            d = normalize(BMO.direction(b))
            # Projection onto the original dir should be approximately cos(alpha)
            @test dot(d, dir)≈cos(alpha) atol=1e-3
        end
    end

    @testset "WavefrontBeamletDecomposition Energy Conservation" begin
        # Define a simple Gaussian wavefront
        λ = 1e-6
        w0_input = 1.0e-3 # 1 mm
        k = 2π / λ

        # Grid
        L = 5.0e-3 # 5 mm total width
        N = 41
        x = LinRange(-L / 2, L / 2, N)
        y = LinRange(-L / 2, L / 2, N)
        dx = x[2] - x[1]
        dy = y[2] - y[1]

        X = [xi for xi in x, yi in y]
        Y = [yi for xi in x, yi in y]
        R2 = X .^ 2 .+ Y .^ 2

        # Input field: A = exp(-r^2/w0^2)
        # Total power P = integral |A|^2 dA / (2Z) = pi * w0^2 / 4Z
        amplitude = exp.(-R2 / w0_input^2)
        phase = zeros(size(amplitude)) # Flat wavefront

        dir = [0.0, 0.0, 1.0]

        # Decomposition
        overlap = 1.2
        w0_beamlet = dx * overlap

        # The normalization in the code is for AMPLITUDE reconstruction.
        # E0 = A * (dx * dy) / (pi * w0s^2)
        # Let's see if this leads to power conservation.

        bg = WavefrontBeamletDecomposition(
            x, y, amplitude, phase, dir, λ; overlap = overlap, threshold = 1e-6)

        # Total power of beamlets (sum of individual powers)
        P_beamlets_sum = sum(optical_power(b) for b in bg.beams)

        # Total power of input field (discrete integral)
        # Using Z_vacuum to match library defaults
        P_input = sum(abs2.(amplitude)) * (dx * dy) / (2 * BeamletOptics.Z_vacuum)

        if length(bg.beams) > 0
            b1 = bg.beams[1]
            params1 = BeamletOptics.gauss_parameters(b1, 0.0)
        end

        # Reconstructed field power at z=0
        # Sample on a finer grid to check field reconstruction
        N_fine = 81
        x_fine = LinRange(-L / 2, L / 2, N_fine)
        y_fine = LinRange(-L / 2, L / 2, N_fine)
        dx_f = x_fine[2] - x_fine[1]
        dy_f = y_fine[2] - y_fine[1]

        E_rec = zeros(ComplexF64, N_fine, N_fine)
        for (ii, xf) in enumerate(x_fine)
            for (jj, yf) in enumerate(y_fine)
                r_vec = [xf, yf, 0.0]
                # Sum up contributions from all beamlets (coherent superposition)
                for b in bg.beams
                    p0 = BeamletOptics.position(b)
                    # For flat wavefront at z=0, the polarized field is appropriate
                    E_rec[ii, jj] += norm(BeamletOptics.polarized_field(b, r_vec - p0, 0.0))
                end
            end
        end

        P_rec = sum(abs2.(E_rec)) * (dx_f * dy_f) / (2 * BeamletOptics.Z_vacuum)

        # The coherent reconstruction should ideally match the input power
        @test isapprox(P_rec / P_input, 1.0, atol = 0.05)
    end
end

end # MODULE
