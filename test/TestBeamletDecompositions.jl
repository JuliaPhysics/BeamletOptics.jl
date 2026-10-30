module TestBeamletDecompositions

using BeamletOptics
using Test
using LinearAlgebra
using Random

const BMO = BeamletOptics
const mm = 1e-3
const nm = 1e-9

@testset "Beamlet decompositions" begin
    λ = 532nm
    dir = [0.0, 1.0, 0.0]
    pos = [0.0, 0.0, 0.0]

    @testset "GaussianBeamletDecomposition" begin
        w0 = 1mm
        # A 21x21 grid will trim beamlets with amp < 1e-4, so length will be < 441
        bg = GaussianBeamletDecomposition(pos, dir, λ, w0; n_grid = 21)
        @test bg isa AstigmaticBeamGroup
        beams = BMO.beams(bg)

        @test length(beams) > 0
        @test length(beams) <= 441

        # Find the beamlet with the maximum amplitude
        center_beam = argmax([norm(BMO.polarization(b)) for b in bg])

        # The peak amplitude beamlet should be at the center (0, 0, 0)
        @test norm(BMO.position(bg[center_beam]) - pos)≈0.0 atol=1e-5
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

        bg = WavefrontBeamletDecomposition(x, y, amp, phase, pos, dir, λ)
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
            x, y, amplitude, phase, zeros(3), dir, λ; overlap = overlap, threshold = 1e-6)

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

@testset "Power and profile of beamlet decompositions" begin
    P0 = 1e-3
    @testset "GaussianBeamletDecomposition is the Gaussian beam of w0 and P0" begin
        λ, w0, L = 1e-6, 2e-3, 0.1
        for n_grid in (10, 20)
            g = GaussianBeamletDecomposition([0.0, 0, 0], [0.0, 1, 0], λ, w0; n_grid, P0)
            pd = Detector(12w0)
            translate3d!(pd, [0, L, 0])
            solve_system!(System([pd]), g; progress = false)
            s = RegularGrid(num_points = 201, x_min = -3w0, x_max = 3w0, z_min = -3w0, z_max = 3w0)
            x, z, I = intensity(pd, s; progress = false)
            w = w0 * sqrt(1 + (L * λ / (π * w0^2))^2)
            @test optical_power(pd, s; progress = false) ≈ P0 rtol = 1e-3
            # the peak intensity 2 P0 / (π w²) of the Gaussian beam, i.e. also its radius
            @test maximum(I) ≈ 2P0 / (π * w^2) rtol = 1e-3
        end
        # the sub-waist of the beamlets must be smaller than the waist
        @test_throws ArgumentError GaussianBeamletDecomposition([0.0, 0, 0], [0.0, 1, 0], λ, w0; n_grid = 4)
    end
end

@testset "Grid position of the decompositions" begin
    λ = 1e-6
    p = [1e-3, -2e-3, 5e-3]
    dir = [0.3, 1.0, -0.2]
    x = collect(LinRange(-1e-3, 1e-3, 9))
    y = collect(LinRange(-0.8e-3, 0.8e-3, 7))
    amp = [exp(-(xi^2 + yi^2) / 0.7e-3^2) for xi in x, yi in y]
    phase = [2π / λ * (2e-4 * xi - 1e-4 * yi) for xi in x, yi in y]
    Eu = amp .* cis.(phase)
    Ev = 0.3im .* Eu
    beam_positions(g) = [BMO.position(BMO.first_ray(b)) for b in BMO.beams(g)]
    # the threads of a decomposition push the beamlets in any order
    same_positions(a, b) = length(a) == length(b) && all(any(isapprox(q, r; rtol = 1e-12) for r in b) for q in a)

    @testset "$name" for (name, F1, F2) in (("scalar", amp, phase), ("vector", Eu, Ev))
        g = WavefrontBeamletDecomposition(x, y, F1, F2, p, dir, λ)
        @test BMO.position(g) == p
        ref = WavefrontBeamletDecomposition(x, y, F1, F2, zeros(3), dir, λ)
        translate3d!(ref, p)
        @test same_positions(beam_positions(g), beam_positions(ref))
    end

    @testset "Float32 grid, Float64 position" begin
        g = WavefrontBeamletDecomposition(Float32.(x), Float32.(y), Float32.(amp),
            Float32.(phase), p, Float32.(dir), 1.0f-6)
        @test g isa AstigmaticBeamGroup{Float64}
        @test BMO.position(g) == p
        # a Float32 direction that fails without the conversion of dir to Float64
        @test WavefrontBeamletDecomposition(Float32.(x), Float32.(y), Float32.(amp), Float32.(phase),
            zeros(3), Float32[-0.75, -0.75, -0.25], 1.0f-6) isa AstigmaticBeamGroup{Float64}
    end

    @testset "Rational phase" begin
        g = WavefrontBeamletDecomposition(x, y, amp, rationalize.(phase ./ 1e3; tol = 1e-6), p, dir, λ)
        @test length(BMO.beams(g)) > 0
    end

    @testset "GaussianModeDecomposition" begin
        xs = collect(LinRange(-1.2e-3, 1.2e-3, 25))
        ys = collect(LinRange(-1e-3, 1e-3, 21))
        xc, yc = 1e-4, -5e-5
        ψ = [complex(exp(-((xi - xc)^2 + (yi - yc)^2) / 0.4e-3^2)) for xi in xs, yi in ys]
        e1, e2, n = [1.0, 0, 0], [0.0, 0, -1], [0.0, 1, 0]
        m = GaussianModeDecomposition(xs, ys, ψ, zero(ψ), p, n, λ; basis = (e1, e2))
        @test BMO.position(BMO.first_ray(m)) ≈ p + xc * e1 + yc * e2 atol = 1e-9
    end
end

@testset "Decomposition helpers" begin
    dir = normalize([0.3, 1.0, -0.2])
    e1, e2 = BMO.grid_axes(dir, nothing)
    @test norm(e1) ≈ 1 && norm(e2) ≈ 1
    @test abs(dot(e1, e2)) < 1e-15 && abs(dot(e1, dir)) < 1e-15 && abs(dot(e2, dir)) < 1e-15
    a, b = [2.0, 0, 0], [0.0, 0, -3]
    @test BMO.grid_axes(dir, (a, b)) == (normalize(a), normalize(b))

    s = BMO.random_support(dir, Random.Xoshiro(1))
    @test norm(s) ≈ 1
    @test abs(dot(s, dir)) < 1e-15

    @test BMO.sphere_phase(1e-3, -2e-3, Inf, 2π / 1e-6) == 0
end

end # MODULE
