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

    @testset "SphericalGaussianBeamletSource" begin
        θ = 0.1 # radians
        bg = SphericalGaussianBeamletSource(pos, dir, θ, λ; num_rings = 5, num_rays = 500)
        @test bg isa AstigmaticBeamGroup

        @test length(bg) == 500
        # Check they all originate from pos
        @test all(BMO.position(b) ≈ pos for b in bg)

        @testset "basis kwarg" begin
            dir_n = normalize(dir)
            b0 = BMO.normal3d(dir_n)
            rotated = SphericalGaussianBeamletSource(pos, dir, θ, λ; num_rings = 5,
                num_rays = 500, basis = BMO.rotate3d(dir_n, deg2rad(30)) * b0)

            dirs = [BMO.direction(b) for b in bg]
            dirs_rot = [BMO.direction(b) for b in rotated]
            # spinning the basis about dir leaves the polar angle distribution intact
            @test sort(round.(BMO.angle3d.(Ref(dir_n), dirs), digits = 11)) ≈
                  sort(round.(BMO.angle3d.(Ref(dir_n), dirs_rot), digits = 11))
            @test !all(dirs .≈ dirs_rot)
            # passing the default basis reproduces the default
            @test all(dirs .≈ [BMO.direction(b) for b in SphericalGaussianBeamletSource(
                pos, dir, θ, λ; num_rings = 5, num_rays = 500, basis = b0)])
            # a degenerate basis would silently yield NaN directions, so it throws
            @test_throws ErrorException SphericalGaussianBeamletSource(pos, dir, θ, λ;
                num_rings = 5, num_rays = 500, basis = dir)
        end
    end

    @testset "EllipticalGaussianBeamletSource" begin
        θ_x = 0.2
        θ_y = 0.1
        bg = EllipticalGaussianBeamletSource(pos, dir, θ_x, θ_y, λ; num_rings = 5, num_rays = 500)
        @test bg isa AstigmaticBeamGroup

        @test length(bg) == 500
        # Check they all originate from pos
        @test all(BMO.position(b) ≈ pos for b in bg)

        # Check directions are within the elliptical cone
        dir_n = normalize(dir)
        b1 = BMO.normal3d(dir_n)
        b2 = cross(dir_n, b1)
        
        for b in bg
            d = BMO.direction(b)
            # project on tangent plane using local basis
            u = dot(d, b1) / dot(d, dir_n)
            v = dot(d, b2) / dot(d, dir_n)
            @test (u / tan(θ_x))^2 + (v / tan(θ_y))^2 <= 1.05
        end

        @testset "basis kwarg" begin
            b0 = BMO.normal3d(dir_n)
            # the elliptical cone is not rotationally symmetric, so a rotated basis
            # genuinely reorients the pattern rather than just relabelling it
            rotated = EllipticalGaussianBeamletSource(pos, dir, θ_x, θ_y, λ; num_rings = 5,
                num_rays = 500, basis = BMO.rotate3d(dir_n, deg2rad(30)) * b0)
            dirs = [BMO.direction(b) for b in bg]
            dirs_rot = [BMO.direction(b) for b in rotated]
            @test !all(dirs .≈ dirs_rot)
            # passing the default basis reproduces the default
            @test all(dirs .≈ [BMO.direction(b) for b in EllipticalGaussianBeamletSource(
                pos, dir, θ_x, θ_y, λ; num_rings = 5, num_rays = 500, basis = b0)])
            # a degenerate basis would silently yield NaN directions, so it throws
            @test_throws ErrorException EllipticalGaussianBeamletSource(pos, dir, θ_x, θ_y, λ;
                num_rings = 5, num_rays = 500, basis = dir)
        end
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

@testset "Power and profile of beamlet sources" begin
    # the power of a source on a detector that captures all of it
    function power(source, L, half_width)
        pd = Detector(2half_width)
        translate3d!(pd, [0, L, 0])
        solve_system!(System([pd]), source; progress = false)
        return BeamletOptics.optical_power(pd; progress = false, n = 200)
    end
    P0 = 1e-3
    @testset "GaussianBeamletDecomposition is the Gaussian beam of w0 and P0" begin
        λ, w0, L = 1e-6, 2e-3, 0.1
        for n_grid in (10, 20)
            g = GaussianBeamletDecomposition([0.0, 0, 0], [0.0, 1, 0], λ, w0; n_grid, P0)
            pd = Detector(12w0)
            translate3d!(pd, [0, L, 0])
            solve_system!(System([pd]), g; progress = false)
            x, z, I = intensity(pd; progress = false, n = 201, x_min = -3w0, x_max = 3w0, z_min = -3w0, z_max = 3w0)
            w = w0 * sqrt(1 + (L * λ / (π * w0^2))^2)
            @test BeamletOptics.trapz((x, z), I) ≈ P0 rtol = 1e-3
            # the peak intensity 2 P0 / (π w²) of the Gaussian beam, i.e. also its radius
            @test maximum(I) ≈ 2P0 / (π * w^2) rtol = 1e-3
        end
        # the sub-waist of the beamlets must be smaller than the waist
        @test_throws ArgumentError GaussianBeamletDecomposition([0.0, 0, 0], [0.0, 1, 0], λ, w0; n_grid = 4)
    end
    @testset "SphericalGaussianBeamletSource emits P0" begin
        for (θ, num_rings, num_rays) in ((0.02, 5, 200), (0.02, 10, 1000))
            source = SphericalGaussianBeamletSource([0.0, 0, 0], [0.0, 1, 0], θ, 1e-6; num_rings, num_rays, P0)
            @test power(source, 0.2, 3θ * 0.2) ≈ P0 rtol = 1e-2
        end
    end
    @testset "EllipticalGaussianBeamletSource emits P0" begin
        source = EllipticalGaussianBeamletSource([0.0, 0, 0], [0.0, 1, 0], 0.02, 0.01, 1e-6; num_rings = 10, num_rays = 1000, P0)
        @test power(source, 0.2, 3 * 0.02 * 0.2) ≈ P0 rtol = 1e-2
    end
    @testset "CollimatedGaussianBeamletSource emits P0" begin
        for n_grid in (5, 10)
            D = 10e-3
            source = CollimatedGaussianBeamletSource([0.0, 0, 0], [0.0, 1, 0], D, 1e-6, 1.2D / n_grid; n_grid, P0)
            @test power(source, 0.05, D) ≈ P0 rtol = 1e-2
        end
    end
    @testset "E0 sets the polarization, P0 the power" begin
        # E0 of any magnitude, linear or elliptical
        for E0 in ([3.0, 0, 0], [7.0, 0, 3.5im])
            spherical = SphericalGaussianBeamletSource([0.0, 0, 0], [0.0, 1, 0], 0.02, 1e-6; num_rings = 10, num_rays = 1000, P0, E0)
            @test power(spherical, 0.2, 0.012) ≈ P0 rtol = 1e-2
            elliptical = EllipticalGaussianBeamletSource([0.0, 0, 0], [0.0, 1, 0], 0.02, 0.01, 1e-6; num_rings = 10, num_rays = 1000, P0, E0)
            @test power(elliptical, 0.2, 0.012) ≈ P0 rtol = 1e-2
            decomposition = GaussianBeamletDecomposition([0.0, 0, 0], [0.0, 1, 0], 1e-6, 2e-3; P0, E0)
            @test power(decomposition, 0.1, 0.012) ≈ P0 rtol = 1e-2
            # the polarization of a beamlet is that of E0
            p = BeamletOptics.electric_field(first(BeamletOptics.beams(spherical)))
            @test p / norm(p) ≈ E0 / norm(E0)
        end
    end
    @testset "num_rings counts the center as a ring" begin
        # the center and num_rings - 1 rings around it, as for PointSource
        dir = [0.0, 1, 0]
        spherical = SphericalGaussianBeamletSource([0.0, 0, 0], dir, 0.02, 1e-6; num_rings = 5, num_rays = 500)
        angles = unique(round.([acos(clamp(dot(BeamletOptics.direction(b), dir), -1, 1)) for b in BeamletOptics.beams(spherical)], digits = 6))
        @test length(angles) == 5
        θx, θy = 0.02, 0.01
        elliptical = EllipticalGaussianBeamletSource([0.0, 0, 0], dir, θx, θy, 1e-6; num_rings = 5, num_rays = 500)
        b1 = BeamletOptics.orientation(elliptical)[:, 1]
        b2 = normalize(cross(dir, b1))
        ρ(d) = (u = dot(d, b1) / dot(d, dir); v = dot(d, b2) / dot(d, dir); sqrt((u / tan(θx))^2 + (v / tan(θy))^2))
        radii = unique(round.([ρ(BeamletOptics.direction(b)) for b in BeamletOptics.beams(elliptical)], digits = 6))
        @test sort(radii) ≈ collect(LinRange(0, 1, 5))
    end
end

end # MODULE
