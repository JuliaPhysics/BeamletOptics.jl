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
        bg = CollimatedGaussianBeamletSource(pos, dir, D, λ, w0s; n_grid=10)
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
        bg = GaussianBeamletDecomposition(pos, dir, w0, λ; n_grid=21)
        @test bg isa AstigmaticBeamGroup
        beams = BMO.beams(bg)
        
        @test length(beams) > 0
        @test length(beams) <= 441
        
        # Find the beamlet with the maximum amplitude
        center_beam = argmax([norm(BMO.polarization(b)) for b in bg])
        
        # The peak amplitude beamlet should be at the center (0, 0, 0)
        @test norm(BMO.position(bg[center_beam]) - pos) ≈ 0.0 atol=1e-5
    end

    @testset "SphericalGaussianBeamletSource" begin
        θ = 0.1 # radians
        bg = SphericalGaussianBeamletSource(pos, dir, θ, λ; num_rings=5, num_rays=500)
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
            @test dot(d, dir) ≈ cos(alpha) atol=1e-3
        end
    end
end

end # MODULE
