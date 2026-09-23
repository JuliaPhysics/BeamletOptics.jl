module TestFraunhofer

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics
const mm = 1e-3
const nm = 1e-9

@testset "Fraunhofer Diffraction (Astigmatic Gaussian Beamlets)" begin
    # D = 2 mm, λ = 632.8 nm, evaluated at 100 m
    D = 2mm
    λ = 632.8nm
    L = 100.0 # meters

    # We represent the square aperture by an array of 15 x 15 Gaussian beamlets for speed
    n_beams = 15
    Δ = D / n_beams
    w0s = Δ

    beams = CollimatedGaussianBeamletSource([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], D, λ, w0s; n_grid=n_beams)

    # Set electric field amplitude to point along z-axis to match previous test
    for b in beams
        BeamletOptics.polarization!(b.c.rays[1], [0.0, 0.0, 1.0])
    end

    # Physical detector large enough to catch all divergence rays
    pd = Detector(3000mm)
    translate_to3d!(pd, [0, L, 0])

    system = System([pd])

    # Trace all beams to the detector
    solve_system!.(Ref(system), beams)

    # All beams should have hit the detector
    @test length(pd.hits) == n_beams^2

    # Use a low resolution grid (n=50) for fast testing
    x, z, I = BMO.intensity(
        pd, n = 50, x_min = -100mm, x_max = 100mm, z_min = -100mm, z_max = 100mm)

    # Compare with analytical solution: I(x,z) ∝ sinc²(D*x/(λ*L)) * sinc²(D*z/(λ*L))
    x_norm = x .* D ./ (λ * L)
    z_norm = z .* D ./ (λ * L)
    sinc_sq(u) = (u == 0) ? 1.0 : (sin(π * u) / (π * u))^2
    analy_I = [sinc_sq(u_x) * sinc_sq(u_z) for u_x in x_norm, u_z in z_norm]

    # Normalize numerical intensity to match analytical peak (1.0)
    I_norm = I ./ maximum(I)

    # Calculate mean absolute error between normalized numerical and analytical intensities
    mae = sum(abs.(I_norm .- analy_I)) / length(I_norm)

    # With a 15x15 beamlet grid, the MAE should be very small (< 0.05)
    @test mae < 0.05

    # Check peak is exactly at the center
    mid_idx_x = argmin(abs.(x))
    mid_idx_z = argmin(abs.(z))
    @test I_norm[mid_idx_x, mid_idx_z] ≈ 1.0
end

end # MODULE
