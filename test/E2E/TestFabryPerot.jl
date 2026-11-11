module TestFabryPerot

using BeamletOptics
using Test
using LinearAlgebra
using AbstractTrees

const BMO = BeamletOptics
const mm = 1e-3
const μm = 1e-6

@testset "Beamlets split at a coated face" begin
    plate = RectangularCompensatorPlate(20mm, 20mm, 2mm, ConstantRefractiveIndex(1.5))
    surface_model!(plate, FlatFace([0, -1, 0]), BeamsplitterCoating(0.3))
    zrotate3d!(plate, deg2rad(10))
    translate3d!(plate, [0, 20mm, 0])
    for beamlet in (GaussianBeamlet([0, 0, 0], [0, 1.0, 0], 1e-6, 0.5mm),
            AstigmaticGaussianBeamlet([0, 0, 0], [0, 1.0, 0], 1e-6, 0.5mm, 0.3mm))
        solve_system!(System([plate]), beamlet; depth_max = 1)
        @test length(children(beamlet)) == 2
        @test optical_power(children(beamlet)[2]) ≈ 0.3 * optical_power(beamlet)
        @test sum(optical_power, children(beamlet)) ≈ optical_power(beamlet)
    end
end

@testset "Fabry-Pérot etalon: Airy transmission" begin
    # A glass etalon of 10 µm with lossless 80 % mirrors on both faces, Born & Wolf §7.6
    n, d, R = 1.5, 10μm, 0.8
    etalon = SphericalLens(Inf, Inf, d, 2mm, ConstantRefractiveIndex(n))
    surface_model!(etalon, FrontFace(), BeamsplitterCoating(R))
    surface_model!(etalon, BackFace(), BeamsplitterCoating(R))
    detector = Detector(1mm)
    translate3d!(detector, [0, 1mm, 0])
    system = System([etalon, detector])
    airy(λ) = 1 / (1 + 4R / (1 - R)^2 * sin(2π * n * d / λ)^2)
    for λ in (620e-9, 622.5e-9, 625e-9, 627.5e-9)
        empty!(detector)
        beamlet = AstigmaticGaussianBeamlet([0, -1mm, 0], [0, 1.0, 0], λ, 100μm, 100μm; P0 = 1.0)
        # 40 splits are 20 round trips: R⁴⁰ ≈ 1e-4 of the field is cut off
        solve_system!(system, beamlet; depth_max = 40)
        @test optical_power(detector) ≈ airy(λ) atol = 0.03
    end
end

@testset "Fizeau wedge: fringe period λ/(2α)" begin
    R, α, λ = 0.8, deg2rad(0.1), 632.8e-9
    function mirror_plate(inner)
        plate = RectangularCompensatorPlate(5mm, 5mm, 0.5mm, ConstantRefractiveIndex(1.5))
        surface_model!(plate, FlatFace([0, inner, 0]), BeamsplitterCoating(R))
        surface_model!(plate, FlatFace([0, -inner, 0]), AntireflectionCoating())
        return plate
    end
    p1, p2 = mirror_plate(1), mirror_plate(-1)
    translate3d!(p1, [0, -0.5mm, 0])
    zrotate3d!(p1, α / 2)
    translate3d!(p2, [0, 20μm, 0])
    zrotate3d!(p2, -α / 2)
    detector = Detector(1.5mm)
    translate3d!(detector, [0, 2mm, 0])
    beamlet = AstigmaticGaussianBeamlet([0, -2mm, 0], [0, 1.0, 0], λ, 300μm, 300μm; P0 = 1.0)
    solve_system!(System([p1, p2, detector]), beamlet; depth_max = 40)
    x, _, I = intensity(detector, RegularGrid(num_points = 400))
    profile = vec(sum(I; dims = 2))
    peaks = [x[i] for i in 2:(length(profile) - 1)
             if profile[i] > profile[i - 1] && profile[i] ≥ profile[i + 1] && profile[i] > 0.2maximum(profile)]
    @test length(peaks) ≥ 3
    @test sum(diff(peaks)) / (length(peaks) - 1) ≈ λ / (2α) rtol = 0.05
end

# The field of the k-th echo of a tilted plate (thickness t, index n, angle θ) on a screen at the
# distance zdet from the waist w0, as a Gaussian beam: shifted by k s across the beam, delayed by
# k 2nt cos θ′, at the reduced distance zdet - t + (2k + 1) t / n and weakened by R^k
function echo_field(x, k; λ, n, t, θ, w0, zdet, x0, side)
    θ′ = asin(sin(θ) / n)
    s = 2t * cos(θ) * tan(θ′)
    R = (((cos(θ) - n * cos(θ′)) / (cos(θ) + n * cos(θ′)))^2 +
         ((n * cos(θ) - cos(θ′)) / (n * cos(θ) + cos(θ′)))^2) / 2
    zR = π * w0^2 / λ
    D = zdet - t + (2k + 1) * t / n
    w, Rc, ψ = w0 * sqrt(1 + (D / zR)^2), D * (1 + (zR / D)^2), atan(D / zR)
    xk = x0 + k * side * s
    phase = 2π / λ * (k * 2n * t * cos(θ′) + (x - xk)^2 / (2Rc)) - ψ
    return R^k * (w0 / w) * exp(-(x - xk)^2 / w^2) * cis(phase)
end

@testset "Ghost fringes of a tilted window" begin
    # A beamlet through an uncoated window, 10 mm thick and tilted by 2°: its echoes interfere with
    # the beam on the screen. Profile, offsets and phases agree with Gaussian beams of the echoes.
    λ, n, t, θ, w0, zdet = 632.8e-9, 1.515, 10mm, deg2rad(2), 20μm, 150mm
    θ′ = asin(sin(θ) / n)
    xs = collect(range(-3mm, 3mm; length = 601))
    # The setup as built and moved: placement does not matter
    for (shift, φ) in (([0, 0, 0], 0.0), ([12mm, -30mm, 5mm], deg2rad(35)))
        window = Lens(CircularFlatSurface(25mm), CircularFlatSurface(25mm), t, ConstantRefractiveIndex(n))
        zrotate3d!(window, θ)
        translate3d!(window, [0, 50mm, 0])
        screen = Detector(8mm)
        translate3d!(screen, [0, zdet, 0])
        beam = GaussianBeamlet([0, 0, 0], [0, 1, 0], λ, w0)
        for part in (window, screen, beam)
            rotate3d!(part, [0, 0, 1], φ, [0, 0, 0])
            translate3d!(part, shift)
        end
        solve_system!(System([window, screen]; surface_model = FresnelGhosts()), beam; depth_max = 6)
        # Beam and two echoes reach the screen, the first one shifted by the plate
        pts, _ = ray_hits(screen)
        x = sort(first.(pts); by = abs)
        @test length(x) == 3
        @test abs(x[1]) ≈ t * (sin(θ) - cos(θ) * tan(θ′)) rtol = 1e-6
        @test abs(x[2] - x[1]) ≈ 2t * cos(θ) * tan(θ′) rtol = 1e-6
        @test abs(x[3] - x[1]) ≈ 4t * cos(θ) * tan(θ′) rtol = 1e-6
        I = vec(last(intensity(screen, RectilinearGrid(xs, [0.0]))))
        I_echoes = [abs2(sum(echo_field(xi, k; λ, n, t, θ, w0, zdet, x0 = x[1], side = sign(x[2] - x[1]))
                    for k in 0:2)) for xi in xs]
        I_echoes .*= dot(I_echoes, I) / dot(I_echoes, I_echoes)
        # Up to the brightness, the fringes agree to 1e-3 of the peak, against a modulation of 7 %
        @test maximum(abs.(I .- I_echoes)) < 1e-3 * maximum(I)
        I_beam = [abs2(echo_field(xi, 0; λ, n, t, θ, w0, zdet, x0 = x[1], side = 1)) for xi in xs]
        I_beam .*= dot(I_beam, I) / dot(I_beam, I_beam)
        @test maximum(abs.(I .- I_beam)) > 0.05 * maximum(I)
    end
end

end # module
