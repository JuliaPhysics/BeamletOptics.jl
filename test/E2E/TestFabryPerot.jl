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

end # module
