module TestPolarizingBeamsplitters

using BeamletOptics
using Test
using LinearAlgebra
using AbstractTrees

const BMO = BeamletOptics
const mm = 1e-3

@testset "Polarizing Beamsplitters" begin
    @testset "Polarizing PlateBeamsplitter Constructor Type Promotion" begin
        n(λ) = 1.5
        # Passing integer dimensions to RectangularPolarizingPlateBeamsplitter
        rppbs = RectangularPolarizingPlateBeamsplitter(10, 10, 2, n)
        @test rppbs isa RectangularPolarizingPlateBeamsplitter{Float64}
        @test BMO.shape(rppbs) isa BMO.BoxSDF{Float64}
        @test BMO.coatings(rppbs)[1].second isa BMO.JonesCoating

        # Passing integer dimensions to RoundPolarizingPlateBeamsplitter
        round_ppbs = RoundPolarizingPlateBeamsplitter(20, 2, n)
        @test round_ppbs isa RoundPolarizingPlateBeamsplitter{Float64}
        @test BMO.shape(round_ppbs) isa BMO.PlanoSurfaceSDF{Float64}
        @test BMO.coatings(round_ppbs)[1].second isa BMO.JonesCoating
    end

    @testset "RectangularPolarizingPlateBeamsplitter (n > 1)" begin
        # Substrate n=1.5
        n_val = 1.5
        n(λ) = n_val
        width = 10mm
        height = 10mm
        thickness = 2mm

        # Create PPBS
        ppbs = RectangularPolarizingPlateBeamsplitter(width, height, thickness, n)

        system = StaticSystem([ppbs])

        angle_deg = 15
        # Coming from -y, moving +y. Tilted in +x.
        ki = normalize([sind(angle_deg), cosd(angle_deg), 0])

        # P-pol (in plane of incidence, xy-plane) -> should be transmitted
        E0_p = normalize([cosd(angle_deg), -sind(angle_deg), 0])

        # Start at y = -10mm
        ray_p = PolarizedRay([0, -10mm, 0], ki, 1000e-9, E0_p)
        beam_p = Beam(ray_p)
        solve_system!(system, beam_p)

        trans_ray_p = first(BMO.rays(beam_p.children[1]))
        refl_ray_p = first(BMO.rays(beam_p.children[2]))

        th2 = asin(sind(angle_deg) / 1.5)
        expected_dir = [sin(th2), cos(th2), 0]
        @test BMO.direction(trans_ray_p)≈expected_dir atol=1e-6

        # Check orthogonality
        @test abs(dot(BMO.direction(trans_ray_p), BMO.polarization(trans_ray_p))) < 1e-9

        # Power Conservation Test
        leaves = collect(AbstractTrees.Leaves(beam_p))
        total_power_out = sum(BMO.refractive_index(last(BMO.rays(leaf))) *
                              norm(BMO.polarization(last(BMO.rays(leaf))))^2
        for leaf in leaves)
        power_in = BMO.refractive_index(ray_p) * norm(E0_p)^2

        @test total_power_out < power_in
        @test total_power_out > 0.9 * power_in
    end

    @testset "PolarizingCubeBeamsplitter (n > 1)" begin
        n_val = 1.5
        n(λ) = n_val
        L = 20mm

        pcbs = PolarizingCubeBeamsplitter(L, n)

        system = StaticSystem([pcbs])

        angle_deg = 10
        ki = normalize([sind(angle_deg), -cosd(angle_deg), 0])

        # P-pol (in plane)
        E0_p = normalize([-ki[2], ki[1], 0])

        ray = PolarizedRay([0, 2 * L, 0], ki, 1000e-9, E0_p)
        beam = Beam(ray)
        solve_system!(system, beam)

        # Check children (Transmission)
        leaves = collect(AbstractTrees.Leaves(beam))
        total_power_out = sum(BMO.refractive_index(last(BMO.rays(leaf))) *
                              norm(BMO.polarization(last(BMO.rays(leaf))))^2
        for leaf in leaves)
        power_in = BMO.refractive_index(ray) * norm(E0_p)^2
        @test total_power_out < power_in
        @test total_power_out > 0.9 * power_in
    end

    @testset "Pi phase jump for PCBS reflection" begin
        pcbs = PolarizingCubeBeamsplitter(20mm, n -> 1.5)
        system = StaticSystem([pcbs])

        # Ray along y, z-polarized
        ki = [0.0, 1.0, 0.0]
        E0_z = [0.0, 0.0, 1.0]
        ray = PolarizedRay([0, -20mm, 0], ki, 1000e-9, E0_z)
        beam = Beam(ray)
        solve_system!(system, beam)

        leaves = collect(AbstractTrees.Leaves(beam))

        refl_rays = filter(l -> abs(BMO.direction(last(BMO.rays(l)))[1]) > 0.9, leaves)
        @test !isempty(refl_rays)
        refl_ray = last(BMO.rays(first(refl_rays)))

        # Check that Ez is flipped
        @test real(BMO.polarization(refl_ray)[3])≈0.96 atol=1e-6
    end

    @testset "AstigmaticGaussianBeamlet with Polarizing Beamsplitters" begin
        # Thin PolarizingBeamSplitter
        pbs = PolarizingBeamSplitter(20mm, 20mm)
        system1 = System([pbs])
        agb1 = AstigmaticGaussianBeamlet([0.0, -50mm, 0.0], [0.0, 1.0, 0.0], 1000e-9, 1mm;
            E0 = [1.0, 0.0, 1.0] / sqrt(2), support = [1.0, 0.0, 0.0])
        solve_system!(system1, agb1)
        @test length(children(agb1)) == 2
        @test isa(children(agb1)[1], AstigmaticGaussianBeamlet)
        @test isa(children(agb1)[2], AstigmaticGaussianBeamlet)

        # PolarizingCubeBeamsplitter
        pcbs = PolarizingCubeBeamsplitter(20mm, n -> 1.0)
        system2 = System([pcbs])
        agb2 = AstigmaticGaussianBeamlet([0.0, -50mm, 0.0], [0.0, 1.0, 0.0], 1000e-9, 1mm;
            E0 = [1.0, 0.0, 1.0] / sqrt(2), support = [1.0, 0.0, 0.0])
        solve_system!(system2, agb2)
        @test length(children(agb2)) == 2

        # RectangularPolarizingPlateBeamsplitter
        ppbs = RectangularPolarizingPlateBeamsplitter(20mm, 20mm, 5mm, n -> 1.0)
        system3 = System([ppbs])
        agb3 = AstigmaticGaussianBeamlet([0.0, -50mm, 0.0], [0.0, 1.0, 0.0], 1000e-9, 1mm;
            E0 = [1.0, 0.0, 1.0] / sqrt(2), support = [1.0, 0.0, 0.0])
        solve_system!(system3, agb3)
        @test length(children(agb3)) == 2
    end
end
end
