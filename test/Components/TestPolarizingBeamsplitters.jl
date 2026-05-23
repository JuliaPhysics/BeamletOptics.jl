module TestPolarizingBeamsplitters

using BeamletOptics
using Test
using LinearAlgebra
using AbstractTrees

const BMO = BeamletOptics
const mm = 1e-3

@testset "Polarizing Beamsplitters" begin
    @testset "RectangularPolarizingPlateBeamsplitter (n > 1)" begin
        # Substrate n=1.5
        n_val = 1.5
        n(λ) = n_val
        width = 10mm
        height = 10mm
        thickness = 2mm

        # Create PPBS
        ppbs = RectangularPolarizingPlateBeamsplitter(width, height, thickness, n)
        # It is aligned with negative y-axis. Coating is at origin (after translation in constructor?)
        # Logic in constructor:
        # substrate translated [0, thickness/2, 0]
        # coating interact sets T=ray, R=reflection
        # interact3d refracts T into substrate (entering) -> n=1 to n=1.5

        system = StaticSystem([ppbs])

        # Ray coming from -y direction (to hit coating at y=0 first)
        # Normal of coating (rotated pi) is [0, -1, 0].
        # Ray dir should represent incidence from outside.

        angle_deg = 15
        # Coming from -y, moving +y. Tilted in +x.
        # Dir = [sin, cos, 0]
        ki = normalize([sind(angle_deg), cosd(angle_deg), 0])

        # P-pol (in plane of incidence, xy-plane) -> should be transmitted
        # Parallel to plane. Orthogonal to ki.
        # ki = [sin, cos, 0].
        # p = [cos, -sin, 0]. dot = sc - cs = 0.
        E0_p = normalize([cosd(angle_deg), -sind(angle_deg), 0])

        # Start at y = -10mm
        ray_p = PolarizedRay([0, -10mm, 0], ki, 1000e-9, E0_p)
        beam_p = Beam(ray_p)
        solve_system!(system, beam_p)

        trans_ray_p = first(BMO.rays(beam_p.children[1]))
        refl_ray_p = first(BMO.rays(beam_p.children[2]))

        # Snell: n1=1 (air), n2=1.5 (substrate).
        # sin(th1) = sind(angle_deg).
        # sin(th2) = sin(th1) / 1.5.
        th2 = asin(sind(angle_deg) / 1.5)
        expected_dir = [sin(th2), cos(th2), 0]
        @test BMO.direction(trans_ray_p)≈expected_dir atol=1e-6

        # Check orthogonality (THIS WAS THE BUG)
        @test abs(dot(BMO.direction(trans_ray_p), BMO.polarization(trans_ray_p))) < 1e-9

        # Power Conservation Test
        # I_in = I_reflected_front + I_transmitted + I_reflected
        # But our system just treats entering substrate:
        # Beam hits front face -> trans + refl.
        # Trans hits coating -> trans + refl (internal).

        # Let's just trace the whole system and sum up the intensity of all leaves
        leaves = collect(AbstractTrees.Leaves(beam_p))
        total_power_out = sum(BMO.refractive_index(last(BMO.rays(leaf))) *
                              norm(BMO.polarization(last(BMO.rays(leaf))))^2
        for leaf in leaves)
        power_in = BMO.refractive_index(ray_p) * norm(E0_p)^2

        # With n=1.5, we expect Fresnel losses at entry/exit.
        # Total power out should be approx 0.9216 * power_in if both branches exit.
        # But in this test, one branch reflects from air-coating interface.
        @test total_power_out < power_in
        @test total_power_out > 0.9 * power_in
    end

    @testset "PolarizingCubeBeamsplitter (n > 1)" begin
        n_val = 1.5
        n(λ) = n_val
        L = 20mm

        pcbs = PolarizingCubeBeamsplitter(L, n)
        # Cube centered at origin.
        # Coating at 45 deg.
        # Front face at y = L/2 ?
        # RightAnglePrism front: legs L, L.
        # Geometry roughly: usually input from -x or -y or +x etc.
        # Standard interactions usually assumes input along optical axis.

        # Let's try simple on-axis input (normal incidence to front face)
        # If normal incidence, no refraction direction change, so defects might be hidden.
        # We need off-axis input to trigger refraction direction change @ front face.

        system = StaticSystem([pcbs])

        # Front face is likely at some position.
        # Let's shoot ray from far away at angle.
        angle_deg = 10
        ki = normalize([sind(angle_deg), -cosd(angle_deg), 0])
        # Assume front face is facing +y (plane at y=L/2)
        # Intersect front face.

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

        # Check Pi phase jump for reflection
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

        # Reflected ray should be along -x or +x.
        # Orientation of PCBS coating is 135 deg (normal along [1,1,0]/sqrt(2)).
        # Reflected dir should be [-1, 0, 0].
        refl_rays = filter(l -> abs(BMO.direction(last(BMO.rays(l)))[1]) > 0.9, leaves)
        @test !isempty(refl_rays)
        refl_ray = last(BMO.rays(first(refl_rays)))

        # Check that Ez is flipped
        # With n=1.5, amplitude is -0.96
        @test real(BMO.polarization(refl_ray)[3])≈-0.96 atol=1e-6
    end

    @testset "AstigmaticGaussianBeamlet with Polarizing Beamsplitters" begin
        # 1. Thin PolarizingBeamSplitter
        pbs = PolarizingBeamSplitter(20mm, 20mm)
        system1 = System([pbs])
        agb1 = AstigmaticGaussianBeamlet([0.0, -50mm, 0.0], [0.0, 1.0, 0.0], 1000e-9, 1mm;
            E0=[1.0, 0.0, 1.0]/sqrt(2), support=[1.0, 0.0, 0.0])
        solve_system!(system1, agb1)
        @test length(children(agb1)) == 2
        @test isa(children(agb1)[1], AstigmaticGaussianBeamlet)
        @test isa(children(agb1)[2], AstigmaticGaussianBeamlet)
        
        # 2. PolarizingCubeBeamsplitter
        pcbs = PolarizingCubeBeamsplitter(20mm, n -> 1.0)
        system2 = System([pcbs])
        agb2 = AstigmaticGaussianBeamlet([0.0, -50mm, 0.0], [0.0, 1.0, 0.0], 1000e-9, 1mm;
            E0=[1.0, 0.0, 1.0]/sqrt(2), support=[1.0, 0.0, 0.0])
        solve_system!(system2, agb2)
        @test length(children(agb2)) == 2

        # 3. RectangularPolarizingPlateBeamsplitter
        ppbs = RectangularPolarizingPlateBeamsplitter(20mm, 20mm, 5mm, n -> 1.0)
        system3 = System([ppbs])
        agb3 = AstigmaticGaussianBeamlet([0.0, -50mm, 0.0], [0.0, 1.0, 0.0], 1000e-9, 1mm;
            E0=[1.0, 0.0, 1.0]/sqrt(2), support=[1.0, 0.0, 0.0])
        solve_system!(system3, agb3)
        @test length(children(agb3)) == 2
    end
end
end
