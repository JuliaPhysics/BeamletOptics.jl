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
        # It should traverse Front (refract) -> Coating (split) -> Back (refract)
        # We need to find the final transmitted beam.
        # Beam tree: Root -> [Front Interaction] -> [Coating Interaction] -> [Back Interaction] -> Exit?
        # Or Root -> Front -> Coating -> Back?
        # "interact3d" returns a beam system (children).

        # We can just collect all final rays (leaves of beam tree)

        # Wait, interact3d logic for Cube:
        # if Front hit -> interact(front), hint(coating). Returns interaction (new ray).
        # That new ray is child of root.
        # Then System continues tracing child.
        # Child hits Coating -> interact(coating). Returns T and R rays.
        # T ray hits Back -> interact(back). Returns Exit ray.

        # So leaf rays are the ones.
        leaves = collect(AbstractTrees.Leaves(beam))
        # We expect at least one transmitted ray that exited the cube.

        # Filter for transmitted one (roughly same direction as input, but refracted twice -> should be parallel to input if faces parallel)
        # Cube faces are parallel. Plate/Cube with parallel faces -> output parallel to input.
        # So final direction ≈ input direction.

        for leaf in leaves
            r = first(BMO.rays(leaf))
            # Check orthogonality
            @test abs(dot(BMO.direction(r), BMO.polarization(r))) < 1e-9
        end
    end
end
end
