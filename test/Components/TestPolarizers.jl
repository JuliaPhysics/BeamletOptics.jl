module TestPolarizers

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics

const mm = 1e-3

@testset "Polarizing components" begin
    # Test fcts.
    pseudo_I(v) = norm(v)^2
    pseudo_I(ray::BMO.PolarizedRay) = pseudo_I(BMO.polarization(ray))
    malus_law(θ) = cosd(θ)^2
    function generate_linpol_angles(thetas)
        angles = zeros(length(thetas))
        # consider flip in angle quadrant above 90°
        for i in eachindex(thetas)
            theta = thetas[i] - 1
            if theta > 90 && theta <= 270
                angles[i] = 180
            end
        end
        return angles
    end

    @testset "Jones polarization matrices" begin
        @test isdefined(BMO, :AbstractJonesMatrix)
        @test isdefined(BMO, :LocalJonesBasis)
        @test isdefined(BMO, :GlobalJonesBasis)
        # test matrix operations
        J = BMO.XYBasis(1, 2, 3, 4)
        R = [1 2 0; 3 4 0; 0 0 1]
        @test J == R
        @test J * 2 == R * 2
        @test transpose(J) == transpose(R)
        @test inv(J) ≈ inv(R)
    end

    @testset "Polarizing filter" begin
        filter = PolarizationFilter(5mm)
        system = System([filter])
        # rotation angle steps, intensity and angle comparison
        thetas = 1:10:360
        i_n = zeros(length(thetas))
        i_a = similar(i_n)
        angles_n = similar(i_n)
        angles_a = generate_linpol_angles(thetas)

        @testset "Pol. filter normal incidence - rotation around y-axis" begin
            Rfil = orientation(filter)
            ray_dir = Rfil[:, 2]
            ray_pos = position(filter) - 10mm * ray_dir
            local_x = Rfil[:, 1]
            local_z = Rfil[:, 3]
            pol_vec = local_x
            lambda = 1e-6
            beam = Beam(ray_pos, ray_dir, lambda, pol_vec)
            # step ref. vector by angle increment
            RotMat = BMO.rotate3d(ray_dir, deg2rad(step(thetas)))
            for i in eachindex(thetas)
                solve_system!(system, beam)
                E1 = BMO.polarization(BMO.rays(beam)[2])
                i_n[i] = pseudo_I(E1)
                i_a[i] = malus_law(thetas[i] - 1)
                rotate3d!(filter, ray_dir, deg2rad(step(thetas)))
                # test polarization state
                @test all(BMO.islinear.(beam.rays))
                # test projected polarization direction
                v1 = real(E1)
                angles_n[i] = rad2deg(BMO.angle3d(v1, pol_vec))
                # update ref vector
                pol_vec = RotMat * pol_vec
            end
            @test angles_n ≈ angles_a
        end

        # Move and rotate filter
        translate3d!(filter, [0, 10mm, 0])
        xrotate3d!(filter, deg2rad(45))
        zrotate3d!(filter, deg2rad(30))

        @testset "Pol. filter tilted incidence - rotation around ray optical axis" begin
            Rfil = orientation(filter)
            ray_dir = Rfil[:, 2]
            ray_pos = position(filter) - 10mm * ray_dir
            local_x = Rfil[:, 1]
            local_z = Rfil[:, 3]
            pol_vec = local_x
            lambda = 1e-6
            beam = Beam(ray_pos, ray_dir, lambda, pol_vec)
            # tilt filter
            θ_tilt = 45
            rotate3d!(filter, local_x, deg2rad(θ_tilt))
            i_n_tilt = similar(i_n)
            for i in eachindex(thetas)
                solve_system!(system, beam)
                i_n_tilt[i] = pseudo_I(beam.rays[2].E0)
                rotate3d!(filter, ray_dir, deg2rad(step(thetas)))
                @test all(BMO.islinear.(beam.rays))
            end
            # only applicable for θ_tilt = 45°
            i_a_tilt = malus_law.(thetas .- 1) .* 0.75 .+ 0.25
            @test i_n_tilt ≈ i_a_tilt
        end
    end

    @testset "Polarizing Cube Beamsplitter" begin
        # PBS cube of leg length 25mm, index n = 1.5
        pbs = PolarizingCubeBeamsplitter(25e-3, n -> 1.5)
        translate3d!(pbs, [0, 50mm, 0])
        system = System([pbs])

        @testset "s-polarization (reflected)" begin
            # s-polarization is along z-axis, perpendicular to the plane of incidence (x-y)
            beam = Beam([0, 0, 0], [0, 1, 0], 1e-6, [0, 0, 1])
            solve_system!(system, beam)

            # Transmitted beam (children[1]) should have 0 amplitude for s
            # Reflected beam (children[2]) should have full amplitude for s
            @test length(beam.children) == 2
            t_rays = BMO.rays(beam.children[1])
            r_rays = BMO.rays(beam.children[2])
            @test norm(BMO.polarization(last(t_rays)))≈0.0 atol=1e-12
            @test norm(BMO.polarization(last(r_rays)))≈0.96 atol=1e-12
            @test BMO.direction(last(r_rays)) ≈ [-1, 0, 0]
        end

        @testset "p-polarization (transmitted)" begin
            # p-polarization is in the x-y plane, along x-axis
            beam = Beam([0, 0, 0], [0, 1, 0], 1e-6, [1, 0, 0])
            solve_system!(system, beam)

            # Transmitted beam (children[1]) should have full amplitude for p (with interface loss)
            # Reflected beam (children[2]) should have 0 amplitude for p
            @test length(beam.children) == 2
            t_rays = BMO.rays(beam.children[1])
            r_rays = BMO.rays(beam.children[2])
            @test norm(BMO.polarization(last(t_rays)))≈0.96 atol=1e-12
            @test norm(BMO.polarization(last(r_rays)))≈0.0 atol=1e-12
            @test BMO.direction(last(t_rays)) ≈ [0, 1, 0]
        end
    end

    @testset "Waveplates" begin
        @testset "Flat Half-Wave Plate (round)" begin
            # Fast axis oriented at 22.5 degrees around y-axis.
            hwp = HalfWavePlate(5mm; fast_axis_angle = deg2rad(22.5))
            system = System([hwp])
            # Incident polarization along x-axis
            beam = Beam([0, -10mm, 0], [0, 1, 0], 1e-6, [1, 0, 0])
            solve_system!(system, beam)

            # Rotation by 2*theta = 45 degrees, so polarization should be along [1, 0, 1]/sqrt(2)
            E = BMO.polarization(last(beam.rays))
            @test abs(dot(E, [1.0, 0.0, 1.0] / sqrt(2)))≈1.0 atol=1e-6
        end

        @testset "Flat Half-Wave Plate (rectangular)" begin
            # Fast axis oriented at 45 degrees around y-axis.
            hwp = HalfWavePlate(5mm, 5mm; fast_axis_angle = deg2rad(45.0))
            system = System([hwp])
            # Incident polarization along x-axis
            beam = Beam([0, -10mm, 0], [0, 1, 0], 1e-6, [1, 0, 0])
            solve_system!(system, beam)

            # Rotation by 2*theta = 90 degrees, so polarization should be along z-axis [0, 0, 1]
            E = BMO.polarization(last(beam.rays))
            @test abs(dot(E, [0.0, 0.0, 1.0]))≈1.0 atol=1e-6
        end

        @testset "Flat Quarter-Wave Plate" begin
            # Fast axis along x-axis (0 degrees)
            qwp = QuarterWavePlate(5mm; fast_axis_angle = 0.0)
            system = System([qwp])
            # Incident polarization at 45 degrees in x-z plane
            beam = Beam([0, -10mm, 0], [0, 1, 0], 1e-6, [1, 0, 1] / sqrt(2))
            solve_system!(system, beam)

            # Output should be circularly polarized
            E = BMO.polarization(last(beam.rays))
            @test abs(E[1])≈1 / sqrt(2) atol=1e-6
            @test abs(E[3])≈1 / sqrt(2) atol=1e-6
            @test abs(dot(E, [1, 0, im] / sqrt(2)))≈1.0 atol=1e-6
        end

        @testset "Thick/Plate Half-Wave Plate" begin
            # 1mm thick substrate of index n = 1.5, fast axis at 22.5 degrees
            hwp_plate = HalfWavePlate(5mm, 1mm, n -> 1.5; fast_axis_angle = deg2rad(22.5))
            system = System([hwp_plate])

            # Incident polarization along x-axis at normal incidence
            beam = Beam([0, -10mm, 0], [0, 1, 0], 1e-6, [1, 0, 0])
            solve_system!(system, beam)

            # Output polarization should still be rotated by 45 degrees
            E = BMO.polarization(last(beam.rays))
            @test abs(dot(normalize(E), [1.0, 0.0, 1.0] / sqrt(2)))≈1.0 atol=1e-6
            # Substrate thickness is 1mm, so total propagation distance from -10mm should exit at y = 1mm (i.e. total delta 11mm)
            # Since substrate front face is at y = 0, back face is at y = 1mm.
            @test BMO.position(last(beam.rays))[2]≈1mm atol=1e-6
        end

        @testset "Round waveplate uses circular mesh" begin
            hwp = BMO.HalfWaveplate(0.01)
            @test size(BMO.faces(BMO.shape(hwp)), 1) == 30
            qwp = BMO.QuarterWaveplate(0.01)
            @test size(BMO.faces(BMO.shape(qwp)), 1) == 30
        end

        @testset "Half-waveplate rotates linear polarization" begin
            hwp = BMO.HalfWaveplate(0.01, 0.01)
            BMO.yrotate3d!(hwp, deg2rad(30))
            system = StaticSystem([hwp])
            E0 = BMO.electric_field(1)
            ray = PolarizedRay([0, -0.05, 0], [0, 1.0, 0], 1000e-9, [E0, 0, 0])
            beam = Beam(ray)
            solve_system!(system, beam)
            E = BMO.polarization(last(BMO.rays(beam)))
            @test abs2(E[1]) / abs2(E0)≈cosd(60)^2 atol=1e-6
            @test abs2(E[3]) / abs2(E0)≈sind(60)^2 atol=1e-6
        end

        @testset "Wavelength dependent retardance" begin
            λ0 = 1000e-9
            retardance(λ) = π * (λ / λ0)
            hwp = BMO.Waveplate(0.01, 0.01, retardance)
            BMO.yrotate3d!(hwp, deg2rad(30))
            system = StaticSystem([hwp])
            E0 = BMO.electric_field(1)
            ray = PolarizedRay([0, -0.05, 0], [0, 1.0, 0], λ0, [E0, 0, 0])
            beam = Beam(ray)
            solve_system!(system, beam)
            E = BMO.polarization(last(BMO.rays(beam)))
            @test abs2(E[1]) / abs2(E0)≈cosd(60)^2 atol=1e-6
            @test abs2(E[3]) / abs2(E0)≈sind(60)^2 atol=1e-6
        end

        @testset "Arbitrary Jones matrix" begin
            J = BMO.XZBasis(0, 1, 1, 0)
            wp = BMO.Waveplate(0.01, 0.01, J)
            system = StaticSystem([wp])
            E0 = BMO.electric_field(1)
            ray = PolarizedRay([0, -0.05, 0], [0, 1.0, 0], 1000e-9, [E0, 0, 0])
            beam = Beam(ray)
            solve_system!(system, beam)
            E = BMO.polarization(last(BMO.rays(beam)))
            @test abs2(E[1]) < 1e-12
            @test abs2(E[3]) / abs2(E0)≈1 atol=1e-6
        end

        @testset "Quarter-waveplate creates circular polarization" begin
            qwp = BMO.QuarterWaveplate(0.01, 0.01)
            BMO.yrotate3d!(qwp, deg2rad(45))
            system = StaticSystem([qwp])
            E0 = BMO.electric_field(1)
            ray = PolarizedRay([0, -0.05, 0], [0, 1.0, 0], 1000e-9, [E0, 0, 0])
            beam = Beam(ray)
            solve_system!(system, beam)
            E = BMO.polarization(last(BMO.rays(beam)))
            @test abs(E[1])≈abs(E[3]) atol=1e-6
            @test angle(E[3]) - angle(E[1])≈π / 2 atol=1e-6
        end

        @testset "Polarizing beam splitter separates linear polarizations" begin
            pbs = BMO.PolarizingBeamSplitter(0.01, 0.01)
            BMO.translate3d!(pbs, [0, 0.01, 0])
            system = StaticSystem([pbs])
            E0 = BMO.electric_field(1)

            rayx = PolarizedRay([0, -0.05, 0], [0, 1.0, 0], 1000e-9, [E0, 0, 0])
            beamx = Beam(rayx)
            solve_system!(system, beamx)
            tx = BMO.polarization(first(BMO.rays(beamx.children[1])))
            rx = BMO.polarization(first(BMO.rays(beamx.children[2])))
            @test abs2(tx[1]) / abs2(E0)≈1 atol=1e-6
            @test abs2(rx[3]) / abs2(E0)≈0 atol=1e-6

            rayz = PolarizedRay([0, -0.05, 0], [0, 1.0, 0], 1000e-9, [0, 0, E0])
            beamz = Beam(rayz)
            solve_system!(system, beamz)
            tz = BMO.polarization(first(BMO.rays(beamz.children[1])))
            rz = BMO.polarization(first(BMO.rays(beamz.children[2])))
            @test abs2(tz[1]) / abs2(E0)≈0 atol=1e-6
            @test abs2(rz[3]) / abs2(E0)≈1 atol=1e-6
        end

        @testset "HWP before polarizing beamsplitter" begin
            hwp = BMO.HalfWaveplate(0.01, 0.01)
            BMO.yrotate3d!(hwp, deg2rad(22.5))
            pbs = BMO.PolarizingBeamSplitter(0.01, 0.01)
            BMO.translate3d!(pbs, [0, 0.01, 0])
            system = StaticSystem([hwp, pbs])
            E0 = BMO.electric_field(1)
            ray = PolarizedRay([0, -0.05, 0], [0, 1.0, 0], 1000e-9, [E0, 0, 0])
            beam = Beam(ray)
            solve_system!(system, beam)
            tE = BMO.polarization(first(BMO.rays(beam.children[1])))
            rE = BMO.polarization(first(BMO.rays(beam.children[2])))
            It = abs2(tE[1]) / abs2(E0)
            Ir = abs2(rE[3]) / abs2(E0)
            @test It + Ir≈1 atol=1e-6
            @test It > 0 && Ir > 0
        end

        @testset "HWP before polarizing cube beamsplitter" begin
            hwp = BMO.HalfWaveplate(0.01, 0.01)
            BMO.yrotate3d!(hwp, deg2rad(22.5))
            pcbs = BMO.PolarizingCubeBeamsplitter(0.01, n -> 1.0)
            BMO.translate3d!(pcbs, [0, 0.015, 0])
            system = StaticSystem([hwp, pcbs])
            E0 = BMO.electric_field(1)
            ray = PolarizedRay([0, -0.05, 0], [0, 1.0, 0], 1000e-9, [E0, 0, 0])
            beam = Beam(ray)
            solve_system!(system, beam)
            @test !isempty(BMO.rays(beam.children[1]))
            @test !isempty(BMO.rays(beam.children[2]))
        end

        @testset "HWP before tilted polarizing cube beamsplitter" begin
            hwp = BMO.HalfWaveplate(0.01, 0.01)
            BMO.yrotate3d!(hwp, deg2rad(22.5))
            pcbs = BMO.PolarizingCubeBeamsplitter(0.01, n -> 1.0)
            BMO.translate3d!(pcbs, [0, 0.015, 0])
            BMO.xrotate3d!(pcbs, deg2rad(10))
            system = StaticSystem([hwp, pcbs])
            E0 = BMO.electric_field(1)
            ray = PolarizedRay([0, -0.05, 0], [0, 1.0, 0], 1000e-9, [E0, 0, 0])
            beam = Beam(ray)
            solve_system!(system, beam)
            @test !isempty(BMO.rays(beam.children[1]))
            @test !isempty(BMO.rays(beam.children[2]))
        end

        @testset "Polarizing plate beamsplitter splits 45° polarization" begin
            ppbs = BMO.RectangularPolarizingPlateBeamsplitter(0.01, 0.01, 0.005, n -> 1.0)
            BMO.translate3d!(ppbs, [0, 0.01, 0])
            system = StaticSystem([ppbs])
            E0 = BMO.electric_field(1)
            E = E0 / sqrt(2)
            ray = PolarizedRay([0, -0.05, 0], [0, 1.0, 0], 1000e-9, [E, 0, E])
            beam = Beam(ray)
            solve_system!(system, beam)
            tE = BMO.polarization(last(BMO.rays(beam.children[1])))
            rE = BMO.polarization(last(BMO.rays(beam.children[2])))
            It = abs2(tE[1]) / abs2(E0)
            Ir = abs2(rE[3]) / abs2(E0)
            @test It≈0.5 atol=1e-6
            @test Ir≈0.5 atol=1e-6
        end

        @testset "QWP mirror roundtrip" begin
            qwp = BMO.QuarterWaveplate(0.01, 0.01)
            BMO.yrotate3d!(qwp, deg2rad(45))
            mir = BMO.Mirror(BMO.RectangularFlatMesh(0.02, 0.02))
            BMO.translate3d!(mir, [0, 0.02, 0])
            system = StaticSystem([qwp, mir])
            E0 = BMO.electric_field(1)
            ray = PolarizedRay([0, -0.05, 0], [0, 1.0, 0], 1000e-9, [E0, 0, 0])
            beam = Beam(ray)
            solve_system!(system, beam)
            finalE = BMO.polarization(last(BMO.rays(beam)))
            Ix = abs2(finalE[1]) / abs2(E0)
            Iz = abs2(finalE[3]) / abs2(E0)
            @test Ix≈0 atol=1e-6
            @test Iz≈1 atol=1e-6
        end

        @testset "Roundtrip with HWP, PBSC, QWP, mirror" begin
            hwp = BMO.HalfWaveplate(0.01, 0.01)
            BMO.yrotate3d!(hwp, deg2rad(45))
            pcbs = BMO.PolarizingCubeBeamsplitter(0.01, n -> 1.0)
            BMO.translate3d!(pcbs, [0, 0.015, 0])
            qwp = BMO.QuarterWaveplate(0.01, 0.01)
            BMO.yrotate3d!(qwp, deg2rad(45))
            BMO.translate3d!(qwp, [0, 0.03, 0])
            mir = BMO.Mirror(BMO.RectangularFlatMesh(0.02, 0.02))
            BMO.translate3d!(mir, [0, 0.05, 0])
            system = StaticSystem([hwp, pcbs, qwp, mir])
            E0 = BMO.electric_field(1)
            ray = PolarizedRay([0, -0.05, 0], [0, 1.0, 0], 1000e-9, [0, 0, E0])
            beam = Beam(ray)
            solve_system!(system, beam)
            rE0 = BMO.polarization(first(BMO.rays(beam.children[2])))
            @test abs2(rE0[3]) / abs2(E0)≈0 atol=1e-6
            rE = BMO.polarization(first(BMO.rays(beam.children[1].children[2])))
            Ir = abs2(rE[3]) / abs2(E0)
            @test Ir≈1 atol=1e-6
        end

        @testset "AstigmaticGaussianBeamlet with Polarizers" begin
            # PolarizationFilter
            filter = PolarizationFilter(20mm)
            system1 = System([filter])
            agb1 = AstigmaticGaussianBeamlet(
                [0.0, -50mm, 0.0], [0.0, 1.0, 0.0], 1000e-9, 1mm;
                E0 = [1.0, 0.0, 0.0], support = [0.0, 0.0, 1.0])
            solve_system!(system1, agb1)
            # Should pass through with polarization preserved
            @test BMO.polarization(last(BMO.rays(agb1.c))) ≈ [1.0, 0.0, 0.0]

            # Waveplate
            hwp = HalfWaveplate(20mm, 20mm)
            BMO.yrotate3d!(hwp, deg2rad(30))
            system2 = System([hwp])
            agb2 = AstigmaticGaussianBeamlet(
                [0.0, -50mm, 0.0], [0.0, 1.0, 0.0], 1000e-9, 1mm;
                E0 = [1.0, 0.0, 0.0], support = [0.0, 0.0, 1.0])
            solve_system!(system2, agb2)
            # Final polarization vector rotated to [cos(60°), 0, -sin(60°)] (with -im global phase)
            final_pol = BMO.polarization(last(BMO.rays(agb2.c)))
            expected = [cosd(60), 0.0, -sind(60)]
            @test final_pol ≈ -im * expected atol=1e-6
        end
    end

    @testset "Power Cutoff in solve_system!" begin
        bs = ThinBeamsplitter(5mm, 5mm; reflectance = 0.5)
        m = SquarePlanoMirror2D(5mm)
        translate3d!(m, [0.0, 20mm, 0.0])
        zrotate3d!(m, deg2rad(10.0))
        system = System([bs, m])

        # Test 1: No power cutoff (default 0.0) -> traces all children
        beam_no_cut = Beam([0.0, -10mm, 0.0], [0.0, 1.0, 0.0], 1e-6, [1.0, 0.0, 0.0])
        solve_system!(system, beam_no_cut)
        @test length(beam_no_cut.children) == 2
        # Transmitted child propagates to mirror, so it has 2 rays
        @test length(BMO.rays(beam_no_cut.children[1])) == 2
        # Reflected child has no obstacles, so it has 1 ray
        @test length(BMO.rays(beam_no_cut.children[2])) == 1

        # Test 2: power_cutoff = 0.6 -> cuts off children since child power (0.5) < 0.6
        beam_cut = Beam([0.0, -10mm, 0.0], [0.0, 1.0, 0.0], 1e-6, [1.0, 0.0, 0.0])
        solve_system!(system, beam_cut; power_cutoff = 0.6)
        @test length(beam_cut.children) == 2
        # Child beams were cut off, so they have only 1 ray (the start/interaction ray)
        @test length(BMO.rays(beam_cut.children[1])) == 1
        @test length(BMO.rays(beam_cut.children[2])) == 1

        # Test 3: Standard unpolarized Ray power cutoff
        beam_unpol = Beam([0.0, -10mm, 0.0], [0.0, 1.0, 0.0], 1e-6)
        solve_system!(system, beam_unpol; power_cutoff = 0.6)
        @test length(beam_unpol.children) == 2
        # Cutoff should also happen here because reflectance=0.5 < 0.6
        @test length(BMO.rays(beam_unpol.children[1])) == 1
        @test length(BMO.rays(beam_unpol.children[2])) == 1
    end
end

end # MODULE
