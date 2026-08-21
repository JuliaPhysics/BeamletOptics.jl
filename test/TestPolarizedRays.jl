module TestPolarizedRays

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics

@testset "Polarized rays" begin
    @testset "Polarization transforms" begin
        # Reflection matrix and lin. x-pol
        J = BMO.SPBasis(-1, 0, 0, 1)
        E0 = [1, 0, 0]

        @testset "90° reflection" begin
            in_dir = [0, 0, 1]
            out_dir = [1, 0, 0]
            nml = normalize([1, 0, -1])
            P90 = BMO._calculate_global_E0(in_dir, out_dir, nml, J)
            @test P90 * E0 ≈ [0, 0, -1]
            @test P90 * in_dir ≈ out_dir
        end

        @testset "0° reflection" begin
            in_dir = [0, 0, 1]
            out_dir = [0, 0, -1]
            nml = normalize([0, 0, -1])
            P00 = BMO._calculate_global_E0(in_dir, out_dir, nml, J)
            @test P00 * E0 ≈ [-1, 0, 0]
            @test P00 * in_dir ≈ out_dir
        end
    end

    @testset "Test error messages" begin
        # Test dir. error msg
        @test_throws ErrorException PolarizedRay(zeros(3), zeros(3))
        @test_throws ErrorException PolarizedRay(zeros(3), ones(3)*eps())
        # Test polarization orthogonal error msg
        @test_throws ErrorException PolarizedRay(zeros(3), [0,1,0], 1e-6, [0,1,1])
    end

    @testset "Mirror reflections" begin
        # Setup system as in https://opg.optica.org/ao/fulltext.cfm?uri=ao-50-18-2855&id=218813
        m1 = SquarePlanoMirror2D(1.0)
        m2 = SquarePlanoMirror2D(1.0)
        m3 = SquarePlanoMirror2D(1.0)
        translate3d!(m2, [2, 0, 0])
        translate3d!(m3, [2, 2, 0])
        zrotate3d!(m1, deg2rad(-90))
        yrotate3d!(m1, deg2rad(45))
        zrotate3d!(m2, deg2rad(45))
        xrotate3d!(m3, deg2rad(135))

        system = StaticSystem([m1, m2, m3])

        I0_1 = 1
        I0_2 = 5
        lin_x_pol = [I0_1, 0, 0]
        lin_y_pol = [0, I0_2, 0]

        # Beam of polarized rays
        ray = PolarizedRay([0.0, 0, -2], [0, 0, 1], 1000e-9, lin_x_pol)
        beam = Beam(ray)

        @testset "x-Polarization" begin
            beam = Beam([0.0, 0, -2], [0, 0, 1], 1000e-9, lin_x_pol)
            # test tracing
            solve_system!(system, beam)
            @test BMO.polarization(beam.rays[1]) ≈ lin_x_pol
            @test BMO.polarization(beam.rays[2]) ≈ [0, 0, -I0_1]
            @test BMO.polarization(beam.rays[3]) ≈ [0, 0, I0_1]
            @test BMO.polarization(beam.rays[4]) ≈ [0, -I0_1, 0]
            @test length(beam) ≈ 6.0
        end

        @testset "y-Polarization" begin
            beam = Beam([0.0, 0, -2], [0, 0, 1], 1000e-9, lin_y_pol)
            translate3d!(m3, [0, 2, 0])
            # test re-solving after a kinematic change
            solve_system!(system, beam)
            @test BMO.polarization(beam.rays[1]) ≈ lin_y_pol
            @test BMO.polarization(beam.rays[2]) ≈ [0, -I0_2, 0]
            @test BMO.polarization(beam.rays[3]) ≈ [I0_2, 0, 0]
            @test BMO.polarization(beam.rays[4]) ≈ [-I0_2, 0, 0]
            @test length(beam) ≈ 8.0
        end
    end


    @testset "Brewster windows" begin
        brewster_angle(n) = atan(n)
        # Testcase based on 5 successive Brewster windows
        n = 1.5
        θb = brewster_angle(n)
        d = 0.1
        # Calculate transmission efficiency
        rs, rp, ts, tp = BMO.fresnel_coefficients(θb, n)
        Ts = 1 - abs2(rs)
        Tp = 1 - abs2(rp)
        # Setup testcase
        s1 = BMO.CuboidMesh(1.0, d, 1.0)
        s2 = BMO.CuboidMesh(1.0, d, 1.0)
        s3 = BMO.CuboidMesh(1.0, d, 1.0)
        s4 = BMO.CuboidMesh(1.0, d, 1.0)
        s5 = BMO.CuboidMesh(1.0, d, 1.0)
        l1 = Lens(s1, x -> n)
        l2 = Lens(s2, x -> n)
        l3 = Lens(s3, x -> n)
        l4 = Lens(s4, x -> n)
        l5 = Lens(s5, x -> n)
        translate3d!.([l1, l2, l3, l4, l5], Ref([-0.5, -d / 2, -0.5]))
        BMO.set_new_origin3d!.(BMO.shape.([l1, l2, l3, l4, l5]))
        translate3d!(l2, [0, 0.5, -1d / 2])
        translate3d!(l3, [0, 1.0, -2d / 2])
        translate3d!(l4, [0, 1.5, -3d / 2])
        translate3d!(l5, [0, 2.0, -4d / 2])
        xrotate3d!.([l1, l2, l3, l4, l5], -θb)
        # Solve system of s- and p-polarized beams
        system = StaticSystem([l1, l2, l3, l4, l5])
        x_pol_ray = PolarizedRay(
            [-0.1, -1, 0], [0, 1.0, 0], 1000e-9, [BMO.electric_field(1), 0, 0])
        z_pol_ray = PolarizedRay(
            [+0.1, -1, 0], [0, 1.0, 0], 1000e-9, [0, 0, BMO.electric_field(1)])
        s_beam = Beam(x_pol_ray)
        p_beam = Beam(z_pol_ray)
        solve_system!(system, s_beam)
        solve_system!(system, p_beam)
        # Since system is non-focussing, calculate pseudo-intensity
        pseudo_Is = abs2(BMO.polarization(last(BMO.rays(s_beam)))[1]) /
                    (2 * BMO.Z_vacuum)
        pseudo_Ip = abs2(BMO.polarization(last(BMO.rays(p_beam)))[3]) /
                    (2 * BMO.Z_vacuum)
        # Test against m interfaces
        m = length(system.objects) * 2
        @test pseudo_Is ≈ Ts^m
        @test pseudo_Ip ≈ Tp^m
    end

    @testset "Fresnel rhomb" begin
        # Create Fresnel rhomb with n=1.5 and θ=53.3° for quarter-wave plate effect
        n = 1.5
        s1 = BMO.CuboidMesh(0.5, 1.25, 0.5, deg2rad(53.3))
        l1 = Lens(s1, x -> n)
        translate3d!(l1, [-0.25, 0, -0.25])
        BMO.set_new_origin3d!(s1)
        # Rotate prism to obtain 45° beam input polarization
        yrotate3d!(l1, deg2rad(135))
        # Solve system
        system = StaticSystem([l1])
        ray = PolarizedRay(
            [0, -1, 0], [0, 1.0, 0], 1000e-9, [0, 0, BMO.electric_field(1)])
        beam = Beam(ray)
        solve_system!(system, beam)
        # Assumes propagation along the y-axis after rhomb, calculate polarization state
        Ex = getindex.(BMO.polarization.(beam.rays), 1)
        Ey = getindex.(BMO.polarization.(beam.rays), 2)
        Ez = getindex.(BMO.polarization.(beam.rays), 3)
        # Test for circular polarization and Ey error
        phi = angle(last(Ez)) - angle(last(Ex))
        @test phi ≈ π / 2
        @test abs(last(Ey)) < 2e-14
    end
end

end # MODULE