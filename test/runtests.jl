using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics
using AbstractTrees

const BMO = BeamletOptics

const mm = 1e-3
const inch = 25.4mm

include(joinpath(@__DIR__, "TestUtils.jl"))
include(joinpath(@__DIR__, "TestAbstractTypes.jl"))
include(joinpath(@__DIR__, "TestRays.jl"))
include(joinpath(@__DIR__, "TestBeams.jl"))
include(joinpath(@__DIR__, "TestBeamGroups.jl"))
include(joinpath(@__DIR__, "Geometry", "TestMesh.jl"))
include(joinpath(@__DIR__, "Geometry", "TestSDFs.jl"))
include(joinpath(@__DIR__, "TestSystem.jl"))
include(joinpath(@__DIR__, "TestObjectGroups.jl"))
include(joinpath(@__DIR__, "Lenses", "TestSphericalLenses.jl"))
include(joinpath(@__DIR__, "Lenses", "TestSurfaces.jl"))
include(joinpath(@__DIR__, "Lenses", "TestAsphericalLenses.jl"))
include(joinpath(@__DIR__, "Lenses", "TestCylindricalLenses.jl"))
include(joinpath(@__DIR__, "Lenses", "TestDoubleGaussLens.jl"))
include(joinpath(@__DIR__, "TestGaussianBeamlet.jl"))

@testset "Detectors" begin
    @testset "Testing Spotdetector" begin
        # Set up tilted spot detection screens
        α = 45
        sd = Detector(1.0)
        system = System([sd])
        translate3d!(sd, [0, 1, 0])
        zrotate3d!(sd, deg2rad(α))
        beam = Beam([0, 0, 0], [0, 1, 0], 1e-6)
        # Trace beams in x-y-plane
        xs = LinRange(-0.25, 0.25, 10)
        pts = Vector{Point3{Float64}}(undef, length(xs))
        for (i, x) in enumerate(xs)
            BMO.position!(first(beam.rays), Point3{Float64}(x, 0, 0))
            solve_system!(system, beam)
            # compare ray intersection to stored data
            hit = last(BMO.hits(sd))
            pos_dta = BMO.hit_point(hit)
            ray = last(beam.rays)
            pos_ray = position(ray) +
                      length(BMO.intersection(ray)) * BMO.direction(ray)
            @test pos_ray ≈ pos_dta
        end
        # Test spot diagram projection
        pts = spot_diagram(sd)
        # Test reset function
        empty!(sd)
        @test isnothing(BMO.hits(sd))
    end
end

@testset "Interference" begin
    @testset "Pre-Beamsplitter tests with seperate beams" begin
        # Gauss beam parameters (selected for ring fringes)
        w0 = 0.01e-3
        λ = 1000e-9
        M2 = 1
        P0 = 1e-3
        I0 = 2 * P0 / (π * w0^2)
        E0 = BMO.electric_field(I0)
        zR = BMO.rayleigh_range(λ, w0, M2)
        # Detector parameters
        z = 0.1     # distance to detector
        l = 1e-2    # detector size
        n = 1000    # detector grid resolution
        # Lens parameters
        R1 = R2 = d = 0.01
        nl = 1.5
        f = BMO.lensmakers_eq(R1, -R2, nl)
        # Raytracing system (for all tests)
        pd_l = Detector(l) # n
        pd_s = Detector(l / 10) # n ÷ 10
        ln = ThinLens(R1, R2, d, nl)
        translate3d!(pd_l, [0, z, 0])
        translate3d!(pd_s, [0, z, 0])
        translate3d!(ln, [0, z - f - thickness(ln) / 2, 0])

        @testset "Testing fringe pattern" begin
            system = System(pd_l)
            Δz = 5e-3   # arm length difference
            # Analytic solution
            xs = ys = LinRange(-l / 2, l / 2, n)
            screen = zeros(ComplexF64, length(xs), length(ys))
            for (j, y) in enumerate(ys)
                for (i, x) in enumerate(xs)
                    r = sqrt(x^2 + y^2)
                    screen[i, j] += BMO.electric_field(r, z, E0, w0, λ, M2)
                    screen[i, j] += BMO.electric_field(r, z + Δz, E0, w0, λ, M2)
                end
            end
            # Numerical solution
            empty!(pd_l)
            g_1 = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0],
                λ,
                w0,
                M2 = M2,
                P0 = P0)
            g_2 = GaussianBeamlet([0.0, -Δz, 0], [0.0, 1, 0],
                λ,
                w0,
                M2 = M2,
                P0 = P0)
            solve_system!(system, g_1)
            solve_system!(system, g_2)

            # Compare solutions
            I_analytical = intensity.(screen)
            ~, ~, I_numerical = intensity(pd_l; n, x_min=-l/2, x_max=l/2, z_min=-l/2, z_max=l/2)
            Pt = optical_power(pd_l; n, x_min=-l/2, x_max=l/2, z_min=-l/2, z_max=l/2)
            @test all(isapprox.(I_analytical, I_numerical, atol = 2e-1))
            @test isapprox(Pt, 2 * P0, atol = 3e-5)
        end

        @testset "Testing λ phase shift" begin
            system = System([pd_s, ln])
            # Numerical solution
            Δz = LinRange(0, λ, 50)
            Pt_numerical = zeros(length(Δz))
            for (i, z_i) in enumerate(Δz)
                empty!(pd_s)
                g_1 = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0],
                    λ,
                    w0,
                    M2 = M2,
                    P0 = P0)
                g_2 = GaussianBeamlet([0.0, z_i, 0], [0.0, 1, 0],
                    λ,
                    w0,
                    M2 = M2,
                    P0 = P0)
                solve_system!(system, g_1)
                solve_system!(system, g_2)
                Pt_numerical[i] = BMO.optical_power(pd_s)
                # Test length/opl function
                @test length(g_1) ≈ z
                @test length(g_1) ≈ BMO.optical_path_length(g_1) - thickness(ln) * (nl - 1)
                @test length(g_2) ≈ length(g_1) - z_i
            end
            # Analytical solution (cosine over Δz), ref. power is 4*P0 since beamsplitter is missing
            Pt_analytical = 4 * P0 * [(cos(2π * z / (maximum(Δz))) + 1) / 2 for z in Δz]

            # Compare detectors (this also tests correct behavior when focussing the beam)
            @test all(isapprox.(Pt_numerical, Pt_analytical, atol = 1e-4))
        end
    end

    @testset "Michelson Interferometer" begin
        # setup Michelson Interferometer
        l_0 = 0.1
        pd_size = BMO.inch / 5
        pd_resolution = 100
        m1 = SquarePlanoMirror2D(BMO.inch)
        m2 = SquarePlanoMirror2D(BMO.inch)
        bs = ThinBeamsplitter(BMO.inch, reflectance = 0.5)
        pd = Detector(pd_size)
        translate3d!(m1, [l_0, 0, 0])
        translate3d!(m2, [0, l_0, 0])
        translate3d!(pd, [-l_0, 0, 0])
        zrotate3d!(bs, deg2rad(45))
        zrotate3d!(m1, deg2rad(90))
        zrotate3d!(pd, deg2rad(90))

        system = System([m1, m2, bs, pd])

        # Test correct values for reflectivity/transmission
        @test isvalid(bs)

        @testset "Equal armlength MI - integrated power" begin
            # setup 635 nm laser with 0.1 mm waist for fast divergence
            λ = 635e-9
            P_0 = 5e-3
            beam = GaussianBeamlet([0, -l_0, 0], [0, 1.0, 0], λ, 1e-4, P0 = P_0)

            # Shift mirror #2 by -λ to +λ
            lambdas = LinRange(-λ, λ, 200)

            path_length_numerical = zeros(length(lambdas))
            optical_pwr_numerical = zeros(length(lambdas))

            for (i, lambda) in enumerate(lambdas)
                translate_to3d!(m2, [0, l_0, 0] + [0, lambda, 0])
                empty!(pd)
                solve_system!(system, beam)

                # Moving mirror path length
                path_length_numerical[i] = length(beam.children[1].children[2])
                optical_pwr_numerical[i] = BMO.optical_power(
                    pd;
                    # restore pre-v0.11 behavior, e.g. no autolims
                    n=pd_resolution,
                    x_min=-pd_size/2,
                    x_max=pd_size/2,
                    z_min=-pd_size/2,
                    z_max=pd_size/2
                )
            end

            path_length_analytical = @. 2 * lambdas + 4l_0
            optical_pwr_analytical = @. P_0 * (1 / 2 * cos(2π * (2lambdas / λ) + π) + 1 / 2)

            # Compare correct PD signal and λ shift in moving arm
            @test all(isapprox.(optical_pwr_analytical, optical_pwr_numerical, atol = 5e-6))
            @test all(isapprox.(path_length_analytical, path_length_numerical))
        end

        @testset "Unequal armlength MI - electrical field" begin
            λ = 635e-9
            w0 = 1e-4
            P0 = 1e-3
            M2 = 1
            I0 = 2 * P0 / (π * w0^2)
            E0 = BMO.electric_field(I0) * 1 / sqrt(2)^2
            zR = BMO.rayleigh_range(λ, w0, M2)

            beam = GaussianBeamlet([0, -l_0, 0], [0, 1.0, 0], λ, w0, P0 = P0, M2 = M2)

            # arm length diff
            Δl = 1 * l_0
            translate_to3d!(m2, [0, l_0 + Δl, 0])

            # numerical solution
            empty!(pd)
            solve_system!(system, beam)

            xs, ys, Efield = electric_field(pd; n=pd_resolution)

            # analytical solution
            short_arm = 4l_0
            long_arm = short_arm + 2Δl
            # xs = ys = LinRange(-pd_size / 2, pd_size / 2, pd_resolution)
            screen = zeros(ComplexF64, length(xs), length(ys))
            for (j, y) in enumerate(ys)
                for (i, x) in enumerate(xs)
                    r = sqrt(x^2 + y^2)
                    screen[i, j] += BMO.electric_field(
                        r, short_arm, E0, w0, λ, M2)
                    screen[i, j] += BMO.electric_field(
                        r, long_arm, E0, w0, λ, M2) * exp(im * pi)
                end
            end

            Re_analytical = real.(screen)
            Im_analytical = imag.(screen)

            Re_numerical = real(Efield)
            Im_numerical = imag(Efield)

            # Compare solutions, units V/m
            @test all(isapprox.(Re_analytical, Re_numerical, atol = 5e-2))
            @test all(isapprox.(Im_analytical, Im_numerical, atol = 5e-2))
        end
    end

    @testset "Testing power conservation" begin
        # variables
        P0 = 0.5 # W
        l0 = 0.1 # m
        w0 = 0.5e-3
        λ = 1064e-9

        bs = ThinBeamsplitter(10e-3)
        pd_resolution = 100
        pd_1 = Detector(10e-3)
        pd_2 = Detector(10e-3)

        zrotate3d!(bs, deg2rad(45))
        translate3d!(pd_1, [0, l0, 0])
        zrotate3d!(pd_1, deg2rad(180))

        translate3d!(pd_2, [l0, 0, 0])
        zrotate3d!(pd_2, deg2rad(90))

        # add BS and PD orientation error
        zrotate3d!(bs, deg2rad(0.017))
        zrotate3d!(pd_1, deg2rad(10))
        xrotate3d!(pd_1, deg2rad(15))

        # define system and beams -> solve
        system = System([bs, pd_1, pd_2])

        phis = LinRange(0, 2pi, 25)
        p1 = similar(phis)
        p2 = similar(phis)

        l1 = GaussianBeamlet([0, -l0, 0], [0, 1.0, 0], λ, w0; P0)
        l2 = GaussianBeamlet([-l0, 0, 0], [1.0, 0, 0], λ, w0; P0)

        E0_buffer = l1.E0

        for (i, phi) in enumerate(phis)
            # Iterate over relative phase shifts, use retracing
            l1.E0 = E0_buffer * exp(im * phi)
            empty!(pd_1)
            empty!(pd_2)
            solve_system!(system, l1)
            solve_system!(system, l2)
            p1[i] = BMO.optical_power(pd_1)
            p2[i] = BMO.optical_power(pd_2)
            # Test power conservation
            @test p1[i] + p2[i] - 2P0 < 1e-4 # W
        end
    end
end

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
            BMO.polarization!(ray, lin_x_pol)
            # test tracing
            solve_system!(system, beam)
            @test BMO.polarization(beam.rays[1]) ≈ lin_x_pol
            @test BMO.polarization(beam.rays[2]) ≈ [0, 0, -I0_1]
            @test BMO.polarization(beam.rays[3]) ≈ [0, 0, I0_1]
            @test BMO.polarization(beam.rays[4]) ≈ [0, -I0_1, 0]
            @test length(beam) == 6.0
        end

        @testset "y-Polarization" begin
            BMO.polarization!(ray, lin_y_pol)
            translate3d!(m3, [0, 2, 0])
            # test retracing
            solve_system!(system, beam)
            @test BMO.polarization(beam.rays[1]) ≈ lin_y_pol
            @test BMO.polarization(beam.rays[2]) ≈ [0, -I0_2, 0]
            @test BMO.polarization(beam.rays[3]) ≈ [I0_2, 0, 0]
            @test BMO.polarization(beam.rays[4]) ≈ [-I0_2, 0, 0]
            @test length(beam) == 8.0
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

    @testset "Mach-Zehnder Interferometer" begin
        # setup MZI
        m1 = SquarePlanoMirror2D(BMO.inch)
        m2 = SquarePlanoMirror2D(BMO.inch)
        b1 = ThinBeamsplitter(BMO.inch, reflectance = 0.5)
        b2 = ThinBeamsplitter(BMO.inch, reflectance = 0.5)

        system = StaticSystem([m1, m2, b1, b2])

        translate3d!(b1, [0 * BMO.inch, 0 * BMO.inch, 0])
        translate3d!(b2, [2 * BMO.inch, 2 * BMO.inch, 0])
        translate3d!(m1, [0 * BMO.inch, 2 * BMO.inch, 0])
        translate3d!(m2, [2 * BMO.inch, 0 * BMO.inch, 0])

        # Rotate with consideration to mirror/bs normal
        zrotate3d!(b1, deg2rad(360 - 135))
        zrotate3d!(b2, deg2rad(45))
        zrotate3d!(m1, deg2rad(360 - 135))
        zrotate3d!(m2, deg2rad(45))

        ray = PolarizedRay([0, -0.1, 0], [0.0, 1.0, 0], 1000e-9, [0, 0, 1])
        beam = Beam(ray)

        @testset "z-polarized ray along y-axis" begin
            # Solve with z-polarized ray along y-axis
            BMO.polarization!(ray, [0, 0, 1])
            solve_system!(system, beam)

            # Extract E0s: t - transmitted, r - reflected
            t = BMO.polarization(beam.children[1].rays[1])
            r = BMO.polarization(beam.children[2].rays[1])
            tr = BMO.polarization(beam.children[1].rays[2])
            rr = BMO.polarization(beam.children[2].rays[2])
            trt = BMO.polarization(beam.children[1].children[1].rays[1])
            trr = BMO.polarization(beam.children[1].children[2].rays[1])
            rrt = BMO.polarization(beam.children[2].children[1].rays[1])
            rrr = BMO.polarization(beam.children[2].children[2].rays[1])

            # Test phase flips
            @test t[3] ≈ sqrt(2) / 2
            @test r[3] ≈ -sqrt(2) / 2
            @test tr[3] ≈ -sqrt(2) / 2
            @test rr[3] ≈ sqrt(2) / 2
            @test trt ≈ rrr
            @test trr ≈ rrt
        end

        @testset "x-polarized ray along y-axis" begin
            # Test num. of leaves before retracing
            @test length(collect(Leaves(beam))) == 4
            # Retrace with z-polarized ray along y-axis
            BMO.polarization!(ray, [1, 0, 0])
            solve_system!(system, beam)

            # Extract E0s: t - transmitted, r - reflected
            t = BMO.polarization(beam.children[1].rays[1])
            r = BMO.polarization(beam.children[2].rays[1])
            tr = BMO.polarization(beam.children[1].rays[2])
            rr = BMO.polarization(beam.children[2].rays[2])
            trt = BMO.polarization(beam.children[1].children[1].rays[1])
            trr = BMO.polarization(beam.children[1].children[2].rays[1])
            rrt = BMO.polarization(beam.children[2].children[1].rays[1])
            rrr = BMO.polarization(beam.children[2].children[2].rays[1])

            # Test phase flips
            @test t[1] ≈ sqrt(2) / 2
            @test r[2] ≈ -sqrt(2) / 2
            @test tr ≈ r
            @test rr ≈ t
            @test trt ≈ rrr
            @test trr ≈ rrt
        end
    end
end

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
        @test J*2 == R*2
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
            ray_dir = Rfil[:,2]
            ray_pos = position(filter) - 10mm * ray_dir
            local_x = Rfil[:,1]
            local_z = Rfil[:,3]
            pol_vec = local_x
            lambda = 1e-6
            beam = Beam(ray_pos, ray_dir, lambda, pol_vec)            
            # step ref. vector by angle increment 
            RotMat = BMO.rotate3d(ray_dir, deg2rad(step(thetas)))
            for i in eachindex(thetas)
                solve_system!(system, beam)
                E1 = BMO.polarization(BMO.rays(beam)[2])
                i_n[i] = pseudo_I(E1)
                i_a[i] = malus_law(thetas[i]-1)
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
        translate3d!(filter, [0,10mm,0])
        xrotate3d!(filter, deg2rad(45))
        zrotate3d!(filter, deg2rad(30))

        @testset "Pol. filter tilted incidence - rotation around ray optical axis" begin
            Rfil = orientation(filter)
            ray_dir = Rfil[:,2]
            ray_pos = position(filter) - 10mm * ray_dir
            local_x = Rfil[:,1]
            local_z = Rfil[:,3]
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
end

@testset "Beamsplitters" begin
    N0 = 1.5
    @testset "Testing RectangularPlateBeamsplitter with Beam" begin
        # Init splitter
        N0 = 1.5
        pbs = RectangularPlateBeamsplitter(36mm, 25mm, 1mm, n -> N0)
        system = System([pbs])
        beam = Beam([0, -50mm, 0], [0, 1, 0], 1e-6)
        # Trace normally
        zrotate3d!(pbs, deg2rad(45))
        solve_system!(system, beam)

        @testset "Test pos/dir" begin
            @test position(pbs) == zeros(3)
            @test orientation(pbs) ≈ orientation(pbs.substrate)
        end

        @testset "Test children after tracing" begin
            p = beam.rays
            t = beam.children[1].rays
            r = beam.children[2].rays
            # no of rays
            @test length(p) == 1
            @test length(t) == 2
            @test length(r) == 1
            # correct ref. index
            @test all(BMO.refractive_index.(p) .== 1)
            @test all(BMO.refractive_index.(t) .== [N0, 1])
            @test all(BMO.refractive_index.(r) .== 1)
            # correct dir
            @test BMO.direction(first(p)) ≈ BMO.direction(last(t))
            @test BMO.direction(first(r)) ≈ [1, 0, 0]
        end

        # Retrace backside
        zrotate3d!(pbs, π)
        solve_system!(system, beam)

        @testset "Test children after retracing" begin
            p = beam.rays
            t = beam.children[1].rays
            r = beam.children[2].rays
            # no of rays
            @test length(p) == 2
            @test length(t) == 1
            @test length(r) == 2
            # correct ref. index
            @test all(BMO.refractive_index.(p) .== [1, N0])
            @test all(BMO.refractive_index.(t) .== 1)
            @test all(BMO.refractive_index.(r) .== [N0, 1])
            # correct dir
            @test BMO.direction(first(p)) ≈ BMO.direction(last(t))
            @test BMO.direction(last(r)) ≈ [1, 0, 0]
        end
    end

    @testset "Testing CubeBeamsplitter with Beam" begin
        # Init splitter
        cbs = CubeBeamsplitter(25e-3, n -> N0)
        translate3d!(cbs, [0, 50mm, 0])
        system = System([cbs])
        beam = Beam([0, 0, 0], [0, 1, 0], 1e-6)

        @testset "Initial CBS tracing" begin
            # Trace normally
            solve_system!(system, beam)
            # Test correct ray length, ref. indices, dirs
            p = BMO.rays(beam)
            t = BMO.rays(beam.children[1])
            r = BMO.rays(beam.children[2])

            @test length(p) == 2
            @test length(r) == 2
            @test length(t) == 2
            @test BMO.refractive_index.(p) == [1, N0]
            @test BMO.refractive_index.(t) == [N0, 1]
            @test BMO.refractive_index.(r) == [N0, 1]
            @test BMO.direction(last(t)) ≈ BMO.direction(first(p))
            @test BMO.direction(last(r)) ≈ [-1, 0, 0]
        end

        @testset "Retrace after 45° CBS rotation" begin
            # Retrace
            zrotate3d!(cbs, π / 2)
            solve_system!(system, beam)

            # Test correct ray dirs
            p = BMO.rays(beam)
            t = BMO.rays(beam.children[1])
            r = BMO.rays(beam.children[2])

            @test BMO.direction(last(t)) == BMO.direction(first(p))
            @test BMO.direction(last(t)) == [0, 1, 0]
        end

        @testset "Retrace CBS backside" begin
            # Retrace backside
            zrotate3d!(cbs, π / 2)
            solve_system!(system, beam)

            # Test correct ray length, ref. indices, dirs
            p = BMO.rays(beam)
            t = BMO.rays(beam.children[1])
            r = BMO.rays(beam.children[2])

            @test length(p) == 2
            @test length(r) == 2
            @test length(t) == 2
            @test BMO.refractive_index.(p) == [1, N0]
            @test BMO.refractive_index.(t) == [N0, 1]
            @test BMO.refractive_index.(r) == [N0, 1]
            @test BMO.direction(last(t)) ≈ BMO.direction(first(p))
            @test BMO.direction(last(r)) ≈ [-1, 0, 0]
        end
    end
end

@testset "Dummy objects" begin
    # Setup dummy cube and test beam
    cube_shape = BMO.CubeMesh(1)
    translate3d!(cube_shape, -[0.5, 0, 0.5])
    translate3d!(cube_shape, [0, 5, 0])
    @testset "IntersectableObject" begin
        beam = Beam([0, 0, 0], [0, 1, 0], 1e-6)
        intersectable = IntersectableObject(cube_shape)
        system = System([intersectable])
        solve_system!(system, beam)
        # Test nothing interaction
        @test length(BMO.rays(beam)) == 1
        @test BMO.object(BMO.intersection(last(BMO.rays(beam)))) == intersectable
        @test BMO.shape(BMO.intersection(last(BMO.rays(beam)))) == cube_shape
        @test isnothing(BMO.interact3d(system, intersectable, beam, first(BMO.rays(beam))))
    end

    @testset "NonInteractableObject" begin
        beam = Beam([0, 0, 0], [0, 1, 0], 1e-6)
        noninteract = NonInteractableObject(cube_shape)
        system = System([noninteract])
        solve_system!(system, beam)
        # Test nothing interaction and intersection
        @test length(BMO.rays(beam)) == 1
        @test isnothing(BMO.intersection(last(BMO.rays(beam))))
        @test isnothing(BMO.interact3d(system, noninteract, beam, first(BMO.rays(beam))))
    end
end

@testset "Render" begin
    axis = nothing
    cube = BMO.CubeMesh(1)
    @test_throws BMO.MissingBackendError render!(axis, cube)
    @test_throws BMO.MissingBackendError BMO.get_view(axis)
    @test_throws BMO.MissingBackendError BMO.set_view(axis, [1 1; 0 0])
    @test_throws BMO.MissingBackendError BMO.hide_axis(axis, true)
end

@testset "Point-spread functions" begin
    # parameters for an almost thin-lens
    l = 1e-3
    R1 = 100e-3
    R2 = Inf
    d = 25.4e-3
    n = 1.5
    λ = 1e-6    
    D = 15e-3
    num_rays = 1000
    
    # plane wave source
    cs = UniformDiscSource([0, -10e-3, 0], [0, 1, 0], D, λ; num_rays)
    
    # test lens
    lens = SphericalLens(R1, R2, l, d, x -> n)
    
    # PSF detector
    psfd = PSFDetector(10e-3)
    translate3d!(psfd, [0, 200e-3 + 0.13e-3, 0])

    @testset "Airy-disc test" begin
        # build system and solve it
        sys = System([lens, psfd])
        solve_system!(sys, cs)

        # get PSF
        x, y, I_num = intensity(psfd; n=500, crop_factor=5, center=:bbox)

        # get numerical first zero of the Airy-disk
        ci = argmax(I_num)
        ix_ctr, jx_ctr = Tuple(ci)
        num_min = x[argmin(I_num[:, jx_ctr])]   # first zero through the centre column

        # theoretical Airy-disk 1st zero
        airy_min = 1.22*λ*200e-3/D

        @test abs(num_min) ≈ airy_min rtol=1e-2
    end

    @testset "PSFDetector reset" begin
        # Test if data in psfd
        @test length(psfd.data) == num_rays
        # Empty psfd and test
        empty!(psfd)
        @test isempty(psfd.data)
    end
end

@testset "Bug fixes" begin
    @testset "Issue#14" begin
        pd_res = 1000
        pd_size = 10mm
        pd = Detector(pd_size)
        BMO.zrotate3d!(pd, deg2rad(90))
        BMO.translate3d!(pd, [0.46, 0, 0])
        # Setup beam
        y_0 = 0.2
        beam = BMO.GaussianBeamlet(
            [0, y_0, 0], [0.46, -y_0, 0], 532e-9, 2.5e-3, P0 = 10e-3)
        # Solve system
        system = BMO.System([pd])
        empty!(pd)
        BMO.solve_system!(system, beam)

        pd_pwr = optical_power(
            pd;
            # restore pre-v0.11 behavior, e.g. no autolims
            n=pd_res,
            x_min=-pd_size/2,
            x_max=pd_size/2,
            z_min=-pd_size/2,
            z_max=pd_size/2
        )

        @test pd_pwr≈10e-3 atol=1e-5
    end

    @testset "Issue#22 and Issue#23" begin
        # https://github.com/JuliaPhysics/BeamletOptics.jl/issues/22
        # https://github.com/JuliaPhysics/BeamletOptics.jl/issues/23
        
        mutable struct TestSubstrate{T, S <: BMO.AbstractShape{T}, N} <: BMO.AbstractRefractiveOptic{T, S, N}
            const shape::S
            n::N
        end
        
        BMO.refractive_index(ts::TestSubstrate, ::Real) = ts.n
        set_index(ts::TestSubstrate, new) = (ts.n = new)
        get_index(ts::TestSubstrate) = ts.n
        
        "Shifts the phase of the beamlet by a specific amount in [rad]."
        function shift_phase(gb::BMO.GaussianBeamlet, phase::Real)
            BMO.electric_field!(gb, BMO.electric_field(gb) * exp(im*phase))
            return nothing
        end
        
        ref_signal(ϕ, A) = (cos(ϕ)+1)/2 * A
        
        # setup system for tests below
        splitter = CubeBeamsplitter(10mm, n->1)
        substrate_length = 10mm
        substrate = TestSubstrate(BMO.CylinderSDF(5mm, substrate_length/2), 1.5)
        
        pd_size = 10mm
        pd_res = 250
        detector = Detector(pd_size)
        
        translate3d!(substrate, [0, -25mm, 0])
        translate3d!(detector, [0, 40mm, 0])
        
        system = System([substrate, splitter, detector])
        
        start_offset = 50mm
        
        @testset "Testing electric_field calculation - non-imaging ref. index change" begin
            indices = (1, 10, 100, 1000)
            for index in indices
                set_index(substrate, index)
                phi = LinRange(0, 2pi, 30)
                int = zeros(length(phi))
                for (i, p) in enumerate(phi)
                    gb_prb = GaussianBeamlet([0, -start_offset, 0], [0, 1, 0], 1e-6, .5mm)
                    gb_ref = GaussianBeamlet([start_offset, 0, 0], [-1, 0, 0], 1e-6, .5mm)
                    empty!(detector)
                    shift_phase(gb_ref, p)
                    solve_system!(system, gb_prb)
                    solve_system!(system, gb_ref)
                    int[i] = BMO.optical_power(detector)
                end
                @test isapprox(BMO.visibility(int), 1, atol=1e-2)
            end
        end
        
        @testset "Testing electric_field calculation - ref. index based phase shift" begin
            λ = 1e-6
            gb_prb = GaussianBeamlet([0, -start_offset, 0], [0, 1, 0], λ, .5mm)
            gb_ref = GaussianBeamlet([start_offset, 0, 0], [-1, 0, 0], λ, .5mm)
            
            n_lambdas = substrate_length / BMO.wavelength(gb_prb)
            
            n_factors = LinRange(0, 1, 50)
            pwr = zeros(length(n_factors))
            # Increase ref. index of substrate until one additional λ of OPL has been introduced
            for (i, nf) in enumerate(n_factors)
                set_index(substrate, 1 + 1/n_lambdas * nf)        
                empty!(detector)
                solve_system!(system, gb_prb)
                solve_system!(system, gb_ref)
                pd_pwr = optical_power(
                    detector;
                    # restore pre-v0.11 behavior, e.g. no autolims
                    n=pd_res,
                    x_min=-pd_size/2,
                    x_max=pd_size/2,
                    z_min=-pd_size/2,
                    z_max=pd_size/2
                )
                @test isapprox(pd_pwr, ref_signal(2pi*nf, 2e-3), atol=1e-8)
            end
            # Test if opl difference is indeed one λ
            delta = BMO.optical_path_length(gb_prb)
            delta -= BMO.optical_path_length(gb_ref)
            delta /= λ
            @test delta ≈ 1
        end
        
        @testset "Testing electric_field mutation during retracing" begin    
            gb_prb = GaussianBeamlet([0, -start_offset, 0], [0, 1, 0], 1e-6, .5mm)
            gb_ref = GaussianBeamlet([start_offset, 0, 0], [-1, 0, 0], 1e-6, .5mm)    
            phis = LinRange(0, 2pi, 50)
            pwr = zeros(length(phis))
            # Vary starting phase by 0...2pi via retracing    
            for (i, phi) in enumerate(phis)
                empty!(detector)
                solve_system!(system, gb_prb)
                solve_system!(system, gb_ref)
                pd_pwr = optical_power(
                    detector;
                    # restore pre-v0.11 behavior, e.g. no autolims
                    n=pd_res,
                    x_min=-pd_size/2,
                    x_max=pd_size/2,
                    z_min=-pd_size/2,
                    z_max=pd_size/2
                )
                @test isapprox(pd_pwr, ref_signal(phi, 2e-3), atol=1e-8)
                shift_phase(gb_prb, step(phis))
            end
        end
    end
end

@testset "Aqua" begin
    using Aqua
    Aqua.test_all(BMO)
end
