module TestAsphericalLenses

using BeamletOptics
using Test
using GeometryBasics
using LinearAlgebra: dot

const BMO = BeamletOptics

const mm = 1e-3

@testset "Aspherical Lenses" begin
    @testset "Testing type definitions" begin
        @test isdefined(BMO, :ConvexAsphericalSurfaceSDF)
        @test isdefined(BMO, :ConcaveAsphericalSurfaceSDF)
    end

    @testset "Testing aspherical lens SDFs" begin
        # Test that the SDF in combination with ray-marching correctly approximates the
        # aspheric surface. Let's use the Thorlabs AL50100J lens for this test case.

        # radius
        R = 50.3583mm
        # conic constant
        k = -0.789119
        # even aspheric coefficients
        A = [0, 2.10405e-7 * (1e3)^3, 1.76468e-11 * (1e3)^5, 1.02641e-15 * (1e3)^7]
        # center thickness
        ct = 10.2mm
        # diameter
        d = 50mm
        # refractive index of BK-7 @ 1310 nm (design wavelength)
        n = 1.5036

        lens = Lens(
            EvenAsphericalSurface(R, d, k, A),
            ct,
            x -> n
        )

        system = System(lens)

        surf_errors = zeros(100)

        for (i, z) in enumerate(range(-0.02, 0.02, 100))
            ray = Ray(Point3(0.0, -0.1, z), Point3(0.0, 1.0, 0))
            beam = Beam(ray)
            solve_system!(system, beam, r_max = 40)

            surf_errors[i] = (position(beam.rays[begin]) + length(beam.rays[begin]) .* BMO.direction(beam.rays[begin]))[2] -
                             BMO.aspheric_equation(ray.pos[3], 1 / R, k, A)
        end

        # FIXME: The atol is actually derived from the raymarching epsilon. If this is puts
        # into a configurable option, this should be changed as well.
        @test all(x -> isapprox(x, 0.0; atol = 1e-10), surf_errors)

        # test if the working distance is correct
        ray = Ray([0.0, -0.1, 0.02], [0.0, 1.0, 0])
        beam = Beam(ray)
        solve_system!(system, beam, r_max = 40)

        dist = -beam.rays[end].pos[3] / beam.rays[end].dir[3]
        α = asind(beam.rays[end].dir[3])
        wd = cosd(α) * dist

        @test wd≈93.2mm atol=1e-4
    end

    @testset "Near-axis rays through hyperbolic plano-concave and plano-convex lenses" begin
        # A plano-conic lens with k = -n² images a collimated beam that enters through the plane
        # face exactly onto a point: every outgoing ray passes through the (real or virtual)
        # focus at `f` from the conic vertex. Near the vertex the concave SDF is a sliver thinner
        # than the finite-difference step of `numeric_gradient`, which used to give wrong normals
        # for |h| ≲ 10 µm.
        n = 1.458
        f = 50mm
        ct = 4mm
        d = 30mm
        concave = Lens(CircularFlatSurface(d),
            EvenAsphericalSurface((n - 1) * f, d, -n^2, [0.0]), ct, x -> n)
        convex = Lens(CircularFlatSurface(d),
            EvenAsphericalSurface(-(n - 1) * f, d, -n^2, [0.0]), 2ct, x -> n)
        # The conic vertex lies on the back face, `thickness` along the lens axis
        for (lens, thickness, f_signed) in ((concave, ct, -f), (convex, 2ct, f)),
            tilt in (0.0, deg2rad(0.5))

            xrotate3d!(lens, tilt)
            system = System(lens)
            e_x, axis = BMO.orientation(lens)[:, 1], BMO.orientation(lens)[:, 2]
            focus = position(lens) .+ (thickness + f_signed) .* axis
            # Below h ≈ 0.1 µm the sag (h²/2R ≈ 1e-13 m) is smaller than the ray-marching
            # tolerance, the plane face is hit instead; the angle error stays below h / f.
            for h in (1e-6, 1e-5, 1e-4, 1e-3, 5mm)
                beam = Beam(Ray(position(lens) .- 0.1 .* axis .+ h .* e_x, Vector(axis)))
                solve_system!(system, beam)
                p, out = position(last(beam.rays)), BMO.direction(last(beam.rays))
                # lateral miss of the focus by the (extended) outgoing ray, relative to h
                t = dot(focus .- p, axis) / dot(out, axis)
                miss = dot(p .+ t .* out .- focus, e_x)
                @test abs(miss) < 1e-6 * h
            end
            xrotate3d!(lens, -tilt)
        end
    end

    @testset "Complex aspherical imaging system" begin
        # setup system
        L1 = Lens(
            EvenAsphericalSurface(
                1.054mm, # r
                1.333024mm, # d
                -0.14294, # conic
                [0, 0.038162 * (1e3)^3, 0.06317 * (1e3)^5,
                    -0.020792 * (1e3)^7, 0.18432 * (1e3)^9,
                    -0.04827 * (1e3)^11, 0.094529 * (1e3)^13] # coeffs
            ),
            EvenAsphericalSurface(
                2.027mm, # r
                1.216472mm, # d
                8.0226, # conic
                [0, 0.0074974 * (1e3)^3, 0.064686 * (1e3)^5,
                    0.19354 * (1e3)^7, -0.50703 * (1e3)^9,
                    -0.34529 * (1e3)^11, 5.9938 * (1e3)^13] # coeffs
            ),
            0.72mm, # center thickness
            n -> 1.580200
        )

        L2 = Lens(
            EvenAsphericalSurface(
                -3.116mm, # r
                1.4mm, # d
                -49.984, # conic
                [0, -0.31608 * (1e3)^3, 0.34755 * (1e3)^5,
                    -0.17102 * (1e3)^7, -0.41506 * (1e3)^9,
                    -1.342 * (1e3)^11, 5.0594 * (1e3)^13, -2.7483 * (1e3)^15] # coeffs
            ),
            EvenAsphericalSurface(
                -4.835mm, # r
                1.9mm, # d
                1.6674, # conic
                [0, -0.079727 * (1e3)^3, 0.13899 * (1e3)^5, -0.044057 * (1e3)^7,
                    -0.019369 * (1e3)^9, 0.016993 * (1e3)^11, 0.093716 * (1e3)^13,
                    -0.080329 * (1e3)^15] # coeffs
            ),
            0.55mm, # center_thickness
            n -> 1.804700
        )

        translate3d!(L2, [0, thickness(L1) + 0.39mm, 0])

        L3 = Lens(
            EvenAsphericalSurface(
                3.618mm, # r
                3.04mm, # d
                -44.874, # conic
                [0, -0.14756 * (1e3)^3, 0.035194 * (1e3)^5, -0.0032262 * (1e3)^7,
                    0.0018592 * (1e3)^9, 0.00036658 * (1e3)^11, -0.00016039 * (1e3)^13,
                    -3.1846e-5 * (1e3)^15] # coeffs
            ),
            EvenAsphericalSurface(
                2.161mm, # r
                3.7mm, # d
                -10.719, # conic
                [0, -0.096568 * (1e3)^3, 0.026771 * (1e3)^5, -0.011261 * (1e3)^7,
                    0.0019879 * (1e3)^9, 0.00015579 * (1e3)^11, -0.00012433 * (1e3)^13,
                    1.5264e-5 * (1e3)^15] # coeffs
            ),
            0.7mm, # center_thickness
            n -> 1.580200
        )

        translate_to3d!(L3, position(L2))
        translate3d!(L3, [0, thickness(L2) + 0.63mm, 0])

        Filt = Lens(
            CircularFlatSurface(4.2mm),
            0.15mm,
            n -> 1.516800
        )

        translate_to3d!(Filt, position(L3))
        translate3d!(Filt, [0, thickness(L3) + 0.19mm, 0])

        Cover = Lens(
            CircularFlatSurface(4.9mm),
            0.5mm,
            n -> 1.469200
        )
        translate_to3d!(Cover, position(Filt))
        translate3d!(Cover, [0, thickness(Filt) + 0.18mm, 0])

        # test thickness
        @test thickness(L1) ≈ 0.72mm
        @test thickness(L2) ≈ 0.55mm
        @test thickness(L3) ≈ 0.7mm
        @test thickness(Filt) ≈ 0.15mm
        @test thickness(Cover) ≈ 0.5mm

        system = System([L1, L2, L3, Filt, Cover])

        # 0° beams
        beams = [
            Beam([0, -0.5mm, -1.3mm / 2], [0, 1, 0], 0.5876e-6),
            Beam([0, -0.5mm, 0], [0, 1, 0], 0.5876e-6),
            Beam([0, -0.5mm, 1.3mm / 2], [0, 1, 0], 0.5876e-6)
        ]
        for beam in beams
            solve_system!(system, beam, r_max = 50)
            f_pos = last(beam.rays).pos + 0.12mm * last(beam.rays).dir

            # test if the beam is correctly focussed
            @test length(BMO.rays(beam)) == 11
            @test f_pos[3]≈0 atol=1e-7
        end

        # test rings
        @test 2*L1.shape.sdfs[4].hthickness ≈ 0.00060839 atol=1e-6
        @test 2*L2.shape.sdfs[4].hthickness ≈ 0.00057497 atol=1e-6
        @test 2*L3.shape.sdfs[4].hthickness ≈ 0.00048395 atol=1e-6
    end
end

end # MODULE