module TestSurfaces

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics

const mm = 1e-3

@testset "Surfaces" begin
    function working_distance(lens, offset_z)
        ray = Ray([0.0, -1.0, offset_z], [0.0, 1.0, 0])
        system = System([lens])
        beam = Beam(ray)
        solve_system!(system, beam)

        dist = -beam.rays[end].pos[3] / beam.rays[end].dir[3]
        α = asind(beam.rays[end].dir[3])
        wd = cosd(α) * dist

        return wd
    end

    @testset "Lens construction from surfaces" begin
        ## Thorlabs LA4052, plano-convex
        r1 = 16.1mm
        r2 = Inf
        l = 8.2mm
        d = 25.4mm
        lens = Lens(
            SphericalSurface(r1, d),
            l,
            n -> 1.458
        )

        # test lens thickness
        @test thickness(lens) ≈ l
        # test edge thickness
        @test thickness(BMO.shape(lens).sdfs[1])≈2mm atol=1e-4
        # test back focal length
        @test working_distance(lens, 0.05 * d / 2)≈29.5mm atol=1e-4

        ## Thorlabs LB1761, bi-convex
        r1 = 24.5mm
        r2 = -24.5mm
        l = 9.0mm
        d = 25.4mm
        lens = Lens(
            SphericalSurface(r1, d),
            SphericalSurface(r2, d),
            l,
            n -> 1.517
        )
        shape = lens.shape

        # test lens thickness
        @test thickness(lens) ≈ l
        # test edge thickness
        @test thickness(shape.sdfs[1])≈1.9mm atol=1e-4
        # test back focal length
        @test working_distance(lens, 0.05 * d / 2)≈22.2mm atol=1mm

        ## Thorlabs LC1715, plano-concave
        r1 = Inf
        r2 = 25.7mm
        l = 3.5mm
        d = 25.4mm
        lens = Lens(
            SphericalSurface(r1, d),
            SphericalSurface(r2, d),
            l,
            n -> 1.517
        )
        shape = lens.shape

        # test lens thickness
        @test thickness(lens) ≈ l
        # test edge thickness
        @test BMO.sag(shape.sdfs[2]) +
              thickness(shape.sdfs[1])≈0.006858 atol=1e-4

        ## Thorlabs LD2297, bi-concave
        r1 = -39.6mm
        r2 = 39.6mm
        l = 3.0mm
        d = 25.4mm
        lens = Lens(
            SphericalSurface(r1, d),
            SphericalSurface(r2, d),
            l,
            n -> 1.517
        )
        shape = lens.shape

        # test lens thickness
        @test thickness(lens) ≈ l

        @test BMO.sag(shape.sdfs[2]) + BMO.sag(shape.sdfs[3]) +
              thickness(shape.sdfs[1])≈0.0072 atol=1e-4

        ## Thorlabs LBF254-040, best-form
        r1 = 134.6mm
        r2 = -24.0mm
        l = 6.5mm
        d = 25.4mm
        lens = Lens(
            SphericalSurface(r1, d),
            SphericalSurface(r2, d),
            l,
            n -> 1.517
        )
        shape = lens.shape

        # test lens thickness
        @test thickness(lens) ≈ l

        # test edge thickness
        @test thickness(shape.sdfs[1])≈2.286mm atol=1e-4

        ## Thorlabs LE1234, positive meniscus
        r1 = -82.2mm
        r2 = -32.1mm
        l = 3.6mm
        d = 25.4mm
        lens = Lens(
            SphericalSurface(r1, d),
            SphericalSurface(r2, d),
            l,
            n -> 1.517
        )
        shape = lens.shape

        # test lens thickness
        @test thickness(lens) ≈ l

        # test edge thickness
        @test thickness(shape.sdfs[1]) + BMO.sag(shape.sdfs[2])≈2mm atol=1e-4

        ## Thorlabs LF1822, negative meniscus
        r1 = -33.7mm
        r2 = -100.0mm
        l = 3.0mm
        d = 25.4mm
        lens = Lens(
            SphericalSurface(r1, d),
            SphericalSurface(r2, d),
            l,
            n -> 1.517
        )
        shape = lens.shape

        # test lens thickness
        @test thickness(lens) ≈ l

        # test edge thickness
        @test thickness(shape.sdfs[1]) +
              BMO.sag(shape.sdfs[2])≈4.7mm atol=1e-4

        ## Generic "true" meniscus
        r1 = 103.4371mm
        r2 = 61.14925mm
        l = 1.5mm
        d = 55mm
        lens = Lens(
            SphericalSurface(r1, d),
            SphericalSurface(r2, d),
            l,
            n -> 1.517
        )
        shape = lens.shape

        # test lens thickness
        @test thickness(lens) ≈ l

        # test ring generation
        NBK7 = DiscreteRefractiveIndex([532e-9, 1064e-9], [1.5195, 1.5066])

        s1 = Lens(
            SphericalSurface(38.184mm, 2*1.840mm, 2*2.380mm),
            SphericalSurface(3.467mm, 2*2.060mm, 2*2.380mm),
            0.5mm,
            NBK7
        )

        @test 2*s1.shape.sdfs[5].hthickness ≈ 0.001134 atol=1e-6

        s2 = Lens(
            SphericalSurface(3.467mm, 2*2.060mm, 2*2.380mm),
            SphericalSurface(-5.020mm, 2*2.380mm, 2*2.380mm),
            2.5mm,
            NBK7
        )

        @test 2*s2.shape.sdfs[4].hthickness ≈ 0.001221590 atol=1e-8

        # doublet test case
        s1 = SphericalSurface(7.744mm, 2*2.812mm, 2*3mm)
        s2 = SphericalSurface(-3.642mm, 2*3mm)
        s3 = SphericalSurface(-14.413mm, 2*2.812mm, 2*3mm)

        dl21 = Lens(s1, s2, 3.4mm, NBK7)
        dl22 = Lens(s2, s3, 1.0mm, NBK7)

        @test 2*dl21.shape.sdfs[4].hthickness ≈ 0.001294398 atol=1e-6
        @test 2*dl22.shape.sdfs[4].hthickness ≈ 0.000723025 atol=1e-6
    end

    @testset "Doublet and triplet construction from surfaces" begin
        d = 25.4mm
        hs = (0.0, 1e-6, 1mm, 5mm, -7mm, 10mm)
        n1 = λ -> 1.5
        n2 = λ -> 1.7
        n3 = λ -> 1.6
        # Even aspheric coefficients for r², r⁴ and r⁶ in SI units
        A = [0, 2e-7 * (1e3)^3, 1e-11 * (1e3)^5]

        """Trace a ray along +y at height `h` and extrapolate the last ray to y = 0.1 m."""
        function trace(object, h)
            beam = Beam(Ray([0, -0.05, h], [0, 1, 0], 1e-6))
            solve_system!(System([object]), beam)
            r = last(rays(beam))
            p = position(r)
            dir = direction(r)
            return p + (0.1 - p[2]) / dir[2] * dir, dir
        end

        function test_equivalence(a, b, hs, atol)
            for h in hs
                p_a, d_a = trace(a, h)
                p_b, d_b = trace(b, h)
                @test norm(p_a - p_b) ≤ atol
                @test norm(d_a - d_b) ≤ atol
            end
        end

        @testset "Spherical surfaces" begin
            # must reproduce the spherical constructors
            dl = DoubletLens(SphericalSurface(33.3mm, d), SphericalSurface(-22.3mm, d),
                SphericalSurface(-291.1mm, d), 9mm, 2.5mm, n1, n2)
            test_equivalence(dl,
                SphericalDoubletLens(33.3mm, -22.3mm, -291.1mm, 9mm, 2.5mm, d, n1, n2), hs, 1e-12)
            @test thickness(dl) ≈ 11.5mm

            # plano surface as CircularFlatSurface
            dl = DoubletLens(SphericalSurface(87.9mm, d), SphericalSurface(-105.6mm, d),
                CircularFlatSurface(d), 6mm, 3mm, n1, n2)
            test_equivalence(dl,
                SphericalDoubletLens(87.9mm, -105.6mm, Inf, 6mm, 3mm, d, n1, n2), hs, 1e-12)

            tl = TripletLens(SphericalSurface(50mm, d), SphericalSurface(-40mm, d),
                SphericalSurface(100mm, d), SphericalSurface(-60mm, d), 8mm, 3mm, 6mm, n1, n2, n3)
            test_equivalence(tl,
                SphericalTripletLens(50mm, -40mm, 100mm, -60mm, 8mm, 3mm, 6mm, d, n1, n2, n3),
                hs, 1e-12)
            @test thickness(tl) ≈ 17mm
        end

        @testset "Equal-index invariant" begin
            # with identical ref. indices the cemented surfaces must not change the ray path
            s_back = SphericalSurface(-200mm, d)

            # aspherical outer surface
            s1 = EvenAsphericalSurface(50mm, d, -0.8, A)
            dl = DoubletLens(s1, SphericalSurface(-40mm, d), s_back, 8mm, 3mm, n1, n1)
            test_equivalence(dl, Lens(s1, s_back, 11mm, n1), hs, 1e-9)

            # aspherical cemented surface, concave and convex towards the front element
            for s2 in (EvenAsphericalSurface(-40mm, d, -0.5, [0, -1e-7 * (1e3)^3]),
                EvenAsphericalSurface(60mm, d, -0.5, [0, 1e-7 * (1e3)^3]))
                dl = DoubletLens(SphericalSurface(50mm, d), s2, s_back, 6mm, 6mm, n1, n1)
                test_equivalence(dl, SphericalLens(50mm, -200mm, 12mm, d, n1), hs, 1e-9)
            end

            # triplet with aspherical outer and cemented surfaces
            s3 = EvenAsphericalSurface(100mm, d, -0.5, [0, 1e-7 * (1e3)^3])
            s4 = EvenAsphericalSurface(-60mm, d, -1.0, [0, -1e-7 * (1e3)^3])
            tl = TripletLens(s1, SphericalSurface(-40mm, d), s3, s4, 8mm, 3mm, 6mm, n1, n1, n1)
            test_equivalence(tl, Lens(s1, s4, 17mm, n1), hs, 1e-9)
            @test thickness(tl) ≈ 17mm
        end

        @testset "Different clear apertures" begin
            dl = DoubletLens(SphericalSurface(50mm, 30mm), SphericalSurface(-40mm, d),
                SphericalSurface(-200mm, 20mm), 9mm, 3mm, n1, n1)
            sl = Lens(SphericalSurface(50mm, 30mm), SphericalSurface(-200mm, 20mm), 12mm, n1)
            test_equivalence(dl, sl, (0.0, 1mm, 5mm, -7mm, 9.5mm), 1e-9)
        end

        @testset "Refractive index sequence and kinematics" begin
            s1 = EvenAsphericalSurface(50mm, d, -0.8, A)
            s2 = EvenAsphericalSurface(-40mm, d, -0.5, [0, -1e-7 * (1e3)^3])
            s3 = EvenAsphericalSurface(100mm, d, -0.5, [0, 1e-7 * (1e3)^3])
            s4 = SphericalSurface(-60mm, d)
            dl = DoubletLens(s1, s2, s4, 8mm, 3mm, n1, n2)
            tl = TripletLens(s1, s2, s3, s4, 8mm, 3mm, 6mm, n1, n2, n3)
            for (lens, n_seq) in ((dl, [1, 1.5, 1.7, 1]), (tl, [1, 1.5, 1.7, 1.6, 1]))
                translate3d!(lens, [0.05, 0.05, 0.05])
                xrotate3d!(lens, deg2rad(-60))
                zrotate3d!(lens, deg2rad(45))
                system = System([lens])
                dir = orientation(lens)[:, 2]
                pos = position(lens) - 0.05 * dir
                nv = BMO.normal3d(dir)
                for h in hs
                    beam = Beam(pos + h * nv, dir, 1e-6)
                    solve_system!(system, beam)
                    @test BMO.refractive_index.(rays(beam)) == n_seq
                end
            end
        end

        @testset "Unsupported surfaces" begin
            # cemented cylindrical surfaces are not traced correctly, hence no method exists
            c1 = CylindricalSurface(50mm, d, d)
            c2 = CylindricalSurface(-40mm, d, d)
            flat = RectangularFlatSurface(d)
            @test_throws MethodError DoubletLens(c1, c2, flat, 8mm, 3mm, n1, n2)
            @test_throws MethodError TripletLens(c1, c2, flat, flat, 8mm, 3mm, 3mm, n1, n2, n3)
        end
    end
end

end # MODULE