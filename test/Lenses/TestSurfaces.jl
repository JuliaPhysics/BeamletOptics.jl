module TestSurfaces

using BeamletOptics
using Test

const BMO = BeamletOptics

const mm = 1e-3

@testset "Surfaces" begin
    function working_distance(lens, offset_z)
        ray = Ray([0.0, -1.0, offset_z], [0.0, 1.0, 0])
        system = System([lens])
        beam = Beam(ray)
        solve_system!(system, beam)

        last_ray = beam.rays[end]
        dist = -position(last_ray)[3] / direction(last_ray)[3]
        α = asind(direction(last_ray)[3])
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
end

end # MODULE