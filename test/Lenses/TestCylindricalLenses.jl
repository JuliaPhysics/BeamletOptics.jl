module TestCylindricalLenses

using BeamletOptics
using Test

const BMO = BeamletOptics

const mm = 1e-3

@testset "Cylindrical Lenses" begin
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
    @testset "Convex/cav cylinder lenses" begin
        # Thorlabs LJ1878L2, plano-convex
        r = 5.2mm
        d = 10mm
        h = 20mm
        ct = 5.9mm
        lens = Lens(
            CylindricalSurface(r, d, h),
            ct,
            n -> 1.517
        )

        # test lens thickness
        @test thickness(lens) ≈ ct
        # test edge thickness
        @test thickness(BMO.shape(lens).sdfs[1])≈2.12mm atol=1e-4
        # test back focal length
        @test working_distance(lens, 0.05 * d / 2)≈6.1mm atol=1e-4

        # Thorlabs LK1900L1, plano-concave
        r = -13.1mm
        d = 16mm
        h = 18mm
        ct = 2.0mm
        lens = Lens(
            CylindricalSurface(r, d, h),
            ct,
            n -> 1.517
        )

        # test lens thickness
        @test thickness(lens) ≈ ct
    end

    @testset "Convex/cav acylinder lenses" begin
        # Thorlabs AYL2520, plano-convex
        radius = 15.538mm
        diameter = 25mm
        height = 50mm
        conic_constant = -1.0
        ct = 7.5mm
        lens = Lens(
            BMO.AcylindricalSurface(
                radius,
                diameter,
                height,
                conic_constant,
                [0, 1.1926075e-5*(1e3)^3, -2.9323497e-9*(1e3)^5, -1.8718889e-11*(1e3)^7, -1.7009961e-14*(1e3)^9, 3.5481542e-17*(1e3)^11, 6.5241296e-20*(1e3)^13]
            ),
            ct,
            n -> 1.777
        )

        # test lens thickness
        @test thickness(lens) ≈ ct

        # test back focal length
        @test working_distance(lens, 0.05 * diameter / 2) ≈ 15.8mm atol=1e-4

        # Thorlabs AYL2520, plano-convex (inverted)
        radius = -15.538mm
        diameter = 25mm
        height = 50mm
        conic_constant = -1.0
        ct = 7.5mm
        lens = Lens(
            BMO.AcylindricalSurface(
                radius,
                diameter,
                height,
                conic_constant,
                [0, 1.1926075e-5*(1e3)^3, -2.9323497e-9*(1e3)^5, -1.8718889e-11*(1e3)^7, -1.7009961e-14*(1e3)^9, 3.5481542e-17*(1e3)^11, 6.5241296e-20*(1e3)^13]
            ),
            ct,
            n -> 1.777
        )

        # test lens thickness
        @test thickness(lens) ≈ ct
    end
end

end # MODULE