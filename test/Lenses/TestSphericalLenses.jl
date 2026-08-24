module TestSphericalLenses

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics

const mm = 1e-3

@testset "Spherical Lenses" begin
    @testset "Testing type definitions" begin
        @test isdefined(BMO, :AbstractSDF)
        @test isdefined(BMO, :SphereSDF)
        @test isdefined(BMO, :CylinderSDF)
        @test isdefined(BMO, :CutSphereSDF)
        @test isdefined(BMO, :ThinLensSDF)
    end

    @testset "Thin lens focal length" begin
        # define thin lens
        R1 = 1
        R2 = 1
        nl = 1.5
        tl = BMO.ThinLensSDF(R1, R2, 0.1)
        translate3d!(tl, [0, -thickness(tl) / 2, 0])
        p = Lens(tl, x -> 1.5)
        system = System(p)

        # compare numerical and analytical focal length
        f_analytical = BMO.lensmakers_eq(R1, -R2, nl)
        zs = -0.04:0.01:0.04
        for (i, z) in enumerate(zs)
            # skip optical axis ray
            if z ≈ 0
                continue
            end
            xs = 0.1:0.1:1.5
            df = zeros(Float64, length(xs))
            ray = Ray([0, -0.5, z], [0, 1, 0], 1e3)
            beam = Beam(ray)
            solve_system!(system, beam)
            # test if numerical and analytical focal length agree
            for (i, x) in enumerate(xs)
                df[i] = BMO.line_point_distance3d(beam.rays[end], [0, x, 0])
            end
            @test xs[findmin(df)[2]] ≈ f_analytical
        end
    end

    @testset "Testing lens constructor" begin
        # Test against Thorlab spherical lenses
        r1 = 34.9mm
        r2 = -r1
        l = 6.8mm
        LB1811 = SphericalLens(r1, r2, l)
        @test typeof(BMO.shape(LB1811)) <: BMO.UnionSDF
        @test thickness(BMO.shape(LB1811)) == l
        r1 = Inf
        r2 = -15.5mm
        l = 8.6mm
        LA1805 = SphericalLens(r1, r2, l)
        @test typeof(BMO.shape(LA1805)) <: BMO.UnionSDF
        @test thickness(BMO.shape(LA1805)) == l
        r1 = -52.0mm
        r2 = -r1
        l = 3mm
        LD1464 = SphericalLens(r1, r2, l)
        @test typeof(BMO.shape(LD1464)) <: BMO.UnionSDF
        @test thickness(BMO.shape(LD1464)) == l
        r1 = Inf
        r2 = 25.7mm
        l = 3.5mm
        LC1715 = SphericalLens(r1, r2, l)
        @test typeof(BMO.shape(LC1715)) <: BMO.UnionSDF
        @test thickness(BMO.shape(LC1715)) == l
        r1 = -82.2mm
        r2 = -32.1mm
        l = 3.6mm
        LE1234 = SphericalLens(r1, r2, l)
        @test typeof(BMO.shape(LE1234)) <: BMO.UnionSDF
        @test thickness(BMO.shape(LE1234)) == l
    end

    """Test coma for rotated and translated optical system"""
    function test_coma(ray::BMO.AbstractRay, f0::AbstractArray,
            dir::AbstractArray; atol = 7e-5)
        t = BMO.line_plane_distance3d(f0, dir, position(ray), BMO.direction(ray))
        p0 = position(ray) + t * BMO.direction(ray)
        dz = norm(p0 - f0)
        if dz ≤ atol
            return true
        else
            return error("Coma dz=$dz larger than atol=$atol")
        end
    end

    @testset "Testing doublet lenses" begin
        # Define refractive index functions
        λs = [488e-9, 707e-9, 1064e-9]
        NLAK22 = DiscreteRefractiveIndex(λs, [1.6591, 1.6456, 1.6374])
        NSF10 = DiscreteRefractiveIndex(λs, [1.7460, 1.7168, 1.7021])

        function test_doublet(λ, bfl, δf)
            # Thorlabs lens from https://www.thorlabs.com/thorproduct.cfm?partnumber=AC254-150-AB
            AC254_150_AB = SphericalDoubletLens(
                87.9mm, -105.6mm, Inf, 6mm, 3mm, BMO.inch, NLAK22, NSF10)
            # Rotate and translate to test lens kinematics
            translate3d!(AC254_150_AB, [0.05, 0.05, 0.05])
            xrotate3d!(AC254_150_AB, deg2rad(-60))
            zrotate3d!(AC254_150_AB, deg2rad(45))
            # Define system
            system = System([AC254_150_AB])
            # Define semi-diameter for lens ray bundle, selected for min. spherical aberrations
            z0 = 5mm
            zs = LinRange(-z0, z0, 30)
            fs = similar(zs)
            # Beam spawn point
            dir = -orientation(AC254_150_AB.back.shape)[:, 2]       # rotated collimated ray direction
            pos = position(AC254_150_AB.front.shape) + 0.05 * dir  # rotated collimated ray position
            nv = BMO.normal3d(dir)                                     # orthogonal to moved system optical axis
            beam = Beam(pos, -dir, λ)
            # Calculate equivalent back focal length point
            f_z = thickness(AC254_150_AB) + bfl + δf
            f0 = position(AC254_150_AB.front.shape) + f_z * -dir
            for (i, z) in enumerate(zs)
                beam = Beam(pos + z * nv, -dir, λ)
                solve_system!(system, beam)
                @test length(BMO.rays(beam)) == 4
                @test BMO.refractive_index.(beam.rays) ==
                      [1, NLAK22(λ), NSF10(λ), 1]
                fs[i] = test_coma(last(BMO.rays(beam)), f0, dir, atol = 1e-6)
            end
            # Test center ray normal vectors
            beam = Beam(pos + 0 * nv, -dir, λ)
            solve_system!(system, beam)
            for i in 1:(length(beam.segments) - 1)
                @test abs(dot(BMO.normal3d(beam.segments[i]), direction(beam.segments[i]))) ≈ 1
            end
            return true
        end

        # Run tests for AC254_150_AB against plot data at https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=12767
        @test test_doublet(488e-9, 143.68mm, -2.064e-4)
        @test test_doublet(707e-9, 143.68mm, 0)
        @test test_doublet(1064e-9, 143.68mm, +7.466e-4)

        @testset "Polarized ray through doublet" begin
            λ = 707e-9
            AC254 = SphericalDoubletLens(87.9mm, -105.6mm, Inf, 6mm, 3mm, BMO.inch, NLAK22, NSF10)
            system = System([AC254])
            ray = PolarizedRay([0.0, -0.05, 0.0], [0.0, 1.0, 0.0], λ, [1.0, 0.0, 0.0])
            beam = Beam(ray)
            solve_system!(system, beam)
            @test length(BMO.rays(beam)) == 4
            @test BMO.refractive_index.(beam.rays) == [1.0, NLAK22(λ), NSF10(λ), 1.0]
            @test optical_path_length(beam) > length(beam)
        end

        @testset "Gaussian and AGB beamlet through doublet" begin
            λ = 707e-9
            AC254 = SphericalDoubletLens(87.9mm, -105.6mm, Inf, 6mm, 3mm, BMO.inch, NLAK22, NSF10)
            system = System([AC254])
            
            # Gaussian beamlet
            gauss = GaussianBeamlet([0.0, -0.05, 0.0], [0.0, 1.0, 0.0], λ, 1.0mm)
            solve_system!(system, gauss)
            @test length(BMO.rays(gauss.chief)) == 4
            @test optical_path_length(gauss) > 0.0

            # Astigmatic Gaussian beamlet
            agb = AstigmaticGaussianBeamlet([0.0, -0.05, 0.0], [0.0, 1.0, 0.0], λ, 1.0mm)
            solve_system!(system, agb; check_invariant = true)
            @test length(BMO.rays(agb.c)) == 4
            @test optical_path_length(agb) > 0.0
            @test BMO.check_optical_invariant(agb, 4)
        end

        @testset "Reverse and index-matched doublet propagation" begin
            λ = 707e-9
            # Reverse propagation from +y to -y
            AC254 = SphericalDoubletLens(87.9mm, -105.6mm, Inf, 6mm, 3mm, BMO.inch, NLAK22, NSF10)
            sys_rev = System([AC254])
            beam_rev = Beam([0.0, 0.05, 0.0], [0.0, -1.0, 0.0], λ)
            solve_system!(sys_rev, beam_rev)
            @test length(BMO.rays(beam_rev)) == 4
            @test BMO.refractive_index.(beam_rev.rays) == [1.0, NSF10(λ), NLAK22(λ), 1.0]

            # Index-matched doublet (n1 == n2)
            AC_matched = SphericalDoubletLens(87.9mm, -105.6mm, Inf, 6mm, 3mm, BMO.inch, 1.5, 1.5)
            sys_matched = System([AC_matched])
            beam_matched = Beam([0.0, -0.05, 0.0], [0.0, 1.0, 0.0], λ)
            solve_system!(sys_matched, beam_matched)
            @test length(BMO.rays(beam_matched)) == 4
            @test BMO.refractive_index.(beam_matched.rays) == [1.0, 1.5, 1.5, 1.0]

            # Water-immersed doublet
            water = IsotropicMedium(:Water, 1.333)
            sys_water = System([AC254], water)
            @test BMO.ambient_medium(sys_water) === water
            @test BMO.refractive_index(sys_water, λ) ≈ 1.333
            ray_water = Ray([0.0, -0.05, 0.0], [0.0, 1.0, 0.0], λ, 1.333)
            beam_water = Beam(ray_water)
            solve_system!(sys_water, beam_water)
            @test length(BMO.rays(beam_water)) == 4
            @test BMO.refractive_index.(beam_water.rays) == [1.333, NLAK22(λ), NSF10(λ), 1.333]
        end
    end
end

end # MODULE