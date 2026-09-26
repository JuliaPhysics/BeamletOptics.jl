module TestSphericalLenses

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics

const mm = 1e-3

@testset "Spherical Lenses" begin
    @testset "Mechanical ring spans the lens edge" begin
        # Edge sagitta of R = 40 mm over a clear aperture of 20 mm
        R, d, md = 40mm, 20mm, 25.4mm
        s = R - sqrt(R^2 - (d / 2)^2)
        function ring_extent(lens)
            rings = [p for p in BMO.shape(lens).sdfs if p isa BMO.RingSDF]
            @test length(rings) == 1
            r = only(rings)
            return BMO.position(r)[2] - r.hthickness, BMO.position(r)[2] + r.hthickness
        end
        n = λ -> 1.5
        # biconvex: edge from s to T - s
        lens = Lens(SphericalSurface(R, d, md), SphericalSurface(-R, d), 5mm, n)
        lo, hi = ring_extent(lens)
        @test lo ≈ s atol = 1e-9
        @test hi ≈ 5mm - s atol = 1e-9
        # biconcave: edge from -s to T + s
        lens = Lens(SphericalSurface(-R, d, md), SphericalSurface(R, d), 2mm, n)
        lo, hi = ring_extent(lens)
        @test lo ≈ -s atol = 1e-9
        @test hi ≈ 2mm + s atol = 1e-9
        # convex-plano: edge from 0 to T - s
        lens = Lens(SphericalSurface(Inf, d, md), SphericalSurface(-R, d), 5mm, n)
        lo, hi = ring_extent(lens)
        @test lo ≈ 0 atol = 1e-9
        @test hi ≈ 5mm - s atol = 1e-9
    end

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
        is = BMO.intersect3d(f0, dir, ray)
        p0 = position(ray) + length(is) * BMO.direction(ray)
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
                beam.rays[1].pos = pos + z * nv
                solve_system!(system, beam)
                @test length(BMO.rays(beam)) == 4
                @test BMO.refractive_index.(beam.rays) ==
                      [1, NLAK22(λ), NSF10(λ), 1]
                fs[i] = test_coma(last(BMO.rays(beam)), f0, dir, atol = 1e-6)
            end
            # Test center ray normal vectors
            beam.rays[1].pos = pos + 0 * nv
            solve_system!(system, beam)
            for i in 1:(length(beam.rays) - 1)
                @test abs(dot(beam.rays[i].intersection.n, beam.rays[i].dir)) ≈ 1
            end
            return true
        end
        # Run tests for AC254_150_AB against plot data at https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=12767
        @test test_doublet(488e-9, 143.68mm, -2.064e-4)
        @test test_doublet(707e-9, 143.68mm, 0)
        @test test_doublet(1064e-9, 143.68mm, +7.466e-4)
    end

    @testset "Testing doublet lenses ray and beam type compats." begin
        # Every ray/beam type the solver can produce must trace through a DoubletLens.
        # Regression: the chief ray of an AstigmaticGaussianBeamlet is a PolarizedRay,
        # which is not a subtype of Ray, so the doublet interact3d must accept AbstractRay.
        n1 = λ -> 1.5
        n2 = λ -> 1.7
        λ = 1e-6
        w0 = 1mm
        pos = [0.0, 0, 0]
        dir = [0.0, 1, 0]
        support = [1.0, 0, 0]
        n_seq = [1, n1(λ), n2(λ), 1]

        function doublet_system()
            dl = SphericalDoubletLens(33.3mm, -22.3mm, -291.1mm, 9mm, 2.5mm, 25.4mm, n1, n2)
            translate3d!(dl, [0, 50mm, 0])
            return System([dl])
        end

        @testset "Beam{Ray}" begin
            beam = Beam(pos, dir, λ)
            solve_system!(doublet_system(), beam)
            @test length(rays(beam)) == 4
            @test BMO.refractive_index.(rays(beam)) == n_seq
        end

        @testset "Beam{PolarizedRay}" begin
            beam = Beam(pos, dir, λ, [1.0, 0, 0])
            solve_system!(doublet_system(), beam)
            @test length(rays(beam)) == 4
            @test BMO.refractive_index.(rays(beam)) == n_seq
        end

        @testset "GaussianBeamlet" begin
            gauss = GaussianBeamlet(pos, dir, λ, w0; support)
            solve_system!(doublet_system(), gauss)
            @test length(rays(gauss.chief)) == 4
            @test BMO.refractive_index.(rays(gauss.chief)) == n_seq
        end

        @testset "AstigmaticGaussianBeamlet" begin
            agb = AstigmaticGaussianBeamlet(pos, dir, λ, w0; support)
            solve_system!(doublet_system(), agb)
            @test length(rays(agb.c)) == 4
            @test BMO.refractive_index.(rays(agb.c)) == n_seq
        end
    end

    @testset "Testing triplet lenses" begin
        """Extrapolate the last ray of `beam` to y = 0.1 m and return the point and direction."""
        function extrapolate_ray(beam)
            r = last(rays(beam))
            p = position(r)
            d = direction(r)
            p_end = p + (0.1 - p[2]) / d[2] * d
            return p_end, d
        end

        @testset "Equal-index invariant" begin
            tl = SphericalTripletLens(50mm, -40mm, 100mm, -60mm, 8mm, 3mm, 6mm, BMO.inch,
                λ -> 1.5, λ -> 1.5, λ -> 1.5)
            sl = SphericalLens(50mm, -60mm, 17mm, BMO.inch, λ -> 1.5)
            for h in (0.0, 1e-6, 1mm, 5mm, -7mm, 10mm)
                sys_tl = System([tl])
                sys_sl = System([sl])
                beam_tl = Beam(Ray([0, -0.05, h], [0, 1, 0], 1e-6))
                beam_sl = Beam(Ray([0, -0.05, h], [0, 1, 0], 1e-6))
                solve_system!(sys_tl, beam_tl)
                solve_system!(sys_sl, beam_sl)
                p_tl, d_tl = extrapolate_ray(beam_tl)
                p_sl, d_sl = extrapolate_ray(beam_sl)
                @test norm(p_tl - p_sl) ≤ 1e-9
                @test norm(d_tl - d_sl) ≤ 1e-9
            end
        end

        @testset "Doublet equivalence" begin
            n1 = λ -> 1.5
            n2 = λ -> 1.7
            tl = SphericalTripletLens(33.3mm, -22.3mm, -100mm, -291.1mm, 9mm, 1mm, 1.5mm,
                25.4mm, n1, n2, n2)
            dl = SphericalDoubletLens(33.3mm, -22.3mm, -291.1mm, 9mm, 2.5mm, 25.4mm, n1, n2)
            for h in (0.0, 1e-6, 1mm, 5mm, -7mm, 10mm)
                sys_tl = System([tl])
                sys_dl = System([dl])
                beam_tl = Beam(Ray([0, -0.05, h], [0, 1, 0], 1e-6))
                beam_dl = Beam(Ray([0, -0.05, h], [0, 1, 0], 1e-6))
                solve_system!(sys_tl, beam_tl)
                solve_system!(sys_dl, beam_dl)
                p_tl, d_tl = extrapolate_ray(beam_tl)
                p_dl, d_dl = extrapolate_ray(beam_dl)
                @test norm(p_tl - p_dl) ≤ 1e-9
                @test norm(d_tl - d_dl) ≤ 1e-9
            end
        end

        @testset "Kinematics" begin
            tl = SphericalTripletLens(50mm, -40mm, 100mm, -60mm, 8mm, 3mm, 6mm, BMO.inch,
                λ -> 1.5, λ -> 1.7, λ -> 1.6)
            translate3d!(tl, [0.05, 0.05, 0.05])
            xrotate3d!(tl, deg2rad(-60))
            zrotate3d!(tl, deg2rad(45))
            system = System([tl])
            dir = orientation(tl)[:, 2]
            pos = position(tl) - 0.05 * dir
            nv = BMO.normal3d(dir)
            beam = Beam(pos + 3mm * nv, dir, 1e-6)
            solve_system!(system, beam)
            @test length(rays(beam)) == 5
            @test BMO.refractive_index.(rays(beam)) == [1, 1.5, 1.7, 1.6, 1]
        end

        @testset "Composition from individual lenses" begin
            front = SphericalLens(35.86mm, 85.87mm, 11.81mm, 60mm, 1.671)
            middle = SphericalLens(85.87mm, -646.31mm, 7.05mm, 60mm, 1.4892)
            back = Lens(SphericalSurface(-646.31mm, 60mm), SphericalSurface(23.51mm, 40mm),
                1.9mm, λ -> 1.7394)
            translate3d!(middle, [0, thickness(front), 0])
            translate3d!(back, [0, thickness(front) + thickness(middle), 0])
            tl = TripletLens(front, middle, back)
            @test thickness(tl) ≈ 20.76mm
        end
    end

    @testset "Testing triplet lenses ray and beam type compats." begin
        # Every ray/beam type the solver can produce must trace through a TripletLens.
        # Regression: the chief ray of an AstigmaticGaussianBeamlet is a PolarizedRay,
        # which is not a subtype of Ray, so the triplet interact3d must accept AbstractRay.
        n1 = λ -> 1.5
        n2 = λ -> 1.7
        n3 = λ -> 1.6
        λ = 1e-6
        w0 = 1mm
        pos = [0.0, 0, 0]
        dir = [0.0, 1, 0]
        support = [1.0, 0, 0]
        n_seq = [1, n1(λ), n2(λ), n3(λ), 1]

        function triplet_system()
            tl = SphericalTripletLens(50mm, -40mm, 100mm, -60mm, 8mm, 3mm, 6mm, BMO.inch,
                n1, n2, n3)
            translate3d!(tl, [0, 50mm, 0])
            return System([tl])
        end

        @testset "Beam{Ray}" begin
            beam = Beam(pos, dir, λ)
            solve_system!(triplet_system(), beam)
            @test length(rays(beam)) == 5
            @test BMO.refractive_index.(rays(beam)) == n_seq
        end

        @testset "Beam{PolarizedRay}" begin
            beam = Beam(pos, dir, λ, [1.0, 0, 0])
            solve_system!(triplet_system(), beam)
            @test length(rays(beam)) == 5
            @test BMO.refractive_index.(rays(beam)) == n_seq
        end

        @testset "GaussianBeamlet" begin
            gauss = GaussianBeamlet(pos, dir, λ, w0; support)
            solve_system!(triplet_system(), gauss)
            @test length(rays(gauss.chief)) == 5
            @test BMO.refractive_index.(rays(gauss.chief)) == n_seq
        end

        @testset "AstigmaticGaussianBeamlet" begin
            agb = AstigmaticGaussianBeamlet(pos, dir, λ, w0; support)
            solve_system!(triplet_system(), agb)
            @test length(rays(agb.c)) == 5
            @test BMO.refractive_index.(rays(agb.c)) == n_seq
        end
    end
end

end # MODULE