module TestBeamGroups

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics

const BMO = BeamletOptics

@testset "Beam groups" begin
    @testset "AbstractBeamGroup definitions" begin
        @test isdefined(BMO, :AbstractBeamGroup)

        struct BeamTestGroup{T} <: BMO.AbstractBeamGroup{T, Ray{T}}
            central_beam::BMO.Beam{T, Ray{T}}
        end

        struct BeamTestSystem <: BMO.AbstractSystem end
        BMO.objects(::BeamTestSystem) = [nothing]
        BMO.intersect3d(::Nothing, ::BMO.AbstractRay) = nothing
        BMO.interact3d(::Nothing) = nothing

        BMO.beams(tg::BeamTestGroup) = [tg.central_beam]

        pos = [0,0,0]
        dir = [0,1,0]
        lambda = 1064e-9

        tg = BeamTestGroup(BMO.Beam(pos, dir, lambda))

        @test position(tg) == pos
        @test BMO.direction(tg) == dir
        @test BMO.wavelength(tg) == lambda

        ts = BeamTestSystem()

        @test isnothing(solve_system!(ts, tg))
    end

    is_strictly_sorted(v) = all(diff(v) .> 0)

    @testset "Point source" begin
        # define parameters
        lambda = 486.0e-9
        pos = [0, -0.5, 0]
        dir = [0, 1, 1]
        alpha = deg2rad(2)
        NA = BMO.numerical_aperture(alpha)
        num_rays = 1000
        num_rings = 10

        source = PointSource(pos, dir, alpha, lambda; num_rays, num_rings)

        @testset "Testing point source getters" begin
            @test BMO.numerical_aperture(source) == NA
            @test position(source) == pos
            @test BMO.direction(source) ≈ normalize(dir)
        end

        @testset "Testing point source max. spread" begin
            last_ray = first(BMO.rays(last(BMO.beams(source))))
            @test BMO.angle3d(dir, BMO.direction(last_ray)) ≈ alpha atol = 1e-14
            @test position(last_ray) == pos
            @test length(BMO.beams(source)) == num_rays
        end

        @testset "Testing generated angles" begin
            # Test for center ray 0 deg, generated spread angles etc.
            directions = BMO.direction.(first.(BMO.rays.(BMO.beams(source))))
            angles = BMO.angle3d.(Ref(dir), directions)
            generated_angles = unique(round.(angles, digits=11))
            required_angles = LinRange(0, alpha, num_rings)
            @test generated_angles ≈ required_angles
            # Test for only one center ray, strictly increasing rays per ring
            rays_per_angle = zeros(Int, length(required_angles))
            for (i, _required) in enumerate(required_angles)
                rays_per_angle[i] = count(_required .≈ angles)
            end
            @test first(rays_per_angle) == 1
            @test is_strictly_sorted(rays_per_angle)
            # Test if all generated rays are unique
            @test length(unique(directions)) == length(directions)
        end

        @testset "Testing throw errors" begin
            @test_throws ErrorException PointSource(pos, dir, 1.1*π, lambda; num_rays, num_rings)
            @test_throws ErrorException PointSource(pos, dir, alpha, lambda; num_rays=100, num_rings=10)
        end
    end

    @testset "Collimated source" begin
        # define parameters
        lambda = 486.0e-9
        pos = [0, -0.5, 0]
        dir = [0, 1, 0]
        diameter = 2BMO.inch
        num_rays = 500
        num_rings = 5

        source = CollimatedSource(pos, dir, diameter; num_rays, num_rings)

        @testset "Testing coll. source getters" begin
            @test BMO.diameter(source) == diameter
            @test position(source) == pos
            @test BMO.direction(source) == dir
            @test BMO.wavelength(source) == 1e-6
        end

        @testset "Testing coll. source max. spread diameter" begin
            last_ray = first(BMO.rays(last(BMO.beams(source))))
            @test BMO.direction(last_ray) == dir
            @test norm(position(last_ray) - pos) ≈ diameter/2
            @test length(BMO.beams(source)) == num_rays
        end

        @testset "Testing coll. source generated positions" begin
            # Test for center ray 0 offset, generated radii
            positions = position.(first.(BMO.rays.(BMO.beams(source))))
            radii = norm.(positions .- Ref(pos))
            generated_pos = unique(round.(radii, digits=11))
            required_pos = LinRange(0, diameter/2, num_rings)
            @test generated_pos ≈ required_pos
            # Test for only one center ray, strictly increasing rays per ring
            rays_per_radius = zeros(Int, length(required_pos))
            for (i, _required) in enumerate(required_pos)
                rays_per_radius[i] = count(_required .≈ radii)
            end
            @test first(rays_per_radius) == 1
            @test is_strictly_sorted(rays_per_radius)
            # Test if all generated rays are unique
            @test length(unique(positions)) == length(positions)
        end

        @testset "Testing throw errors" begin
            @test_throws ErrorException CollimatedSource(pos, dir, diameter; num_rays=100, num_rings=10)
        end
    end
end

end # MODULE