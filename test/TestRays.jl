module TestRays

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics

const BMO = BeamletOptics

dummy_intersection(t, n) = Intersection(t, [0.0, 0.0, 0.0], Point3(n))

@testset "Rays" begin
    # Testing constructor
    pos = [0, 0, 0]
    dir = [1.0, 1, 0]
    ray = Ray(pos, dir)
    @test !ismutable(ray)
    @test isa(ray, Ray{Float64, OpenSegment{Float64}})
    @test isapprox(norm(BMO.direction(ray)), 1)
    @test BMO.position(ray) == pos
    @test BMO.direction(ray) ≈ normalize(dir)
    @test BMO.wavelength(ray) == 1000e-9
    @test BMO.refractive_index(ray) == 1.0
    # Test Segment length
    @test isinf(length(ray))
    hit_int = dummy_intersection(5.0, [0, 1, 0])
    ray_hit = Ray(ray, LineSegment(pos, dir, hit_int, 5.0))
    @test length(ray_hit) == 5.0
    @test BMO.optical_path_length(ray_hit) == 5.0
    # Test helper functions
    @test BMO.line_point_distance3d(ray, [1, 1, 0]) == 0
    @test BMO.line_point_distance3d(ray, [-1, 1, 0]) == sqrt(2)
    # Test error msg
    @test_throws ErrorException Ray(zeros(3), zeros(3))
    @test_throws ErrorException Ray(zeros(3), ones(3)*eps())

    @testset "Testing isentering" begin
        i1 = dummy_intersection(0.0, [0, 1.0, 0])
        i2 = dummy_intersection(0.0, [0, -1.0, 0])
        r1 = Ray(Ray([0, 0, 0], [0, 1, 0]), LineSegment([0, 0, 0], [0, 1, 0], i1, 0.0))
        r2 = Ray(Ray([0, 0, 0], [0, 1, 0]), LineSegment([0, 0, 0], [0, 1, 0], i2, 0.0))
        r3 = Ray([0, 0, 0], [0, 1, 0])
        @test !BMO.isentering(r1)
        @test BMO.isentering(r2)
        @test !BMO.isentering(r3)
        @test !BMO.isentering(r1, i1)
        @test BMO.isentering(r2, i2)
    end

    @testset "Testing refraction3d" begin
        n1 = 1
        n2 = 1.5
        dir = normalize([0, 1, 0])
        ray = Ray(zeros(3), dir)
        nml = normalize(Point3{Float64}(0, -1, 1))
        int1 = dummy_intersection(1.0, nml)
        ray1 = Ray(ray, LineSegment(zeros(3), dir, int1, 1.0))
        @test BMO.refraction3d(dir, nml, n1, n2) ==
              BMO.refraction3d(ray1, n2)
        # test for correct exit normal flip
        nml *= -1
        int2 = dummy_intersection(1.0, nml)
        ray2 = Ray(ray, LineSegment(zeros(3), dir, int2, 1.0))
        @test BMO.refraction3d(dir, -nml, n1, n2) ==
              BMO.refraction3d(ray2, n2)
    end
end

end # MODULE