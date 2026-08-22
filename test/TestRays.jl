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
    @test isa(ray, Ray{Float64})
    @test isapprox(norm(ray.dir), 1)
    @test BMO.position(ray) == pos
    @test BMO.direction(ray) ≈ normalize(dir)
    @test BMO.wavelength(ray) == 1000e-9
    @test BMO.refractive_index(ray) == 1.0
    # Test RaySegment length
    seg_unhit = RaySegment(ray, nothing)
    @test isinf(length(seg_unhit))
    seg_hit = RaySegment(ray, dummy_intersection(5.0, [0, 1, 0]))
    @test length(seg_hit) == 5.0
    @test BMO.optical_path_length(seg_hit) == 5.0
    # Test helper functions
    @test BMO.line_point_distance3d(ray, [1, 1, 0]) == 0
    @test BMO.line_point_distance3d(ray, [-1, 1, 0]) == sqrt(2)
    # Test error msg
    @test_throws ErrorException Ray(zeros(3), zeros(3))
    @test_throws ErrorException Ray(zeros(3), ones(3)*eps())

    @testset "Testing isentering" begin
        r1 = Ray([0, 0, 0], [0, 1, 0])
        r2 = Ray([0, 0, 0], [0, 1, 0])
        r3 = Ray([0, 0, 0], [0, 1, 0])
        i1 = dummy_intersection(0.0, [0, 1.0, 0])
        i2 = dummy_intersection(0.0, [0, -1.0, 0])
        seg1 = RaySegment(r1, i1)
        seg2 = RaySegment(r2, i2)
        seg3 = RaySegment(r3, nothing)
        @test !BMO.isentering(seg1)
        @test BMO.isentering(seg2)
        @test !BMO.isentering(seg3)
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
        seg1 = RaySegment(ray, int1)
        @test BMO.refraction3d(dir, nml, n1, n2) ==
              BMO.refraction3d(seg1, n2)
        # test for correct exit normal flip
        nml *= -1
        int2 = dummy_intersection(1.0, nml)
        seg2 = RaySegment(ray, int2)
        @test BMO.refraction3d(dir, -nml, n1, n2) ==
              BMO.refraction3d(seg2, n2)
    end
end

end # MODULE