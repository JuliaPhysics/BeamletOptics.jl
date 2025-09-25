module TestRays

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics

const BMO = BeamletOptics

@testset "Rays" begin
    # Testing constructor
    pos = [0, 0, 0]
    dir = [1.0, 1, 0]
    ray = Ray(pos, dir)
    @test ismutable(ray)
    @test isa(ray, Ray{Float64})
    @test isnothing(ray.intersection)
    @test isinf(length(ray))
    @test isapprox(norm(ray.dir), 1)
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
        i1 = BMO.Intersection(nothing, nothing, 0.0, Point3(0, 1.0, 0))
        i2 = BMO.Intersection(nothing, nothing, 0.0, Point3(0, -1.0, 0))
        BMO.intersection!(r1, i1)
        BMO.intersection!(r2, i2)
        @test !BMO.isentering(r1)
        @test BMO.isentering(r2)
        @test !BMO.isentering(r3)
    end

    @testset "Testing refraction3d" begin
        n1 = 1
        n2 = 1.5
        dir = normalize([0, 1, 0])
        ray = Ray(zeros(3), dir)
        nml = normalize(Point3{Float64}(0, -1, 1))
        BMO.intersection!(
            ray, BMO.Intersection(nothing, nothing, 1.0, nml))
        @test BMO.refraction3d(dir, nml, n1, n2) ==
              BMO.refraction3d(ray, n2)
        # test for correct exit normal flip
        nml *= -1
        BMO.intersection!(
            ray, BMO.Intersection(nothing, nothing, 1.0, nml))
        @test BMO.refraction3d(dir, -nml, n1, n2) ==
              BMO.refraction3d(ray, n2)
    end
end

end # MODULE