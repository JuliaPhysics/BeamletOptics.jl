module TestSDFs

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics

const BMO = BeamletOptics

@testset "SDFs" begin
    @testset "Testing type definitions" begin
        @test isdefined(BMO, :AbstractSDF)
        @test isdefined(BMO, :SphereSDF)
        @test isdefined(BMO, :CylinderSDF)
        @test isdefined(BMO, :CutSphereSDF)
        @test isdefined(BMO, :ThinLensSDF)
    end

    # Orientation-less test point sdf
    mutable struct TestPointSDF{T} <: BMO.AbstractSDF{T}
        position::Point3{T}
        orientation::Matrix{T}
    end

    TestPointSDF(p::AbstractArray{T}) where {T} = TestPointSDF{T}(
        Point3{T}(p), Matrix{T}(I, 3, 3))
    TestPointSDF(T = Float64) = TestPointSDF{T}(Point3{T}(0), Matrix{T}(I, 3, 3))

    BMO.position(tps::TestPointSDF) = tps.position
    BMO.position!(tps::TestPointSDF{T}, new::Point3{T}) where {T} = (tps.position = new)

    BMO.orientation(tps::TestPointSDF) = tps.orientation
    BMO.orientation!(tps::TestPointSDF{T}, new::Matrix{T}) where {T} = (tps.orientation = new)

    BMO.transposed_orientation(tps::TestPointSDF) = transpose(tps.orientation)
    BMO.transposed_orientation!(::TestPointSDF, ::Any) = nothing

    function BMO.sdf(tps::TestPointSDF, point)
        p = BMO._world_to_sdf(tps, point)
        return norm(p)
    end

    @testset "Testing kinematics and transforms" begin
        point = TestPointSDF()
        t = 10
        θ = deg2rad(30)
        translate3d!(point, [t, 0, 0])
        rotate3d!(point, [0, 1, 0], θ)
        pt = BMO._world_to_sdf(point, [0, 0, 0])
        @test pt[1] ≈ -t * cos(θ)
        @test pt[2] ≈ 0
        @test pt[3] ≈ -t * sin(θ)
    end

    @testset "Testing intersect3d" begin
        t = 10.0
        point = TestPointSDF(zeros(3))
        translate3d!(point, [t, 0, 0])

        r1 = Ray(zeros(3), [1.0, 0, 0])
        r2 = Ray(zeros(3), [1.0, 1, 0])
        r3 = Ray(zeros(3), [1.0, 0, 1])

        i1 = BMO.intersect3d(point, r1)
        i2 = BMO.intersect3d(point, r2)
        i3 = BMO.intersect3d(point, r3)

        @test length(i1) == t
        @test isnothing(i2)
        @test isnothing(i3)
    end

    @testset "Testing normal3d" begin
        point = TestPointSDF(zeros(3))
        offset = [5, 0, 0]
        translate3d!(point, offset)
        p1 = [1, 0, 0]
        p2 = [0, 1, 0]
        p3 = [0, 0, 1]
        @test BMO.normal3d(point, p1 + offset) == p1
        @test BMO.normal3d(point, p2 + offset) == p2
        @test BMO.normal3d(point, p3 + offset) == p3
    end
end

end # MODULE