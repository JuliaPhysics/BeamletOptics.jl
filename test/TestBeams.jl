module TestBeams

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics

const BMO = BeamletOptics

mutable struct TestShapeless{T} <: BMO.AbstractShape{T}
    pos::Vector{T}
    dir::Matrix{T}
end
TestShapeless() = TestShapeless{Float64}(zeros(3), Matrix{Float64}(I, 3, 3))

struct TestObject{T, S <: BMO.AbstractShape{T}} <: BMO.AbstractObject{T}
    shape::S
end
TestObject() = TestObject(TestShapeless())

@testset "Beams" begin
    is = BMO.ObjectIntersection(TestObject(), BMO.ShapeIntersection(TestShapeless(), 1.0, Point3(zeros(3))))
    r1 = Ray([0.0, 0, 0], [1, 0, 0])
    r2 = Ray([1.0, 0, 0], [0, 1, 0])
    r3 = Ray([1.0, 1, 0], [0, 0, 1])
    r4 = Ray([1.0, 1, 1], [1, 0, 0])
    BMO.intersection!(r1, is)
    BMO.intersection!(r2, is)
    BMO.intersection!(r3, is)
    BMO.intersection!(r4, nothing)
    # Test beam
    beam = Beam(r1)
    push!(beam, r2)
    push!(beam, r3)
    push!(beam, r4)
    @test length(beam) == 3
    @test BMO.point_on_beam(beam, 0) == ([0, 0, 0], 1)
    @test BMO.point_on_beam(beam, 1) == ([1, 0, 0], 2)
    @test BMO.point_on_beam(beam, 2) == ([1, 1, 0], 3)
    @test BMO.point_on_beam(beam, 3) == ([1, 1, 1], 4)
    @test BMO.point_on_beam(beam, 10) == ([8, 1, 1], 4)
    @test BMO.isparentbeam(beam, r2) == true
end

end # MODULE