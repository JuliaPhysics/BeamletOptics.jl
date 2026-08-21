module TestBeams

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics

const BMO = BeamletOptics

@testset "Beams" begin
    is = Intersection(1.0, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    r1 = Ray([0.0, 0, 0], [1, 0, 0])
    r2 = Ray([1.0, 0, 0], [0, 1, 0])
    r3 = Ray([1.0, 1, 0], [0, 0, 1])
    r4 = Ray([1.0, 1, 1], [1, 0, 0])
    seg1 = RaySegment(r1, is)
    seg2 = RaySegment(r2, is)
    seg3 = RaySegment(r3, is)
    seg4 = RaySegment(r4, nothing)
    # Test beam
    beam = Beam(r1)
    push!(beam, seg1)
    push!(beam, seg2)
    push!(beam, seg3)
    push!(beam, seg4)
    @test length(beam) == 3
    @test length(beam.segments) == 4
    @test length(BMO.rays(beam)) == 4
    @test BMO.point_on_beam(beam, 0) == ([0, 0, 0], 1)
    @test BMO.point_on_beam(beam, 1) == ([1, 0, 0], 2)
    @test BMO.point_on_beam(beam, 2) == ([1, 1, 0], 3)
    @test BMO.point_on_beam(beam, 3) == ([1, 1, 1], 4)
    @test BMO.point_on_beam(beam, 10) == ([8, 1, 1], 4)
    @test BMO.isparentbeam(beam, r2) == true
end

end # MODULE