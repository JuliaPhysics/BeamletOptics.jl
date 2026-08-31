module TestBeams

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics

const BMO = BeamletOptics

@testset "Beams" begin
    is = Intersection(1.0, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    r1 = Ray(Ray([0.0, 0, 0], [1, 0, 0]), LineSegment([0.0, 0, 0], [1, 0, 0], is, 1.0))
    r2 = Ray(Ray([1.0, 0, 0], [0, 1, 0]), LineSegment([1.0, 0, 0], [0, 1, 0], is, 1.0))
    r3 = Ray(Ray([1.0, 1, 0], [0, 0, 1]), LineSegment([1.0, 1, 0], [0, 0, 1], is, 1.0))
    r4 = Ray([1.0, 1, 1], [1, 0, 0]) # OpenSegment

    beam = Beam(r1)
    empty!(beam.rays)
    push!(beam, r1)
    push!(beam, r2)
    push!(beam, r3)
    push!(beam, r4)
    @test length(beam) == 3
    @test length(beam.rays) == 4
    @test length(BMO.rays(beam)) == 4
    @test BMO.point_on_beam(beam, 0) == ([0, 0, 0], 1)
    @test BMO.point_on_beam(beam, 1) == ([1, 0, 0], 2)
    @test BMO.point_on_beam(beam, 2) == ([1, 1, 0], 3)
    @test BMO.point_on_beam(beam, 3) == ([1, 1, 1], 4)
    @test BMO.isparentbeam(beam, r1) == true
    @test BMO.isparentbeam(beam, beam.rays[2]) == true
end

end # MODULE