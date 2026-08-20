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