module TestBeams

using BeamletOptics
using Test

const BMO = BeamletOptics

@testset "Beams" begin
    is = BMO.Intersection(1.0, zeros(3))
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

@testset "Beam of PolarizedRays with a complex field" begin
    # a complex `E0`, e.g. circular polarization, gives a beam of real number type
    for E0 in ([1.0, 0, 0], [1, 0, im] / sqrt(2), ComplexF32[1, 0, 0], [1, 0, 0])
        b = Beam([0.0, 0, 0], [0.0, 1, 0], 1e-6, E0)
        @test b isa Beam{Float64, PolarizedRay{Float64}}
        @test BMO.polarization(first(BMO.rays(b))) ≈ E0
    end
end

end # MODULE