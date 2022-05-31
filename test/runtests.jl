using SCDI
using Test
using LinearAlgebra

@testset "Ray tracing utilites" begin
    @test SCDI.orthogonal3d([2,0,0], [0,0,1]) == [0,-1,0]
    @test norm(SCDI.rotate3d([0,0,1], π/2)*[1,0,0] - [0,1,0]) < 1e-9
    @test SCDI.align3d([1,0,0],[1,1,0])*[sqrt(2),0,0] == [1,1,0]
end
