using SCDI
using Test
using LinearAlgebra

@testset "Ray tracing utilites" begin
    @test SCDI.orthogonal3d([2, 0, 0], [0, 0, 1]) == [0, -1, 0]
    @test isapprox(SCDI.rotate3d([0, 0, 1], π / 2) * [1, 0, 0], [0, 1, 0], atol=1e-9)
    @test SCDI.align3d([1, 0, 0], [1, 1, 0]) * [sqrt(2), 0, 0] == [1, 1, 0]
end

@testset "Mesh/geometry backend" begin
    struct Plane{T} <: SCDI.AbstractMesh
        mesh::SCDI.Mesh{T}
    end
    #@test isdefined(SCDI, :Plane)

    vertices = [
        1 1 0
        1 -1 0
        -1 -1 0
        -1 1 0
    ]
    faces = [
        1 2 3
        3 4 1
    ]
    pos = [0, 0, 0]
    dir = Matrix{Int}(I, 3, 3)
    scale = 1
    plane = Plane{Float64}(SCDI.Mesh{Float64}(
        vertices,
        faces,
        dir,
        pos,
        scale
    ))

    SCDI.translate3d!(plane, [0, 0, 1])
    @test all(plane.mesh.pos .== pos + [0, 0, 1])

    SCDI.xrotate3d!(plane, π / 2)
    @test all(plane.mesh.vertices[2:3, 3] .== 2)

    SCDI.yrotate3d!(plane, π / 4)
    @test isapprox(plane.mesh.vertices[2, 3], 1 + sqrt(2), atol=5e-6)

    SCDI.zrotate3d!(plane, π / 2)
    @test isapprox(plane.mesh.vertices[3, 2], sqrt(2), atol=5e-6)

    SCDI.zrotate3d!(plane, -π / 2)
    SCDI.yrotate3d!(plane, -π / 4)
    SCDI.xrotate3d!(plane, -π / 2)
    SCDI.translate3d!(plane, [0, 0, -1])
    @test isapprox(plane.mesh.vertices, vertices, atol=1e-8)

    SCDI.scale3d!(plane, 2)
    @test isapprox(plane.mesh.vertices, vertices .* 2, atol=1e-8)
    @test plane.mesh.scale == 2
end
