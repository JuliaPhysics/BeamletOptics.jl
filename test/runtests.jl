using SCDI
using Test
using LinearAlgebra

@testset "Utilites" begin
    @debug "Testing orthogonal3d for right hand rule and unit length"
    orth = SCDI.orthogonal3d([2, 0, 0], [0, 0, 1])
    @test isapprox(orth, [0, -1, 0])

    @debug "Testing rotate3d for clockwise dir. and conservation of length"
    R = SCDI.rotate3d([0, 0, 1], π / 2)
    @test isapprox(R * [1, 0, 0], [0, 1, 0])

    @debug "Testing align3d for rotation and conservation of length"
    target = [0, 0, 1]
    R = SCDI.align3d([1, 0, 0], target)
    @test isapprox(R * [1, 0, 0], target)

    @debug "Testing angle3d for resulting angle"
    a = SCDI.angle3d([1, 0, 0], [0, 0, 1])
    @test isapprox(a, π / 2)

    @debug "Testing fast_dot3d"
    c = SCDI.fast_dot3d([1, 2, 3], [4, 5, 6])
    @test c == 32

    @debug "Testing fast_cross3d"
    c = SCDI.fast_cross3d([1, 0, 0], [0, 1, 0])
    @test isapprox(c, [0, 0, 1])

    @debug "Testing fast_cross3d!"
    c = zeros(Int, 3)
    SCDI.fast_cross3d!(c, [1, 0, 0], [0, 1, 0])
    @test isapprox(c, [0, 0, 1])

    @debug "Testing fast_sub3d!"
    c = zeros(Int, 3)
    SCDI.fast_sub3d!(c, [4, 5, 6], [1, 2, 3])
    @test isapprox(c, [3, 3, 3])

    @debug "Testing line_point_distance3d"
    pos = [0, 0, 0]
    dir = [1, 0, 0]
    point = [5, 1, 1]
    d = SCDI.line_point_distance3d(pos, dir, point)
    @test isapprox(d, √2)
end

@testset "Rays" begin
    @debug "Testing Ray struct definition"
    @test isdefined(SCDI, :Ray)
    pos = [0.0, 0, 0]
    dir = [1.0, 1, 1]
    ray = SCDI.Ray{Float64}(pos, dir)
    @test ismutable(ray)

    @debug "Testing Ray constructor (dir normalization and init. Inf length)"
    @test norm(ray.dir) == 1
    @test isinf(ray.len)

    @debug "Testing Beam struct definition"
    @test isdefined(SCDI, :Beam)
    beam = SCDI.Beam([ray], 1e3)
    @test ismutable(beam)
end

@testset "Mesh" begin
    @testset "Testing type definitions" begin
        @test isdefined(SCDI, :AbstractEntity)
        @test isdefined(SCDI, :AbstractMesh)

        @debug "Testing Mesh struct definition"
        @test isdefined(SCDI, :Mesh)
    end

    include("Cube.jl")

    @debug "Generating stationary and moving test cubes with scale 1"
    foo = Cube(1) # stationary
    bar = Cube(1) # moving

    @testset "Testing translate3d!" begin
        SCDI.translate3d!(bar, -0.5 * [1, 1, 1]) # move COG to origin
        @test minimum(bar.mesh.vertices[:, 1]) == -0.5
        @test minimum(bar.mesh.vertices[:, 2]) == -0.5
        @test minimum(bar.mesh.vertices[:, 3]) == -0.5
        @test maximum(bar.mesh.vertices[:, 1]) == 0.5
        @test maximum(bar.mesh.vertices[:, 2]) == 0.5
        @test maximum(bar.mesh.vertices[:, 3]) == 0.5
        @test all(bar.mesh.pos .== -0.5)
    end

    @testset "Testing set_new_origin3d!" begin
        SCDI.set_new_origin3d!(bar)
        @test all(iszero.(bar.mesh.pos))
    end

    @testset "Testing x/y/zrotate3d!" begin
        @debug "Testing xrotate3d!"
        SCDI.xrotate3d!(bar, π / 4)
        @test isapprox(minimum(bar.mesh.vertices[:, 1]), -0.5)
        @test isapprox(minimum(bar.mesh.vertices[:, 2]), -√2 / 2)
        @test isapprox(minimum(bar.mesh.vertices[:, 3]), -√2 / 2)
        @test isapprox(maximum(bar.mesh.vertices[:, 1]), 0.5)
        @test isapprox(maximum(bar.mesh.vertices[:, 2]), √2 / 2)
        @test isapprox(maximum(bar.mesh.vertices[:, 3]), √2 / 2)

        @debug "Testing yrotate3d!"
        SCDI.yrotate3d!(bar, π / 2)
        @test isapprox(minimum(bar.mesh.vertices[:, 1]), -√2 / 2)
        @test isapprox(minimum(bar.mesh.vertices[:, 2]), -√2 / 2)
        @test isapprox(minimum(bar.mesh.vertices[:, 3]), -0.5)
        @test isapprox(maximum(bar.mesh.vertices[:, 1]), √2 / 2)
        @test isapprox(maximum(bar.mesh.vertices[:, 2]), √2 / 2)
        @test isapprox(maximum(bar.mesh.vertices[:, 3]), 0.5)

        @debug "Testing zrotate3d!"
        SCDI.zrotate3d!(bar, π / 4)
        @test isapprox(minimum(bar.mesh.vertices[:, 1]), -0.5)
        @test isapprox(minimum(bar.mesh.vertices[:, 2]), -0.5)
        @test isapprox(minimum(bar.mesh.vertices[:, 3]), -0.5)
        @test isapprox(maximum(bar.mesh.vertices[:, 1]), 0.5)
        @test isapprox(maximum(bar.mesh.vertices[:, 2]), 0.5)
        @test isapprox(maximum(bar.mesh.vertices[:, 3]), 0.5)

        @debug "Testing orientation of dir matrix"
        @test all(bar.mesh.dir[[2, 6, 7]] .== 1)
    end

    @testset "Testing reset_rotation3d!" begin
        SCDI.reset_rotation3d!(bar)
        SCDI.reset_translation3d!(bar) # sneaky, useless command for that sweet code coverage
        @test all(bar.mesh.dir .== Matrix{Float64}(I, 3, 3))
    end

    @testset "Testing orthogonal3d" begin
        normal = SCDI.orthogonal3d(bar, 1)
        @test isapprox(normal, [0, 0, -1])
    end

    @testset "Testing scale3d!" begin
        SCDI.scale3d!(bar, 2)
        @test isapprox(minimum(bar.mesh.vertices[:, 1]), -1)
        @test isapprox(minimum(bar.mesh.vertices[:, 2]), -1)
        @test isapprox(minimum(bar.mesh.vertices[:, 3]), -1)
        @test isapprox(maximum(bar.mesh.vertices[:, 1]), 1)
        @test isapprox(maximum(bar.mesh.vertices[:, 2]), 1)
        @test isapprox(maximum(bar.mesh.vertices[:, 3]), 1)
        @test bar.mesh.scale == 2
    end
end
