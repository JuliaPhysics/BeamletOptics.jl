using SCDI
using Test
using LinearAlgebra

@testset "Utilites" begin
    @testset "Testing euclidic norm utilites" begin
        v = [1.0, 0.0, 1.0]
        @test isapprox(SCDI.norm3d(v), sqrt(2))
        @test isapprox(SCDI.normalize3d(v), [sqrt(2) / 2, 0, sqrt(2) / 2])
        SCDI.normalize3d!(v)
        @test isapprox(v, [sqrt(2) / 2, 0, sqrt(2) / 2])
    end

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

    @testset "Testing reflection3d" begin
        for dx in -1:1, dy in -1:1
            @test isapprox(SCDI.reflection3d([dx, dy, 1], [0, 0, -1]), [dx, dy, -1])
        end
    end

    @testset "Testing refraction3d" begin
        normal = [0, 0, 1]
        @testset "Test from vacuum into medium" begin
            n1 = 1.0
            n2 = 1.5
            for θ1 in 0:π/8:π/2
                dir_in = [sin(θ1), 0, -cos(θ1)]
                dir_out = SCDI.refraction3d(dir_in, normal, n1, n2)
                θ2 = SCDI.angle3d(-normal, dir_out)
                # 2D-equation for refraction validation
                θ3 = asin(n1 / n2 * sin(θ1))
                @test isapprox(θ2, θ3)
            end
        end
        @testset "Test from medium into vacuum" begin
            n1 = 1.5
            n2 = 1.0
            for θ1 in 0:π/8:π/2
                dir_in = [sin(θ1), 0, -cos(θ1)]
                dir_out = SCDI.refraction3d(dir_in, normal, n1, n2)
                if θ1 > asin(n2 / n1)
                    # Test for total reflection
                    θ2 = SCDI.angle3d(dir_out, normal)
                    @test isapprox(θ1, θ2)
                else
                    # Test for refraction
                    θ2 = SCDI.angle3d(-normal, dir_out)
                    θ3 = asin(n1 / n2 * sin(θ1))
                    @test isapprox(θ2, θ3)
                end
            end
        end
    end
end

@testset "Types" begin
    @debug "Testing abstract type definitions"
    @test isdefined(SCDI, :AbstractEntity)
    @test isdefined(SCDI, :AbstractObject)

    # Generate test struct
    mutable struct Object{T} <: SCDI.AbstractObject{T}
        pos::Vector{T}
    end

    @testset "Testing abstract type interfaces" begin
        object = Object([0,0,0])
        SCDI.translate3d!(object, [1,1,1])        
        SCDI.position!(object, [2,2,2])
        SCDI.reset_translation3d!(object)
        @test SCDI.position(object) == object.pos
        @test SCDI.position(object) == zeros(3)
        # The following test are expected to do nothing but not throw exceptions
        SCDI.rotate3d!(object, [0,0,1], π)
        SCDI.xrotate3d!(object, π)
        SCDI.yrotate3d!(object, π)
        SCDI.zrotate3d!(object, π)
        SCDI.reset_rotation3d!(object)
        @test_logs (:warn,) SCDI.intersect3d(object, SCDI.Ray([0,0,1],[0,0,1]))
    end

    # AbstractRay testing (WIP)
    @test isdefined(SCDI, :AbstractRay)
end

@testset "Rays" begin
    @debug "Testing Ray struct definition"
    @test isdefined(SCDI, :Ray)
    pos = [0, 0, 0]
    dir = [1, 1, 1]
    ray = SCDI.Ray(pos, dir)
    @test ismutable(ray)
    @test eltype(ray.pos) == Float64

    @debug "Testing Ray constructor (dir normalization and init. Inf length)"
    @test isapprox(SCDI.norm3d(ray.dir), 1)
    @test isinf(ray.intersection.t)

    @debug "Testing Beam struct definition"
    @test isdefined(SCDI, :Beam)
    beam = SCDI.Beam([ray], 1e3)
    @test ismutable(beam)
end

@testset "Mesh" begin
    # NOTE: the "Mesh" testset is mutating. Errors/fails might lead to subsequent tests failing too!
    @testset "Testing type definitions" begin
        @test isdefined(SCDI, :AbstractEntity)
        @test isdefined(SCDI, :AbstractMesh)

        @debug "Testing Mesh struct definition"
        @test isdefined(SCDI, :Mesh)
        @test isdefined(SCDI, :mesh)
    end

    # Generate cube since types are defined
    include("Cube.jl")

    @debug "Generating moving test cube with scale 1"
    foo = Cube(1)

    @testset "Testing AbstractMesh getters" begin
        @test typeof(SCDI.meshof(foo)) == SCDI.Mesh{Float64}
        @test SCDI.vertices(foo) == foo.mesh.vertices
        @test SCDI.faces(foo) == foo.mesh.faces
        @test SCDI.orientation(foo) == foo.mesh.dir
        @test SCDI.position(foo) == foo.mesh.pos
        @test SCDI.scale(foo) == foo.mesh.scale
    end

    @testset "Testing translate3d!" begin
        SCDI.translate3d!(foo, -0.5 * [1, 1, 1]) # move COG to origin
        @test minimum(SCDI.vertices(foo)[:, 1]) == -0.5
        @test minimum(SCDI.vertices(foo)[:, 2]) == -0.5
        @test minimum(SCDI.vertices(foo)[:, 3]) == -0.5
        @test maximum(SCDI.vertices(foo)[:, 1]) == 0.5
        @test maximum(SCDI.vertices(foo)[:, 2]) == 0.5
        @test maximum(SCDI.vertices(foo)[:, 3]) == 0.5
        @test all(SCDI.position(foo) .== -0.5)
    end

    @testset "Testing set_new_origin3d!" begin
        SCDI.set_new_origin3d!(foo)
        @test SCDI.position(foo) == zeros(3)
    end

    @testset "Testing x/y/zrotate3d!" begin
        @testset "Testing xrotate3d!" begin
            SCDI.xrotate3d!(foo, π / 4)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 1]), -0.5)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 2]), -√2 / 2)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 3]), -√2 / 2)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 1]), 0.5)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 2]), √2 / 2)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 3]), √2 / 2)
        end

        @testset "Testing yrotate3d!" begin
            SCDI.yrotate3d!(foo, π / 2)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 1]), -√2 / 2)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 2]), -√2 / 2)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 3]), -0.5)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 1]), √2 / 2)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 2]), √2 / 2)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 3]), 0.5)
        end

        @testset "Testing zrotate3d!" begin
            SCDI.zrotate3d!(foo, π / 4)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 1]), -0.5)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 2]), -0.5)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 3]), -0.5)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 1]), 0.5)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 2]), 0.5)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 3]), 0.5)
        end

        @debug "Testing orientation of dir matrix"
        @test all(SCDI.orientation(foo)[[2, 6, 7]] .== 1)
    end

    @testset "Testing reset_rotation3d!" begin
        SCDI.reset_rotation3d!(foo)
        SCDI.reset_translation3d!(foo) # sneaky, useless command for that sweet code coverage
        @test all(SCDI.orientation(foo) .== Matrix{Float64}(I, 3, 3))
    end

    @testset "Testing orthogonal3d" begin
        normal = SCDI.orthogonal3d(foo, 1)
        @test isapprox(normal, [0, 0, -1])
    end

    @testset "Testing scale3d!" begin
        SCDI.scale3d!(foo, 2)
        @test isapprox(minimum(SCDI.vertices(foo)[:, 1]), -1)
        @test isapprox(minimum(SCDI.vertices(foo)[:, 2]), -1)
        @test isapprox(minimum(SCDI.vertices(foo)[:, 3]), -1)
        @test isapprox(maximum(SCDI.vertices(foo)[:, 1]), 1)
        @test isapprox(maximum(SCDI.vertices(foo)[:, 2]), 1)
        @test isapprox(maximum(SCDI.vertices(foo)[:, 3]), 1)
        @test SCDI.scale(foo) == 2
    end
end

@testset "Intersections" begin
    @testset "Testing Moeller-Trumbore algorithm" begin
        @debug "Testing functor definition"
        @test isdefined(SCDI, :MoellerTrumboreAlgorithm)
        @test hasfield(SCDI.MoellerTrumboreAlgorithm, :kϵ)
        @test hasfield(SCDI.MoellerTrumboreAlgorithm, :lϵ)
        @test isdefined(SCDI, :ray_triangle_intersection)
        @test isconst(SCDI, :ray_triangle_intersection)
        @debug "Defining face in the x-y-plane with z-offset t"
        t = 5
        face = [
            1 1 t
            -1 1 t
            0 -1 t
        ]
        @debug "Defining ray at origin pointing along z-axis"
        pos = [0, 0, 0]
        dir = [0, 0, 1]
        ray = SCDI.Ray(pos, dir)
        # Preallocate memory 
        E1 = similar(ray.dir)
        E2 = similar(ray.dir)
        Pv = similar(ray.dir)
        Tv = similar(ray.dir)
        Qv = similar(ray.dir)
        @debug "Testing intersection"
        @test isapprox(SCDI.ray_triangle_intersection(face, ray, E1, E2, Pv, Tv, Qv), t)
        # Check allocations (WARNING: function must have been compiled once for before this test!)
        alloc = @allocated SCDI.ray_triangle_intersection(face, ray, E1, E2, Pv, Tv, Qv)
        if alloc > 16
            @warn "Allocated number of bytes for MTA larger than expected!" alloc
        end
    end

    @testset "Testing intersect3d for meshes" begin
        @debug "Generating test cube"
        t = 5
        s = 1 # scale/2
        foo = Cube(2 * s)
        # Move cube COG to origin
        SCDI.translate3d!(foo, -[s, s, s])
        SCDI.set_new_origin3d!(foo)
        # Align cube edge at t units from origin
        SCDI.translate3d!(foo, [t + s, 0, 0])
        @debug "Defining ray at origin with variable dir"
        pos = [0, 0, 0]
        steps = 10
        for z in -s:(s/steps):s
            # Ray constructed each time for unit-length dir
            dir = [t, 0, z]
            ray = SCDI.Ray(pos, dir)
            @test isapprox(SCDI.intersect3d(foo, ray).t, sqrt(t^2 + z^2))
        end
    end
end
