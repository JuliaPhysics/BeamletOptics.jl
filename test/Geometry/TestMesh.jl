module TestMesh

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics

@testset "Mesh" begin
    # NOTE: the "Mesh" testset is mutating. Errors/fails might lead to subsequent tests failing too!
    @test isdefined(BMO, :AbstractMesh)
    @test isdefined(BMO, :Mesh)

    # Generate cube since types are defined
    foo = BMO.CubeMesh(1) # test cube
    bar = BMO.CubeMesh(1) # reference cube

    to_origin = -0.5 * [1, 1, 1]

    @testset "Testing AbstractMesh getters" begin
        @test typeof(foo) == BMO.Mesh{Float64}
        @test BMO.vertices(foo) == foo.vertices
        @test BMO.faces(foo) == foo.faces
        @test orientation(foo) == foo.dir
        @test position(foo) == foo.pos
        @test BMO.scale(foo) == foo.scale
    end

    @testset "Testing translate3d!" begin
        translate3d!(foo, to_origin) # move COG to origin
        @test minimum(BMO.vertices(foo)[:, 1]) == -0.5
        @test minimum(BMO.vertices(foo)[:, 2]) == -0.5
        @test minimum(BMO.vertices(foo)[:, 3]) == -0.5
        @test maximum(BMO.vertices(foo)[:, 1]) == 0.5
        @test maximum(BMO.vertices(foo)[:, 2]) == 0.5
        @test maximum(BMO.vertices(foo)[:, 3]) == 0.5
        @test all(position(foo) .== -0.5)
    end

    @testset "Testing set_new_origin3d!" begin
        BMO.set_new_origin3d!(foo)
        @test position(foo) == zeros(3)
    end

    @testset "Testing x/y/zrotate3d!" begin
        @testset "Testing rotate3d!" begin
            rotate3d!(foo, [1, 0, 0], π / 4)
            @test isapprox(minimum(BMO.vertices(foo)[:, 1]), -0.5)
            @test isapprox(minimum(BMO.vertices(foo)[:, 2]), -√2 / 2)
            @test isapprox(minimum(BMO.vertices(foo)[:, 3]), -√2 / 2)
            @test isapprox(maximum(BMO.vertices(foo)[:, 1]), 0.5)
            @test isapprox(maximum(BMO.vertices(foo)[:, 2]), √2 / 2)
            @test isapprox(maximum(BMO.vertices(foo)[:, 3]), √2 / 2)
            # Return to original rotation
            rotate3d!(foo, [1, 0, 0], -π / 4)
        end

        @testset "Testing xrotate3d!" begin
            xrotate3d!(foo, π / 4)
            @test isapprox(minimum(BMO.vertices(foo)[:, 1]), -0.5)
            @test isapprox(minimum(BMO.vertices(foo)[:, 2]), -√2 / 2)
            @test isapprox(minimum(BMO.vertices(foo)[:, 3]), -√2 / 2)
            @test isapprox(maximum(BMO.vertices(foo)[:, 1]), 0.5)
            @test isapprox(maximum(BMO.vertices(foo)[:, 2]), √2 / 2)
            @test isapprox(maximum(BMO.vertices(foo)[:, 3]), √2 / 2)
        end

        @testset "Testing yrotate3d!" begin
            yrotate3d!(foo, π / 2)
            @test isapprox(minimum(BMO.vertices(foo)[:, 1]), -√2 / 2)
            @test isapprox(minimum(BMO.vertices(foo)[:, 2]), -√2 / 2)
            @test isapprox(minimum(BMO.vertices(foo)[:, 3]), -0.5)
            @test isapprox(maximum(BMO.vertices(foo)[:, 1]), √2 / 2)
            @test isapprox(maximum(BMO.vertices(foo)[:, 2]), √2 / 2)
            @test isapprox(maximum(BMO.vertices(foo)[:, 3]), 0.5)
        end

        @testset "Testing zrotate3d!" begin
            zrotate3d!(foo, π / 4)
            @test isapprox(minimum(BMO.vertices(foo)[:, 1]), -0.5)
            @test isapprox(minimum(BMO.vertices(foo)[:, 2]), -0.5)
            @test isapprox(minimum(BMO.vertices(foo)[:, 3]), -0.5)
            @test isapprox(maximum(BMO.vertices(foo)[:, 1]), 0.5)
            @test isapprox(maximum(BMO.vertices(foo)[:, 2]), 0.5)
            @test isapprox(maximum(BMO.vertices(foo)[:, 3]), 0.5)
        end

        # Testing orientation of dir matrix
        @test orientation(foo)[[3, 5, 7]] == [-1, 1, 1]
    end

    # center bar reference cube at origin
    translate3d!(bar, to_origin)
    BMO.set_new_origin3d!(bar)

    @testset "Testing reset_rotation3d!" begin
        translate3d!(foo, [1, 2, 3])
        reset_translation3d!(foo)
        reset_rotation3d!(foo)
        @test position(foo) == zeros(3)
        @test orientation(foo) ≈ orientation(bar)
        @test BMO.vertices(foo) ≈ BMO.vertices(bar)
    end

    @testset "Testing align3d!" begin
        align3d!(foo, normalize([0, 1, 1]))
        @test position(foo) == zeros(3)
        @test orientation(foo)[:, 1] ≈ [1, 0, 0]
        @test orientation(foo)[:, 2] ≈ [0, √2 / 2, √2 / 2]
        @test orientation(foo)[:, 3] ≈ [0, -√2 / 2, √2 / 2]
        reset_rotation3d!(foo)
    end

    @testset "Testing normal" begin
        normal = BMO.normal3d(foo, 1)
        @test isapprox(normal, [0, 0, -1])
    end

    @testset "Testing scale3d!" begin
        BMO.scale3d!(foo, 2)
        @test isapprox(minimum(BMO.vertices(foo)[:, 1]), -1)
        @test isapprox(minimum(BMO.vertices(foo)[:, 2]), -1)
        @test isapprox(minimum(BMO.vertices(foo)[:, 3]), -1)
        @test isapprox(maximum(BMO.vertices(foo)[:, 1]), 1)
        @test isapprox(maximum(BMO.vertices(foo)[:, 2]), 1)
        @test isapprox(maximum(BMO.vertices(foo)[:, 3]), 1)
        @test BMO.scale(foo) == 2
    end

    @testset "Testing Moeller-Trumbore algorithm" begin
        t = 5
        face = [1 1 t
                -1 1 t
                0 -1 t]
        # ray at origin pointing along z-axis
        pos = [0.0, 0, 0]
        dir = [0.0, 0, 1]
        ray = Ray(pos, dir)
        # Preallocate memory
        @test isapprox(BMO.MoellerTrumboreAlgorithm(face, ray), t)
        # Check allocations (WARNING: function must have been compiled once for before this test!)
        alloc = @allocated BMO.MoellerTrumboreAlgorithm(face, ray)
        if alloc > 16
            @warn "Allocated number of bytes for MTA larger than expected!" alloc
        end
    end
    @testset "Testing intersect3d" begin
        # Setup test cube and ray
        cube = BMO.CubeMesh(1)
        translate3d!(cube, -0.5 * [1, 1, 1])
        BMO.set_new_origin3d!(cube)
        ray_pos = zeros(3)
        ray_dir = [1.0, 0, 0]
        ray = Ray(ray_pos, ray_dir)
        # Rotate cube 360°, calculate intersection distance
        θ = 0:1:359
        l = zeros(length(θ))
        for (i, ~) in enumerate(θ)
            intersection = BMO.intersect3d(cube, ray)
            l[i] = length(intersection)
            zrotate3d!(cube, deg2rad(step(θ)))
        end
        # Test if 0/45° distances are correct
        @test all(l[1:90:end] .≈ BMO.scale(cube) * 1 / 2)
        @test all(l[(1 + 45):90:end] .≈ BMO.scale(cube) * sqrt(2) / 2)
    end

    @testset "Testing intersect3d - part 2" begin
        t = 5
        s = 1 # scale/2
        cube = BMO.CubeMesh(2 * s)
        # Move cube COG to origin
        translate3d!(cube, -[s, s, s])
        BMO.set_new_origin3d!(cube)
        # Align cube edge at t units from origin
        translate3d!(cube, [t + s, 0, 0])
        pos = [0, 0, 0]
        steps = 10
        for z in (-s):(s / steps):s
            # Ray constructed each time for unit-length dir
            dir = [t, 0, z]
            ray = Ray(pos, dir)
            @test isapprox(BMO.intersect3d(cube, ray).t, sqrt(t^2 + z^2))
        end
    end

    @testset "Testing constructors" begin
        @testset "Testing RectangularFlatMesh" begin
            rfm = BMO.RectangularFlatMesh(2.0, 1)
            @test BMO.vertices(rfm) == [1 0 0.5; 1 0 -0.5; -1 0 -0.5; -1 0 0.5]
            @test BMO.normal3d(rfm, 1) == [0, 1, 0]
        end

        @testset "Testing QuadraticFlatMesh" begin
            qfm = BMO.QuadraticFlatMesh(4.0)
            @test BMO.vertices(qfm) == [2 0 2; 2 0 -2; -2 0 -2; -2 0 2]
            @test BMO.normal3d(qfm, 1) == [0, 1, 0]
        end
    end

    @testset "Testing CircularFlatMesh" begin
        # Testing constructor
        n = 4
        cm = BMO.CircularFlatMesh(1.0f0, n)
        v = BMO.vertices(cm)
        f = BMO.faces(cm)

        # testing vertices
        @test v[1, :] ≈ zeros(3)
        @test v[2, :] ≈ [1, 0, 0]
        @test v[3, :] ≈ [0, 0, 1]
        @test v[4, :] ≈ [-1, 0, 0]
        @test v[5, :] ≈ [0, 0, -1]

        # testing faces
        @test f[:, 1] == ones(4)
        @test f[:, 2] == [2, 3, 4, 5]
        @test f[:, 3] == [3, 4, 5, 2]

        # testing normal vectors
        for i in 1:n
            @test BMO.normal3d(cm, i) ≈ [0, -1, 0]
        end
    end
end

end # MODULE