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
        @testset "Matrix rotation updates geometry about the mesh position" begin
            mesh = BMO.CubeMesh(1)
            translate3d!(mesh, [2, 3, 4])
            pivot = copy(position(mesh))
            vertex = copy(BMO.vertices(mesh)[1, :])
            normal = BMO.normal3d(mesh, 1)
            R = BMO.rotate3d([0, 1, 0], π / 3)

            rotate3d!(mesh, R)

            @test position(mesh) == pivot
            @test BMO.vertices(mesh)[1, :] ≈ pivot + R * (vertex - pivot)
            @test BMO.normal3d(mesh, 1) ≈ R * normal
            @test orientation(mesh) ≈ R
        end

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

    @testset "align3d! from rotated start orientation" begin
        m = BMO.CubeMesh(1)
        translate3d!(m, [1, 2, 3])
        xrotate3d!(m, π / 5)
        V0 = copy(BMO.vertices(m))
        O0 = orientation(m)
        p0 = copy(position(m))
        t = [1.0, 0, 0]
        R = BMO.align3d(O0[:, 2], t)
        align3d!(m, t)
        @test orientation(m)[:, 2] ≈ t
        @test orientation(m) ≈ R * O0
        @test BMO.vertices(m) ≈ (V0 .- p0') * R' .+ p0'
        @test position(m) == p0
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

@testset "ObjectGroup of meshes" begin
    # Part 1: rotate3d! propagates to member mesh vertices/orientation about the group center
    m1 = SquarePlanoMirror(0.0254, 0.005)
    m2 = SquarePlanoMirror(0.0254, 0.005)
    translate3d!(m1, [0.05, 0, 0])
    translate3d!(m2, [0, 0, 0.05])
    g = ObjectGroup([m1, m2])
    members = (m1, m2)
    pg = copy(position(g))
    V0 = [copy(BMO.vertices(BMO.shape(mi))) for mi in members]
    O0 = [orientation(BMO.shape(mi)) for mi in members]

    rotate3d!(g, [0, 0, 1], π / 3)
    R = BMO.rotate3d([0, 0, 1], π / 3)
    for (i, mi) in enumerate(members)
        @test BMO.vertices(BMO.shape(mi)) ≈ (V0[i] .- pg') * R' .+ pg'
        @test orientation(BMO.shape(mi)) ≈ R * O0[i]
    end

    # Part 2: align3d! propagates the same way and sets the group's local y-axis
    m1b = SquarePlanoMirror(0.0254, 0.005)
    m2b = SquarePlanoMirror(0.0254, 0.005)
    translate3d!(m1b, [0.05, 0, 0])
    translate3d!(m2b, [0, 0, 0.05])
    gb = ObjectGroup([m1b, m2b])
    membersb = (m1b, m2b)
    pgb = copy(position(gb))
    V0b = [copy(BMO.vertices(BMO.shape(mi))) for mi in membersb]
    O0b = [orientation(BMO.shape(mi)) for mi in membersb]

    Og0 = orientation(gb)
    target = normalize([1.0, 1.0, 0.0])
    Rb = BMO.align3d(Og0[:, 2], target)
    align3d!(gb, target)
    for (i, mi) in enumerate(membersb)
        @test BMO.vertices(BMO.shape(mi)) ≈ (V0b[i] .- pgb') * Rb' .+ pgb'
        @test orientation(BMO.shape(mi)) ≈ Rb * O0b[i]
    end
    @test orientation(gb)[:, 2] ≈ target
end

end # MODULE
