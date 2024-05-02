using SCDI
using Test
using LinearAlgebra
using GeometryBasics
using UUIDs
using AbstractTrees

@testset "Utilities" begin
    @testset "Testing normal3d" begin
        @testset "normal3d(v) with Array" begin
            v = [1., 1, 1]
            k = SCDI.normal3d(v)
            @test dot(v, k) ≈ 0 atol=1e-14
            @test norm(k) ≈ 1
        end
    
        @testset "normal3d(v) with Point3" begin
            v = Point3(1., 1, 1)
            k = SCDI.normal3d(v)
            @test dot(v, k) ≈ 0 atol=1e-14
            @test norm(k) ≈ 1
        end
    
        @testset "normal3d(v, w) right hand rule and unit length" begin
            orth = SCDI.normal3d([2, 0, 0], [0, 0, 1])
            @test isapprox(orth, [0, -1, 0])
        end
    end

    @testset "Testing rotate3d for clockwise dir. and conservation of length" begin
        Rot = SCDI.rotate3d([0, 0, 1], π / 2)
        @test isapprox(Rot * [1, 0, 0], [0, 1, 0])
    end

    @testset "Testing align3d for rotation and conservation of length" begin
        # Start vector must have unit length!
        start = [1, 0, 0]
        # Test parallel case
        target = [1, 0, 0]
        T = SCDI.align3d(start, target)
        @test T * start ≈ target
        # Test parallel opposite case
        target = [-1, 0, 0]
        T = SCDI.align3d(start, target)
        @test T * start ≈ target
        # Test norm and 45° rotation
        target = [1.0, 1.0, 0.0]
        T = SCDI.align3d(start, target)
        @test T * start ≈ normalize(target)
    end

    @debug "Testing angle3d for resulting angle"
    a = SCDI.angle3d([1, 0, 0], [0, 0, 1])
    @test isapprox(a, π / 2)

    @testset "Testing line_point_distance3d and isinfrontof" begin
        pos = [0, 0, 0]
        dir = [1, 0, 0]
        point = [5, 1, 1]
        d = SCDI.line_point_distance3d(pos, dir, point)
        @test isapprox(d, √2)

        @test SCDI.isinfrontof(point, pos, dir) == true
    end

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
            for θ1 in 0:(π / 8):(π / 2)
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
            for θ1 in 0:(π / 8):(π / 2)
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

    @testset "Fresnel equations" begin
        @testset "Vacuum-glass: normal incidence" begin
            n = 1.5
            θ = 0.0
            rs, rp, ts, tp = SCDI.fresnel_coefficients(θ, n)
            @test real(rs) ≈ (1-n)/(1+n)
            @test real(rs) ≈ real(rp)
            @test real(tp) ≈ 2/(1+n)
            @test real(tp) ≈ real(ts)
        end
        
        @testset "Vacuum-glass: Brewster angle" begin
            n = 1.5
            θb = atan(n)
            rs, rp, ts, tp = SCDI.fresnel_coefficients(θb, n)
            @test real(rp) ≈ 0
        end
        
        @testset "Vacuum-glass: grazing incidence" begin
            n = 1.5
            θ = π/2
            rs, rp, ts, tp = SCDI.fresnel_coefficients(θ, n)
            @test real(rs) ≈ -1
            @test real(rp) ≈ 1
            @test real(ts) ≈ 0
            @test real(tp) ≈ 0 atol=2e-16
        end
        
        @testset "Glass-vacuum: normal incidence" begin
            n = 1/1.5
            θ = 0.0
            rs, rp, ts, tp = SCDI.fresnel_coefficients(θ, n)
            @test real(rs) ≈ (1-n)/(1+n)
            @test real(rs) ≈ real(rp)
            @test real(tp) ≈ 2/(1+n)
            @test real(tp) ≈ real(ts)
        end
        
        @testset "Glass-vacuum: Brewster angle" begin
            n = 1/1.5
            θb = atan(n)
            rs, rp, ts, tp = SCDI.fresnel_coefficients(θb, n)
            @test real(rp) ≈ 0 atol=2e-16
        end
        
        @testset "Glass-vacuum: Total internal reflection" begin
            n = 1/1.5
            θc = asin(n)
            rs, rp, ts, tp = SCDI.fresnel_coefficients(θc, n)
            @test SCDI.is_internally_reflected(rp, rs)
            @test real(rs) ≈ 1
            @test real(rp) ≈ -1
            @test real(ts) ≈ 2
            @test real(tp) ≈ 3 atol=1e-15
        end
    end
end

@testset "Types" begin
    @test isdefined(SCDI, :AbstractEntity)
    @test isdefined(SCDI, :AbstractSystem)
    @test isdefined(SCDI, :AbstractRay)
    @test isdefined(SCDI, :AbstractBeam)
    @test isdefined(SCDI, :AbstractShape)
    @test isdefined(SCDI, :AbstractObject)
    @test isdefined(SCDI, :AbstractObjectGroup)
    @test isdefined(SCDI, :Intersection)
    @test isdefined(SCDI, :AbstractInteraction)

    # Generate test structs
    struct TestEntity <: SCDI.AbstractEntity
        id::UUID
    end

    @testset "AbstractEntity" begin
        e = TestEntity(uuid4())
        @test SCDI.id(e) == e.id
        @test SCDI.id("Hello World") === nothing
    end

    struct TestSystem <: SCDI.AbstractSystem end

    @testset "AbstractSystem" begin
        # no tests
        sys = TestSystem()
    end

    mutable struct TestRay{T} <: SCDI.AbstractRay{T}
        pos::Vector{T}
        dir::Vector{T}
        λ::T
        n::T
    end

    TestRay(pos::AbstractArray{T}, dir::AbstractArray{T}) where T = TestRay(pos, dir, T(1000e-9), T(1))

    Base.length(::TestRay) = π

    @testset "AbstractRay" begin
        r = TestRay([0.0, 0, 0], [1.0, 0, 0])
        # Test getters
        @test SCDI.position(r) == r.pos
        @test SCDI.direction(r) == r.dir
        @test SCDI.wavelength(r) == r.λ
        @test SCDI.refractive_index(r) == r.n
        # Test setters
        n_pos = [1, 1, 1]
        n_dir = [1.0, 1, 0]
        n_lam = 532e-9
        n_rfi = 1.5
        SCDI.position!(r, n_pos)
        SCDI.direction!(r, n_dir)
        SCDI.wavelength!(r, n_lam)
        SCDI.refractive_index!(r, n_rfi)
        @test SCDI.position(r) == n_pos
        @test SCDI.direction(r) ≈ n_dir .* (sqrt(2) / 2)
        @test SCDI.wavelength(r) == n_lam
        @test SCDI.refractive_index(r) == n_rfi
        @test length(r) == π
        # Test ray-plane intersection
        plane_pos = [1, 0, -1]
        plane_nml_1 = [-1, 0, 1]
        plane_nml_2 = [1, 0, 0]
        plane_nml_3 = [0, 0, 1]
        ray = TestRay([0.0, 0, 0], [1.0, 0, 0])
        is_1 = SCDI.intersect3d(plane_pos, plane_nml_1, ray)
        is_2 = SCDI.intersect3d(plane_pos, plane_nml_2, ray)
        is_3 = SCDI.intersect3d(plane_pos, plane_nml_3, ray)
        @test length(is_1) == 2
        @test length(is_2) == 1
        @test isnothing(is_3)
    end

    mutable struct TestBeam{T} <: SCDI.AbstractBeam{T, TestRay{T}}
        id::UUID
        parent::SCDI.Nullable{TestBeam}
        children::Vector{TestBeam}
    end

    TestBeam() = TestBeam{Float64}(uuid4(), nothing, Vector{TestBeam{Float64}}())

    @testset "AbstractBeam" begin
        # Create beam tree
        root = TestBeam()
        cb1 = TestBeam()
        cb2 = TestBeam()
        group = [cb1, cb2]
        cb3 = TestBeam()
        # Add children to root, cb2
        SCDI.children!(root, group)
        SCDI.children!(cb2, cb3)
        # Test tree structure
        @test treeheight(root) == 2
        @test treebreadth(root) == 2
        @test SCDI.children(root) == group
        @test first(SCDI.children(cb2)) === cb3
        # Test parent connection
        @test AbstractTrees.parent(root) === nothing
        @test AbstractTrees.parent(cb1) === root
        @test AbstractTrees.parent(cb2) === root
        @test AbstractTrees.parent(cb3) === cb2
        # Replace bottom child
        cbr = TestBeam()
        @test_throws "_modify_beam_head not implemented for $(typeof(cb2))" SCDI.children!(cb2,
            cbr)
        # Test child removal
        SCDI._drop_beams!(cb2)
        @test isempty(SCDI.children(cb2))
        # Stuff
        @test_throws "_last_beam_intersection not implemented for $(typeof(cb2))" SCDI._last_beam_intersection(cb2)
    end

    mutable struct TestShapeless{T} <: SCDI.AbstractShape{T}
        pos::Vector{T}
        dir::Matrix{T}
    end

    SCDI.id(::TestShapeless) = nothing

    TestShapeless() = TestShapeless{Float64}(zeros(3), Matrix{Float64}(I, 3, 3))

    @testset "AbstractShape" begin
        pos = zeros(3)
        dir = Matrix{Float64}(I, 3, 3)
        shape = TestShapeless(pos, dir)
        # Test get/set
        @test SCDI.position(shape) == pos
        @test SCDI.orientation(shape) == dir
        n_pos = [1, 1, 1]
        n_dir = SCDI.rotate3d([0, 0, 1], π / 4)
        SCDI.position!(shape, n_pos)
        SCDI.orientation!(shape, n_dir)
        @test SCDI.position(shape) == n_pos
        @test SCDI.orientation(shape) == n_dir
        # Test translation
        SCDI.translate3d!(shape, n_pos)
        @test SCDI.position(shape) == 2 * n_pos
        SCDI.reset_translation3d!(shape)
        @test SCDI.position(shape) == zeros(3)
        # Test rotation
        dir = Matrix{Float64}(I, 3, 3)
        SCDI.orientation!(shape, dir)
        SCDI.rotate3d!(shape, [0, 0, 1], π)
        @test SCDI.orientation(shape)[1:4:9] == [-1, -1, 1]
        SCDI.xrotate3d!(shape, π)
        SCDI.yrotate3d!(shape, π)
        SCDI.zrotate3d!(shape, π)
        SCDI.reset_rotation3d!(shape)
        @test SCDI.orientation(shape) == dir
        # The following test are expected to do nothing but not throw exceptions
        ray = TestRay([0.0, 0, 0], [1.0, 0, 0])
        @test_logs (:warn, "No intersect3d method defined for:") SCDI.intersect3d(shape,
            ray)

        @testset "Testing AbstractRay - AbstractShape" begin
            shape = TestShapeless([1, 0, 0], Matrix{Int}(I, 3, 3))
            ray = TestRay([0.0, 0, 0], [1.0, 0, 0])
            @test SCDI.isinfrontof(shape, ray) == true
            SCDI.direction!(ray, -[1, 0, 0])
            @test SCDI.isinfrontof(shape, ray) == false
            SCDI.direction!(ray, [0, 1, 0])
            @test SCDI.isinfrontof(shape, ray) == false
            SCDI.direction!(ray, [1.0, 1, 0])
            @test SCDI.isinfrontof(shape, ray) == true
            @test norm(SCDI.direction(ray)) ≈ 1
        end
    end

    struct TestObject <: SCDI.AbstractObject
        id::UUID
        shape::TestShapeless{<:Real}
    end

    TestObject() = TestObject(uuid4(), TestShapeless())

    @testset "AbstractObject" begin
        object = TestObject()
        @test isa(SCDI.shape(object), TestShapeless)
        # Test forwarding of kin. API to object shape
        @test SCDI.position(object) == SCDI.position(SCDI.shape(object))
        @test SCDI.position(object) == SCDI.position(SCDI.shape(object))
        SCDI.translate3d!(object, ones(3))
        SCDI.rotate3d!(object, [0, 0, 1], π)
        @test SCDI.position(object) == ones(3)
        @test SCDI.orientation(object)[1:4:9] == [-1, -1, 1]
        SCDI.reset_translation3d!(object)
        SCDI.reset_rotation3d!(object)
        @test SCDI.position(object) == zeros(3)
        @test SCDI.orientation(object)[1:4:9] == ones(3)
        # Test translate_to3d
        target_pos = [1, 3, 9]
        SCDI.translate_to3d!(object, target_pos)
        @test SCDI.position(object) == target_pos

        @testset "Testing interact3d" begin
            sys = TestSystem()
            obj = TestObject()
            ray = TestRay(zeros(3), ones(3))
            beam = TestBeam()
            @test_logs (:warn, "No interact3d method defined for:") SCDI.interact3d(sys,
                obj,
                beam,
                ray)===nothing
        end
    end
end

@testset "Rays" begin
    # Testing constructor
    pos = [0, 0, 0]
    dir = [1.0, 1, 0]
    ray = SCDI.Ray(pos, dir)
    @test ismutable(ray)
    @test isa(ray, SCDI.Ray{Float64})
    @test isnothing(ray.intersection)
    @test isinf(length(ray))
    @test isapprox(norm(ray.dir), 1)
    # Test helper functions
    @test SCDI.line_point_distance3d(ray, [1, 1, 0]) == 0
    @test SCDI.line_point_distance3d(ray, [-1, 1, 0]) == sqrt(2)
end

@testset "Beams" begin
    is = SCDI.Intersection(1.0, zeros(3), uuid4())
    r1 = SCDI.Ray([0.0, 0, 0], [1, 0, 0])
    r2 = SCDI.Ray([1.0, 0, 0], [0, 1, 0])
    r3 = SCDI.Ray([1.0, 1, 0], [0, 0, 1])
    r4 = SCDI.Ray([1.0, 1, 1], [1, 0, 0])
    SCDI.intersection!(r1, is)
    SCDI.intersection!(r2, is)
    SCDI.intersection!(r3, is)
    SCDI.intersection!(r4, nothing)
    # Test beam
    beam = SCDI.Beam(r1)
    push!(beam, r2)
    push!(beam, r3)
    push!(beam, r4)
    @test length(beam) == 3
    @test SCDI.point_on_beam(beam, 0) == ([0, 0, 0], 1)
    @test SCDI.point_on_beam(beam, 1) == ([1, 0, 0], 2)
    @test SCDI.point_on_beam(beam, 2) == ([1, 1, 0], 3)
    @test SCDI.point_on_beam(beam, 3) == ([1, 1, 1], 4)
    @test SCDI.point_on_beam(beam, 10) == ([8, 1, 1], 4)
    @test SCDI.isparentbeam(beam, r2) == true
end

@testset "Mesh" begin
    # NOTE: the "Mesh" testset is mutating. Errors/fails might lead to subsequent tests failing too!
    @test isdefined(SCDI, :AbstractMesh)
    @test isdefined(SCDI, :Mesh)

    # Generate cube since types are defined
    foo = SCDI.CubeMesh(1) # test cube

    @testset "Testing AbstractMesh getters" begin
        @test typeof(foo) == SCDI.Mesh{Float64}
        @test SCDI.vertices(foo) == foo.vertices
        @test SCDI.faces(foo) == foo.faces
        @test SCDI.orientation(foo) == foo.dir
        @test SCDI.position(foo) == foo.pos
        @test SCDI.scale(foo) == foo.scale
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
        @testset "Testing rotate3d!" begin
            SCDI.rotate3d!(foo, [1, 0, 0], π / 4)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 1]), -0.5)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 2]), -√2 / 2)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 3]), -√2 / 2)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 1]), 0.5)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 2]), √2 / 2)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 3]), √2 / 2)
            # Return to original rotation
            SCDI.rotate3d!(foo, [1, 0, 0], -π / 4)
        end

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

        # Testing orientation of dir matrix
        @test all(SCDI.orientation(foo)[[2, 6, 7]] .== 1)
    end

    @testset "Testing reset_rotation3d!" begin
        SCDI.reset_rotation3d!(foo)
        SCDI.reset_translation3d!(foo) # sneaky, useless command for that sweet code coverage
        @test all(SCDI.orientation(foo) .== Matrix{Float64}(I, 3, 3))
    end

    @testset "Testing normal" begin
        normal = SCDI.normal3d(foo, 1)
        @test isapprox(normal, [0, -1, 0])
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

    @testset "Testing Moeller-Trumbore algorithm" begin
        t = 5
        face = [1 1 t
            -1 1 t
            0 -1 t]
        # ray at origin pointing along z-axis
        pos = [0.0, 0, 0]
        dir = [0.0, 0, 1]
        ray = SCDI.Ray(pos, dir)
        # Preallocate memory
        @test isapprox(SCDI.MoellerTrumboreAlgorithm(face, ray), t)
        # Check allocations (WARNING: function must have been compiled once for before this test!)
        alloc = @allocated SCDI.MoellerTrumboreAlgorithm(face, ray)
        if alloc > 16
            @warn "Allocated number of bytes for MTA larger than expected!" alloc
        end
    end
    @testset "Testing intersect3d" begin
        # Setup test cube and ray
        cube = SCDI.CubeMesh(1)
        SCDI.translate3d!(cube, -0.5 * [1, 1, 1])
        SCDI.set_new_origin3d!(cube)
        ray_pos = zeros(3)
        ray_dir = [1.0, 0, 0]
        ray = SCDI.Ray(ray_pos, ray_dir)
        # Rotate cube 360°, calculate intersection distance
        θ = 0:1:359
        l = zeros(length(θ))
        for (i, ~) in enumerate(θ)
            intersection = SCDI.intersect3d(cube, ray)
            l[i] = length(intersection)
            SCDI.zrotate3d!(cube, deg2rad(step(θ)))
        end
        # Test if 0/45° distances are correct
        @test all(l[1:90:end] .≈ SCDI.scale(cube) * 1 / 2)
        @test all(l[(1 + 45):90:end] .≈ SCDI.scale(cube) * sqrt(2) / 2)
    end

    @testset "Testing intersect3d - part 2" begin
        t = 5
        s = 1 # scale/2
        cube = SCDI.CubeMesh(2 * s)
        # Move cube COG to origin
        SCDI.translate3d!(cube, -[s, s, s])
        SCDI.set_new_origin3d!(cube)
        # Align cube edge at t units from origin
        SCDI.translate3d!(cube, [t + s, 0, 0])
        pos = [0, 0, 0]
        steps = 10
        for z in (-s):(s / steps):s
            # Ray constructed each time for unit-length dir
            dir = [t, 0, z]
            ray = SCDI.Ray(pos, dir)
            @test isapprox(SCDI.intersect3d(cube, ray).t, sqrt(t^2 + z^2))
        end
    end
end

@testset "SDFs" begin
    @testset "Testing type definitions" begin
        @test isdefined(SCDI, :AbstractSDF)
        @test isdefined(SCDI, :SphereSDF)
        @test isdefined(SCDI, :CylinderSDF)
        @test isdefined(SCDI, :CutSphereSDF)
        @test isdefined(SCDI, :ThinLensSDF)
    end

    # Orientation-less test point sdf
    mutable struct TestPointSDF{T} <: SCDI.AbstractSDF{T}
        position::Point3{T}
        orientation::Matrix{T}
    end
    
    TestPointSDF(p::AbstractArray{T}) where {T} = TestPointSDF{T}(Point3{T}(p), Matrix{T}(I, 3, 3))
    TestPointSDF(T = Float64) = TestPointSDF{T}(Point3{T}(0), Matrix{T}(I, 3, 3))
    
    SCDI.position(tps::TestPointSDF) = tps.position
    SCDI.position!(tps::TestPointSDF{T}, new::Point3{T}) where T = (tps.position = new)
    
    SCDI.orientation(tps::TestPointSDF) = tps.orientation
    SCDI.orientation!(tps::TestPointSDF{T}, new::Matrix{T}) where T = (tps.orientation = new)
    
    SCDI.transposed_orientation(tps::TestPointSDF) = transpose(tps.orientation)
    SCDI.transposed_orientation!(::TestPointSDF, ::Any) = nothing

    orientation(::TestPointSDF{T}) where {T} = Matrix{T}(I, 3, 3)
    orientation!(::TestPointSDF, ::Any) = nothing

    function SCDI.sdf(tps::TestPointSDF, point)
        p = SCDI._world_to_sdf(tps, point)
        return norm(p)
    end

    @testset "Testing kinematics and transforms" begin
        point = TestPointSDF()
        t = 10
        θ = deg2rad(30)
        SCDI.translate3d!(point, [t, 0, 0])
        SCDI.rotate3d!(point, [0,1,0], θ)
        pt = SCDI._world_to_sdf(point, [0,0,0])
        @test pt[1] ≈ -t * cos(θ)
        @test pt[2] ≈ 0
        @test pt[3] ≈ -t * sin(θ)
    end

    @testset "Testing intersect3d" begin
        t = 10.0
        point = TestPointSDF(zeros(3))
        SCDI.translate3d!(point, [t, 0, 0])

        r1 = SCDI.Ray(zeros(3), [1.0, 0, 0])
        r2 = SCDI.Ray(zeros(3), [1.0, 1, 0])
        r3 = SCDI.Ray(zeros(3), [1.0, 0, 1])

        i1 = SCDI.intersect3d(point, r1)
        i2 = SCDI.intersect3d(point, r2)
        i3 = SCDI.intersect3d(point, r3)

        @test length(i1) == t
        @test isnothing(i2)
        @test isnothing(i3)
    end

    @testset "Testing normal3d" begin
        point = TestPointSDF(zeros(3))
        offset = [5, 0, 0]
        SCDI.translate3d!(point, offset)
        p1 = [1, 0, 0]
        p2 = [0, 1, 0]
        p3 = [0, 0, 1]
        @test SCDI.normal3d(point, p1 + offset) == p1
        @test SCDI.normal3d(point, p2 + offset) == p2
        @test SCDI.normal3d(point, p3 + offset) == p3
    end
end

@testset "System" begin
    @testset "Testing implementation" begin
        struct SystemTestBeam{T} <: SCDI.AbstractBeam{T, SCDI.Ray{T}} end
        struct SystemTestObject <: SCDI.AbstractObject
            id::UUID
        end
        o1 = SystemTestObject(uuid4())
        o2 = SystemTestObject(uuid4())
        system = SCDI.System(o1)
        beam = SystemTestBeam{Real}()
        # Test missing implementation warnings
        @test_logs (:warn, "Tracing for $(typeof(beam)) not implemented") SCDI.trace_system!(system,
            beam)
        @test_logs (:warn, "Retracing for $(typeof(beam)) not implemented") SCDI.retrace_system!(system,
            beam)
        # Test getters
        @test SCDI.object(system, SCDI.id(o1)) == o1
        @test isnothing(SCDI.object(system, SCDI.id(o2)))
    end

    # Setup circular multipass cell with flat mirrors
    n_mirrors = 101
    radius = 1
    L = 6 * radius / n_mirrors
    Δθ = 360 / (n_mirrors + 1)
    mirrors = [SCDI.RectangularPlanoMirror2D(L) for _ in 1:n_mirrors]
    θ = 1 * Δθ
    for m in mirrors
        point = radius * [cos(deg2rad(θ)), sin(deg2rad(θ)), 0]
        SCDI.zrotate3d!(m, deg2rad(θ))
        SCDI.translate3d!(m, point)
        θ += Δθ
    end
    SCDI.zrotate3d!.(mirrors, deg2rad(90))

    # Initial ray orientation and position
    dir = [-1, 0, 0]
    Rot = SCDI.rotate3d([0, 0, 1], deg2rad(Δθ * 1))
    dir = Vector(Rot * dir)
    origin = [radius, 0, 0] + -1 * dir

    @testset "Testing tracing subroutines" begin
        system = SCDI.System(mirrors)
        ray = SCDI.Ray(origin, dir)
        first_id = SCDI.id(mirrors[(n_mirrors + 1) ÷ 2 + 2])
        false_id = SCDI.id(mirrors[(n_mirrors + 1) ÷ 2 + 2 + 1])
        # trace_all
        @test SCDI.hasid(SCDI.trace_all(system, ray), first_id)
        # trace_one
        @test SCDI.hasid(SCDI.trace_one(system, ray, first_id), first_id)
        @test SCDI.hasid(SCDI.trace_one(system, ray, false_id), first_id)
        # tracing step
        SCDI.tracing_step!(system, ray, nothing)
        @test SCDI.hasid(SCDI.intersection(ray), first_id)
    end

    @testset "Testing system tracing" begin
        system = SCDI.System(mirrors)
        first_ray = SCDI.Ray(origin, dir)
        beam = SCDI.Beam(first_ray)
        # Test trace_system!
        nmax = 10
        SCDI.trace_system!(system, beam, r_max = nmax)
        @test length(SCDI.rays(beam)) == nmax
        SCDI.trace_system!(system, beam, r_max = 1000000)
        @test length(SCDI.rays(beam)) == n_mirrors + 1
        first_ray_dir = SCDI.direction(first_ray)
        last_ray_dir = SCDI.direction(last(SCDI.rays(beam)))
        @test 180 - rad2deg(SCDI.angle3d(first_ray_dir, last_ray_dir)) ≈ 2 * Δθ
        @test SCDI.id(SCDI.intersection(first_ray)) ==
              SCDI.id(mirrors[(n_mirrors + 1) ÷ 2 + 2])
    end

    @testset "Testing StaticSystem tracing" begin
        # same testset as before
        system = SCDI.StaticSystem(mirrors)
        first_ray = SCDI.Ray(origin, dir)
        beam = SCDI.Beam(first_ray)
        # Test trace_system!
        nmax = 10
        SCDI.trace_system!(system, beam, r_max = nmax)
        @test length(SCDI.rays(beam)) == nmax
        SCDI.trace_system!(system, beam, r_max = 1000000)
        @test length(SCDI.rays(beam)) == n_mirrors + 1
        first_ray_dir = SCDI.direction(first_ray)
        last_ray_dir = SCDI.direction(last(SCDI.rays(beam)))
        @test 180 - rad2deg(SCDI.angle3d(first_ray_dir, last_ray_dir)) ≈ 2 * Δθ
        @test SCDI.id(SCDI.intersection(first_ray)) ==
              SCDI.id(mirrors[(n_mirrors + 1) ÷ 2 + 2])
    end

    @testset "Testing system retracing" begin
        system = SCDI.System(mirrors)
        first_ray = SCDI.Ray(origin, dir)
        beam = SCDI.Beam(first_ray)
        t1 = @timed SCDI.trace_system!(system, beam, r_max = 1000000)
        t2 = @timed SCDI.retrace_system!(system, beam) # for precompilation
        t2 = @timed SCDI.retrace_system!(system, beam)
        if t1.time < t2.time
            @warn "Retracing took longer than tracing, something might be bugged...\n   Tracing: $(t1.time) s\n   Retracing: $(t2.time) s"
        end
    end
end

@testset "Object groups" begin
    mutable struct TestPoint{T} <: SCDI.AbstractShape{T}
        pos::Point3{T}
        dir::Matrix{T}
    end

    TestPoint(position::AbstractArray{T}) where {T <: Real} = TestPoint{T}(Point3{T}(position),
        Matrix{T}(I, 3, 3))

    struct GroupTestObject <: SCDI.AbstractObject
        id::UUID
        shape::TestPoint{<:Real}
    end

    GroupTestObject(position::AbstractArray) = GroupTestObject(uuid4(), TestPoint(position))

    n = 8
    xs = [cos(x) for x in LinRange(0, 2pi * (n - 1) / n, n)]
    ys = [sin(x) for x in LinRange(0, 2pi * (n - 1) / n, n)]

    center = GroupTestObject(zeros(3))
    circle = SCDI.ObjectGroup([GroupTestObject([xs[i], ys[i], 0]) for i in eachindex(xs)])

    objects = SCDI.ObjectGroup([center, circle])

    # Translation test
    target = [3, 0, 0]
    SCDI.translate_to3d!(objects, target)

    @testset "translate3d" begin
        # Test if all objects/subgroups have been translated
        @test SCDI.position(objects) == target
        @test SCDI.position(center) == target
        @test SCDI.position(circle) == target
        for (i, obj) in enumerate(SCDI.objects(circle))
            @test SCDI.position(obj) == [xs[i], ys[i], 0] + target
        end
    end

    # Rotation test
    angle = 2π / n
    SCDI.rotate3d!(objects, [0, 0, 1], angle)

    @testset "rotate3d" begin
        # Test if all objects/subgroups have been rotated relative to the origin
        Rt = SCDI.rotate3d([0, 0, 1], angle)
        xt = circshift(xs, -1)
        yt = circshift(ys, -1)
        @test SCDI.orientation(objects) == Rt
        @test SCDI.orientation(center) == Rt
        @test SCDI.orientation(circle) == Rt
        for (i, obj) in enumerate(SCDI.objects(circle))
            @test SCDI.orientation(obj) == Rt
            @test SCDI.position(obj) ≈ [xt[i], yt[i], 0] + target
        end
    end

    # Reset test
    SCDI.reset_translation3d!(objects)
    SCDI.reset_rotation3d!(objects)

    @testset "reset functions" begin
        Ri = Matrix{Float64}(I, 3, 3)
        # Test if objects are reset correctly
        @test SCDI.position(objects) == zeros(3)
        @test SCDI.position(center) == zeros(3)
        @test SCDI.position(circle) == zeros(3)
        @test SCDI.orientation(objects) == Ri
        @test SCDI.orientation(center) == Ri
        @test SCDI.orientation(circle) == Ri
        for (i, obj) in enumerate(SCDI.objects(circle))
            @test SCDI.orientation(obj) == Ri
            @test SCDI.position(obj) == zeros(3)
        end
    end

    @testset "System compatibility" begin
        # Test if objects in ObjectGroup are exposed correctly when iterating
        system = SCDI.System(objects)
        ctr = 0
        # Only the objects within the groups should be exposed
        for obj in SCDI.objects(system)
            @test isa(obj, GroupTestObject)
            ctr += 1
        end
        @test ctr == n + 1
    end
end

@testset "Spherical Lenses" begin
    @testset "Testing type definitions" begin
        @test isdefined(SCDI, :AbstractSDF)
        @test isdefined(SCDI, :SphereSDF)
        @test isdefined(SCDI, :CylinderSDF)
        @test isdefined(SCDI, :CutSphereSDF)
        @test isdefined(SCDI, :ThinLensSDF)
    end

    @testset "Thin lens focal length" begin
        # define thin lens
        R1 = 1
        R2 = 1
        nl = 1.5
        tl = SCDI.ThinLensSDF(R1, R2, 0.1)
        p = SCDI.Lens(uuid4(), tl, x -> 1.5)
        system = SCDI.System(p)

        # compare numerical and analytical focal length
        f_analytical = SCDI.lensmakers_eq(R1, -R2, nl)
        zs = -0.04:0.01:0.04
        for (i, z) in enumerate(zs)
            # skip optical axis ray
            if z ≈ 0
                continue
            end
            xs = 0.1:0.1:1.5
            df = zeros(Float64, length(xs))
            ray = SCDI.Ray([0, -0.5, z], [0, 1, 0], 1e3)
            beam = SCDI.Beam(ray)
            SCDI.solve_system!(system, beam)
            # test if numerical and analytical focal length agree
            for (i, x) in enumerate(xs)
                df[i] = SCDI.line_point_distance3d(beam.rays[end], [0, x, 0])
            end
            @test xs[findmin(df)[2]] ≈ f_analytical
        end
    end

    @testset "Testing lens constructor" begin
        # Test against Thorlab spherical lenses
        r1 = r2 = 34.9e-3
        l = 6.8e-3
        LB1811 = SCDI.BiConvexLensSDF(r1, r2, l)
        lens = SCDI.SphericalLens(r1, r2, l)
        @test typeof(SCDI.shape(lens)) == typeof(LB1811)
        r1 = 15.5e-3
        r2 = Inf
        l = 8.6e-3
        LA1805 = SCDI.PlanoConvexLensSDF(r1, l)
        lens = SCDI.SphericalLens(r1, r2, l)
        @test typeof(SCDI.shape(lens)) == typeof(LA1805)
        r1 = r2 = 52.0e-3
        l = 3e-3
        LD1464 = SCDI.BiConcaveLensSDF(r1, r2, l)
        lens = SCDI.SphericalLens(-r1, -r2, l)
        @test typeof(SCDI.shape(lens)) == typeof(LD1464)
        r1 = 25.7e-3
        r2 = Inf
        l = 3.5e-3
        LC1715 = SCDI.PlanoConcaveLensSDF(r1, l)
        lens = SCDI.SphericalLens(-r1, r2, l)
        @test typeof(SCDI.shape(lens)) == typeof(LC1715)
        r1 = 32.1e-3
        r2 = 82.2e-3
        l = 3.6e-3
        LE1234 = SCDI.ConvexConcaveLensSDF(r1, r2, l)
        lens = SCDI.SphericalLens(r1, -r2, l)
        @test typeof(SCDI.shape(lens)) == typeof(LE1234)
    end

    @testset "Testing spherical lens SDFs" begin
        # Based on https://www.pencilofrays.com/double-gauss-sonnar-comparison/
        l1 = SCDI.SphericalLens(48.88e-3, -182.96e-3, 8.89e-3, 52.3e-3, λ -> 1.62286)
        l2 = SCDI.SphericalLens(36.92e-3, Inf, 15.11e-3, 45.11e-3, λ -> 1.58565)
        l3 = SCDI.SphericalLens(-23.06e-3, Inf, 2.31e-3, 45.11e-3, λ -> 1.67764)
        l4 = SCDI.SphericalLens(-23.91e-3, Inf, 1.92e-3, 40.01e-3, λ -> 1.57046)
        l5 = SCDI.SphericalLens(36.92e-3, Inf, 7.77e-3, 40.01e-3, λ -> 1.64128)
        l6 = SCDI.SphericalLens(48.88e-3, 1063.24e-3, 6.73e-3, 45.11e-3, λ -> 1.62286)
        SCDI.zrotate3d!(SCDI.shape(l1), π)
        SCDI.zrotate3d!(SCDI.shape(l2), π)
        SCDI.zrotate3d!(SCDI.shape(l4), π)
        SCDI.translate3d!(SCDI.shape(l2), [0, 11.495e-3, 0])
        SCDI.translate3d!(SCDI.shape(l3), [0, 25.491e-3, 0])
        SCDI.translate3d!(SCDI.shape(l4), [0, 35.568e-3, 0])
        SCDI.translate3d!(SCDI.shape(l5), [0, 42.876e-3, 0])
        SCDI.translate3d!(SCDI.shape(l6), [0, 50.813e-3, 0])
        system = SCDI.System([l1, l2, l3, l4, l5, l6])
        # Test against back focal length as per source above
        λ = 486.0 # nm
        f0 = [0, 0.11602585097812582, 0] # corresponds to f=59.21 mm
        zs = -0.02:1e-3:0.02
        for (i, z) in enumerate(zs)
            beam = SCDI.Beam(SCDI.Ray([0, -0.05, z], [0, 1, 0], λ))
            SCDI.solve_system!(system, beam)
            ray = beam.rays[end]
            intersection = SCDI.intersect3d(f0, [0, 1, 0], ray)
            f = SCDI.position(ray) + length(intersection) * SCDI.direction(ray)
            dz = f0[3] - f[3]
            # Test correct amount of beam sections
            @test length(beam.rays) == 13
            # Test coma
            @test dz < 7e-5
        end
    end
end

@testset "Gaussian beamlet" begin
    @testset "Testing type definitions" begin
        @test isdefined(SCDI, :GaussianBeamlet)
    end

    @testset "Testing analytical equations" begin
        λ = 500e-9
        w0 = 1e-3
        M2 = 1
        zR = SCDI.rayleigh_range(λ, w0, M2)
        # Test Rayleigh range and div. angle against Paschotta (https://www.rp-photonics.com/gaussian_beams.html)
        @test isapprox(zR, 6.28, atol = 1e-2)
        @test isapprox(SCDI.beam_waist(zR, w0, zR), sqrt(2) * w0)
        @test isapprox(SCDI.gouy_phase(zR, zR), -π / 4)
        @test isapprox(SCDI.wavefront_curvature(zR, zR), 1 / (2 * zR))
        @test isapprox(SCDI.divergence_angle(λ, w0, M2), 159e-6, atol = 1e-6)
    end

    @testset "Testing parameter correctness" begin
        # Gauss beam parameters
        y = -5:0.01:5       # m
        λ_1 = 500e-9        # m
        λ_2 = 1000e-9       # m
        P0 = 1              # W
        r = 0               # m
        w0_1 = 1e-3         # m
        w0_2 = 2e-3         # m
        M2_1 = 1e-3         # m
        M2_2 = 2e-3         # m
        E0_1 = SCDI.electric_field(2 * P0 / (π * w0_1^2))
        E0_2 = SCDI.electric_field(2 * P0 / (π * w0_2^2))
        gauss_1 = SCDI.GaussianBeamlet(SCDI.Ray([0.0, 0, 0], [0.0, 1, 0]),
            λ_1,
            w0_1,
            M2 = M2_1,
            P0 = P0)
        gauss_2 = SCDI.GaussianBeamlet(SCDI.Ray([0.0, 0, 0], [0.0, 1, 0]),
            λ_2,
            w0_2,
            M2 = M2_2,
            P0 = P0)
        # Calculate analytical values
        zr_1 = SCDI.rayleigh_range(λ_1, w0_1, M2_1)
        zr_2 = SCDI.rayleigh_range(λ_2, w0_2, M2_2)
        wa_1 = SCDI.beam_waist.(y, w0_1, zr_1)
        wa_2 = SCDI.beam_waist.(y, w0_2, zr_2)
        Ra_1 = SCDI.wavefront_curvature.(y, zr_1)
        Ra_2 = SCDI.wavefront_curvature.(y, zr_2)
        ψa_1 = SCDI.gouy_phase.(y, zr_1)
        ψa_2 = SCDI.gouy_phase.(y, zr_2)
        Ea_1 = SCDI.electric_field.(r, y, E0_1, w0_1, λ_1, M2_1)
        Ea_2 = SCDI.electric_field.(r, y, E0_2, w0_2, λ_2, M2_2)
        # Calculate numerical values
        wn_1, Rn_1, ψn_1, w0n_1 = SCDI.gauss_parameters(gauss_1, y)
        wn_2, Rn_2, ψn_2, w0n_2 = SCDI.gauss_parameters(gauss_2, y)
        En_1 = [SCDI.electric_field(gauss_1, r, yi) for yi in y]
        En_2 = [SCDI.electric_field(gauss_2, r, yi) for yi in y]
        # Compare beam diameter within 0.1 nm
        @test all(isapprox.(wa_1, wn_1, atol = 1e-10))
        @test all(isapprox.(wa_2, wn_2, atol = 1e-10))
        # Compare wavefront curvature
        @test all(isapprox.(Ra_1, Rn_1, atol = 5e-9))
        @test all(isapprox.(Ra_2, Rn_2, atol = 5e-9))
        # Compare Gouy phase within 0.1 μrad
        @test all(isapprox.(ψa_1, ψn_1, atol = 1e-7))
        @test all(isapprox.(ψa_2, ψn_2, atol = 1e-7))
        # Compare calculated waist radius
        @test all(isapprox.(w0_1, w0n_1))
        @test all(isapprox.(w0_2, w0n_2))
        # Compare calculate electric field at r
        @test all(isapprox.(Ea_1, En_1, atol = 1e-8))
        @test all(isapprox.(Ea_2, En_2, atol = 1e-7))
    end

    @testset "Testing propagation correctness" begin
        # Analytical result using complex q factor
        q(q0::Complex, M::Matrix) = (M[1] * q0 + M[3]) / (M[2] * q0 + M[4])
        R(q::Complex) = real(1 / q)
        w(q::Complex, λ, n = 1) = sqrt(-λ / (π * n * imag(1 / q)))
        propagate_ABCD(d) = [1 d; 0 1]
        lensmaker_ABCD(f) = [1 0; -1/f 1]
        # Beam parameters
        λ = 1000e-9
        w0 = 1e-3
        M2 = 1
        zr = SCDI.rayleigh_range(λ, w0, M2)
        # Lens parameters
        R1 = 1
        R2 = 1
        lens_y_location = 0.1
        nl = 1.5
        f = SCDI.lensmakers_eq(R1, -R2, nl)
        # Stuff
        dy = 0.001
        ys = 0:dy:1.5
        w_analytical = Vector{Float64}(undef, length(ys))
        R_analytical = Vector{Float64}(undef, length(ys))
        # Propagate using ABCD formalism
        q0 = 0 + zr * im
        for i in 1:length(ys)
            w_analytical[i] = w(q0, λ)
            R_analytical[i] = R(q0)
            # catch first lens
            if i * dy == lens_y_location
                q0 = q(q0, lensmaker_ABCD(f))
                continue
            end
            q0 = q(q0, propagate_ABCD(dy))
        end

        # Numerical result
        tl = SCDI.ThinLensSDF(R1, R2, 0.025)
        lens = SCDI.Lens(uuid4(), tl, x -> nl)
        system = SCDI.System(lens)
        SCDI.translate3d!(lens, [0, lens_y_location, 0])
        # Create and solve beam, calculate beam parameters
        gauss = SCDI.GaussianBeamlet(SCDI.Ray([0.0, 0, 0], [0.0, 1, 0]),
            λ,
            w0,
            support = [1, 0, 0],
            M2 = 1)
        SCDI.solve_system!(system, gauss)
        w_numerical, R_numerical, ψ_numerical, w0_numerical = SCDI.gauss_parameters(gauss,
            ys)
        # Compare beam radius to within 1 μm
        @test all(isapprox.(w_analytical, w_numerical, atol = 1e-6))
        # Compare if radius of curvature agreement to within 1 cm above 95% and any NaN
        temp = isapprox.(R_analytical, R_numerical, atol = 1e-2)
        @test sum(temp) / length(temp) > 0.95
        @test !any(isnan.(R_numerical))
        # Compare if Gouy phase zero at waists
        w0, i = findmin(w_analytical)
        @test isapprox(ψ_numerical[1], 0, atol = 1e-3)
        @test isapprox(ψ_numerical[i], 0, atol = 1e-3)
        # Compare calculated waist after lens
        @test isapprox(w0_numerical[i], w0, atol = 1e-7)

        @testset "Testing isparaxial and istilted" begin
            # Before lens rotation
            @test SCDI.istilted(system, gauss) == false
            @test SCDI.isparaxial(system, gauss) == true
            # Tilt lens, test again with 30° threshold for paraxial approx.
            SCDI.zrotate3d!(lens, deg2rad(45))
            SCDI.solve_system!(system, gauss)
            @test SCDI.istilted(system, gauss) == true
            @test SCDI.isparaxial(system, gauss, deg2rad(30)) == false
        end
    end
end

@testset "Interference" begin
    @testset "Pre-Beamsplitter tests with seperate beams" begin
        # Gauss beam parameters (selected for ring fringes)
        w0 = 0.01e-3
        λ = 1000e-9
        M2 = 1
        P0 = 1e-3
        I0 = 2 * P0 / (π * w0^2)
        E0 = SCDI.electric_field(I0)
        zR = SCDI.rayleigh_range(λ, w0, M2)
        # Detector parameters
        z = 0.1     # distance to detector
        l = 1e-2    # detector size
        n = 1000    # detector grid resolution
        # Lens parameters
        R1 = R2 = d = 0.01
        nl = 1.5
        f = SCDI.lensmakers_eq(R1, -R2, nl)
        # Raytracing system (for all tests)
        pd_l = SCDI.Photodetector(l, n)
        pd_s = SCDI.Photodetector(l/10, n÷10)
        ln = SCDI.ThinLens(R1, R2, d, nl)
        SCDI.translate3d!(pd_l, [0, z, 0])
        SCDI.translate3d!(pd_s, [0, z, 0])
        SCDI.translate3d!(ln, [0, z - f, 0])
        
        @testset "Testing fringe pattern" begin
            system = SCDI.System(pd_l)
            Δz = 5e-3   # arm length difference
            # Analytic solution
            xs = ys = LinRange(-l / 2, l / 2, n)
            screen = zeros(ComplexF64, length(xs), length(ys))
            for (j, y) in enumerate(ys)
                for (i, x) in enumerate(xs)
                    r = sqrt(x^2 + y^2)
                    screen[i, j] += SCDI.electric_field(r, z, E0, w0, λ, M2)
                    screen[i, j] += SCDI.electric_field(r, z + Δz, E0, w0, λ, M2)
                end
            end
            # Numerical solution
            SCDI.reset_photodetector!(pd_l)
            g_1 = SCDI.GaussianBeamlet(SCDI.Ray([0.0, 0, 0], [0.0, 1, 0]),
                λ,
                w0,
                M2 = M2,
                P0 = P0)
            g_2 = SCDI.GaussianBeamlet(SCDI.Ray([0.0, -Δz, 0], [0.0, 1, 0]),
                λ,
                w0,
                M2 = M2,
                P0 = P0)
            SCDI.solve_system!(system, g_1)
            SCDI.solve_system!(system, g_2)

            # Compare solutions
            I_analytical = SCDI.intensity.(screen)
            I_numerical = SCDI.intensity.(pd_l.field)
            Pt = SCDI.optical_power(pd_l)
            @test all(isapprox.(I_analytical, I_numerical, atol = 2e-1))
            @test isapprox(Pt, 2 * P0, atol = 3e-5)
        end

        @testset "Testing λ phase shift" begin
            system = SCDI.System([pd_s, ln])
            # Numerical solution
            Δz = LinRange(0, λ, 50)
            Pt_numerical = zeros(length(Δz))
            for (i, z_i) in enumerate(Δz)
                SCDI.reset_photodetector!(pd_s)
                g_1 = SCDI.GaussianBeamlet(SCDI.Ray([0.0, 0, 0], [0.0, 1, 0]),
                    λ,
                    w0,
                    M2 = M2,
                    P0 = P0)
                g_2 = SCDI.GaussianBeamlet(SCDI.Ray([0.0, z_i, 0], [0.0, 1, 0]),
                    λ,
                    w0,
                    M2 = M2,
                    P0 = P0)
                SCDI.solve_system!(system, g_1)
                SCDI.solve_system!(system, g_2)
                Pt_numerical[i] = SCDI.optical_power(pd_s)
                # Test length/opl function
                @test length(g_1) == z
                @test length(g_1) < length(g_1, opl=true)
                @test length(g_2) == length(g_1) - z_i
            end
            # Analytical solution (cosine over Δz), ref. power is 4*P0 since beam splitter is missing
            Pt_analytical = 4 * P0 * [(cos(2π * z / (maximum(Δz))) + 1) / 2 for z in Δz]

            # Compare detectors (this also tests correct behavior when focussing the beam)
            @test all(isapprox.(Pt_numerical, Pt_analytical, atol = 1e-4))
        end
    end

    @testset "Michelson Interferometer" begin
        # setup Michelson Interferometer
        l_0 = 0.1
        pd_size = SCDI.inch / 5
        pd_resolution = 100
        m1 = SCDI.RectangularPlanoMirror2D(SCDI.inch)
        m2 = SCDI.RectangularPlanoMirror2D(SCDI.inch)
        bs = SCDI.ThinBeamSplitter(SCDI.inch, 0.5)
        pd = SCDI.Photodetector(pd_size, pd_resolution)
        SCDI.translate3d!(m1, [l_0, 0, 0])
        SCDI.translate3d!(m2, [0, l_0, 0])
        SCDI.translate3d!(pd, [-l_0, 0, 0])
        SCDI.zrotate3d!(bs, deg2rad(45))
        SCDI.zrotate3d!(m1, deg2rad(90))
        SCDI.zrotate3d!(pd, deg2rad(90))

        system = SCDI.System([m1, m2, bs, pd])

        # Test correct values for reflectivity/transmission
        @test isvalid(bs)

        @testset "Equal armlength MI - integrated power" begin
            # setup 635 nm laser with 0.1 mm waist for fast divergence
            λ = 635e-9
            P_0 = 5e-3
            ray = SCDI.Ray([0, -l_0, 0], [0, 1.0, 0])
            beam = SCDI.GaussianBeamlet(ray, λ, 1e-4, P0 = P_0)

            # Shift mirror #2 by -λ to +λ
            lambdas = LinRange(-λ, λ, 200)

            path_length_numerical = zeros(length(lambdas))
            optical_pwr_numerical = zeros(length(lambdas))

            for (i, lambda) in enumerate(lambdas)
                SCDI.translate_to3d!(m2, [0, l_0, 0] + [0, lambda, 0])
                SCDI.reset_photodetector!(pd)
                SCDI.solve_system!(system, beam)

                # Moving mirror path length
                path_length_numerical[i] = length(beam.children[1].children[2])
                optical_pwr_numerical[i] = SCDI.optical_power(pd)
            end

            path_length_analytical = @. 2 * lambdas + 4l_0
            optical_pwr_analytical = @. P_0 * (1 / 2 * cos(2π * (2lambdas / λ) + π) + 1 / 2)

            # Compare correct PD signal and λ shift in moving arm
            @test all(isapprox.(optical_pwr_analytical, optical_pwr_numerical, atol = 5e-6))
            @test all(isapprox.(path_length_analytical, path_length_numerical))
        end

        @testset "Unequal armlength MI - electrical field" begin
            λ = 635e-9
            w0 = 1e-4
            P0 = 1e-3
            M2 = 1
            I0 = 2 * P0 / (π * w0^2)
            E0 = SCDI.electric_field(I0) * 1 / sqrt(2)^2
            zR = SCDI.rayleigh_range(λ, w0, M2)

            ray = SCDI.Ray([0, -l_0, 0], [0, 1.0, 0])
            beam = SCDI.GaussianBeamlet(ray, λ, w0, P0 = P0, M2 = M2)

            # arm length diff
            Δl = 1 * l_0
            SCDI.translate_to3d!(m2, [0, l_0 + Δl, 0])

            # numerical solution
            SCDI.reset_photodetector!(pd)
            SCDI.solve_system!(system, beam)

            # analytical solution
            short_arm = 4l_0
            long_arm = short_arm + 2Δl
            xs = ys = LinRange(-pd_size / 2, pd_size / 2, pd_resolution)
            screen = zeros(ComplexF64, length(xs), length(ys))
            for (j, y) in enumerate(ys)
                for (i, x) in enumerate(xs)
                    r = sqrt(x^2 + y^2)
                    screen[i, j] += SCDI.electric_field(r, short_arm, E0, w0, λ, M2)
                    screen[i, j] += SCDI.electric_field(r, long_arm, E0, w0, λ, M2) * exp(im*pi)
                end
            end

            Re_analytical = real.(screen)
            Im_analytical = imag.(screen)

            Re_numerical = real(pd.field)
            Im_numerical = imag(pd.field)

            # Compare solutions, units V/m
            @test all(isapprox.(Re_analytical, Re_numerical, atol = 5e-2))
            @test all(isapprox.(Im_analytical, Im_numerical, atol = 5e-2))
        end
    end

    @testset "Testing power conservation" begin
        # variables
        P0 = 0.5 # W
        l0 = 0.1 # m
        w0 = 0.5e-3
        λ = 1064e-9
        
        bs = SCDI.ThinBeamSplitter(10e-3);
        pd_1 = SCDI.Photodetector(10e-3, 100);
        pd_2 = SCDI.Photodetector(10e-3, 100);
        
        SCDI.zrotate3d!(bs, deg2rad(45))
        SCDI.translate3d!(pd_1, [0, l0, 0])
        SCDI.zrotate3d!(pd_1, deg2rad(180))
        
        SCDI.translate3d!(pd_2, [l0, 0, 0])
        SCDI.zrotate3d!(pd_2, deg2rad(90))
        
        # add BS and PD orientation error
        SCDI.zrotate3d!(bs, deg2rad(0.017))
        SCDI.zrotate3d!(pd_1, deg2rad(10))
        SCDI.xrotate3d!(pd_1, deg2rad(15))
        
        # define system and beams -> solve
        system = SCDI.System([bs, pd_1, pd_2]);
        
        phis = LinRange(0, 2pi, 25)
        p1 = similar(phis)
        p2 = similar(phis)
        
        for (i, phi) in enumerate(phis)
            l1 = SCDI.GaussianBeamlet(SCDI.Ray([0, -l0, 0], [0, 1., 0]), λ, w0; P0);
            l2 = SCDI.GaussianBeamlet(SCDI.Ray([-l0, 0, 0], [1., 0, 0]), λ, w0; P0);
            l1.E0 *= exp(im*phi)
            SCDI.reset_photodetector!(pd_1)
            SCDI.reset_photodetector!(pd_2)
            SCDI.solve_system!(system, l1)
            SCDI.solve_system!(system, l2)
            p1[i] = SCDI.optical_power(pd_1)
            p2[i] = SCDI.optical_power(pd_2)
            # Test power conservation
            @test p1[i] + p2[i] - 2P0 < 1e-4 # W
        end 
    end
end

@testset "Polarized rays" begin
    @testset "Polarization transforms" begin
        # Reflection matrix and lin. x-pol
        J = [-1 0 0; 0 1 0; 0 0 1]
        E0 = [1, 0, 0]
        
        @testset "90° reflection" begin
            in_dir = [0,0,1]
            out_dir = [1,0,0]
            @test SCDI._calculate_global_E0(in_dir, out_dir, J, E0) ≈ [0,0,-1]
        end

        @testset "0° reflection" begin
            in_dir = [0,0,1]
            out_dir = [0,0,-1]
            @test SCDI._calculate_global_E0(in_dir, out_dir, J, E0) ≈ [-1,0,0]
        end
    end

    @testset "Mirror reflections" begin
        # Setup system as in https://opg.optica.org/ao/fulltext.cfm?uri=ao-50-18-2855&id=218813
        m1 = SCDI.RectangularPlanoMirror2D(1.)
        m2 = SCDI.RectangularPlanoMirror2D(1.)
        m3 = SCDI.RectangularPlanoMirror2D(1.)
        SCDI.translate3d!(m2, [2,0,0])
        SCDI.translate3d!(m3, [2,2,0])
        SCDI.zrotate3d!(m1, deg2rad(-90))
        SCDI.yrotate3d!(m1, deg2rad(45))
        SCDI.zrotate3d!(m2, deg2rad(45))
        SCDI.xrotate3d!(m3, deg2rad(135))
        
        system = SCDI.StaticSystem([m1, m2, m3])
        
        I0_1 = 1
        I0_2 = 5
        lin_x_pol = [I0_1,0,0]
        lin_y_pol = [0,I0_2,0]
        
        # Beam of polarized rays
        ray = SCDI.PolarizedRay([0.,0,-2], [0,0,1], 1000e-9, lin_x_pol)
        beam = SCDI.Beam(ray)
        
        @testset "x-Polarization" begin
            SCDI.polarization!(ray, lin_x_pol)
            # test tracing
            SCDI.solve_system!(system, beam)
            @test SCDI.polarization(beam.rays[1]) ≈ lin_x_pol
            @test SCDI.polarization(beam.rays[2]) ≈ [0,0,-I0_1]
            @test SCDI.polarization(beam.rays[3]) ≈ [0,0, I0_1]
            @test SCDI.polarization(beam.rays[4]) ≈ [0,-I0_1,0]
            @test length(beam) == 6.0
        end
        
        @testset "y-Polarization" begin
            SCDI.polarization!(ray, lin_y_pol)
            SCDI.translate3d!(m3, [0,2,0])
            # test retracing
            SCDI.solve_system!(system, beam)
            @test SCDI.polarization(beam.rays[1]) ≈ lin_y_pol
            @test SCDI.polarization(beam.rays[2]) ≈ [0,-I0_2,0]
            @test SCDI.polarization(beam.rays[3]) ≈ [ I0_2,0,0]
            @test SCDI.polarization(beam.rays[4]) ≈ [-I0_2,0,0]
            @test length(beam) == 8.0
        end
    end

    @testset "Brewster windows" begin
        brewster_angle(n) = atan(n)
        # Testcase based on 5 successive Brewster windows
        n = 1.5
        θb = brewster_angle(n)
        d = 0.1
        # Calculate transmission efficiency
        rs, rp, ts, tp = SCDI.fresnel_coefficients(θb, n)
        Ts = 1 - abs2(rs)
        Tp = 1 - abs2(rp)
        # Setup testcase
        s1 = SCDI.CuboidMesh((1., d, 1.))
        s2 = SCDI.CuboidMesh((1., d, 1.))
        s3 = SCDI.CuboidMesh((1., d, 1.))
        s4 = SCDI.CuboidMesh((1., d, 1.))
        s5 = SCDI.CuboidMesh((1., d, 1.))
        l1 = SCDI.Lens(uuid4(), s1, x->n)
        l2 = SCDI.Lens(uuid4(), s2, x->n)
        l3 = SCDI.Lens(uuid4(), s3, x->n)
        l4 = SCDI.Lens(uuid4(), s4, x->n)
        l5 = SCDI.Lens(uuid4(), s5, x->n)
        SCDI.translate3d!.([l1, l2, l3, l4, l5], Ref([-0.5,-d/2,-0.5]))
        SCDI.set_new_origin3d!.(SCDI.shape.([l1, l2, l3, l4, l5]))
        SCDI.translate3d!(l2, [0,0.5, -1d/2])
        SCDI.translate3d!(l3, [0,1.0, -2d/2])
        SCDI.translate3d!(l4, [0,1.5, -3d/2])
        SCDI.translate3d!(l5, [0,2.0, -4d/2])
        SCDI.xrotate3d!.([l1, l2, l3, l4, l5], -θb)
        # Solve system of s- and p-polarized beams
        system = SCDI.StaticSystem([l1, l2, l3, l4, l5])
        x_pol_ray = SCDI.PolarizedRay([-0.1, -1, 0], [0, 1., 0], 1000e-9, [SCDI.electric_field(1), 0, 0])
        z_pol_ray = SCDI.PolarizedRay([+0.1, -1, 0], [0, 1., 0], 1000e-9, [0, 0, SCDI.electric_field(1)])
        s_beam = SCDI.Beam(x_pol_ray)
        p_beam = SCDI.Beam(z_pol_ray)
        SCDI.solve_system!(system, s_beam)
        SCDI.solve_system!(system, p_beam)
        # Since system is non-focussing, calculate pseudo-intensity
        pseudo_Is = abs2(SCDI.polarization(last(SCDI.rays(s_beam)))[1]) / (2*SCDI.Z_vacuum)
        pseudo_Ip = abs2(SCDI.polarization(last(SCDI.rays(p_beam)))[3]) / (2*SCDI.Z_vacuum)
        # Test against m interfaces
        m = length(system.objects) * 2
        @test pseudo_Is ≈ Ts^m
        @test pseudo_Ip ≈ Tp^m
    end

    @testset "Fresnel rhomb" begin
        # Create Fresnel rhomb with n=1.5 and θ=53.3° for quarter-wave plate effect
        n = 1.5
        s1 = SCDI.CuboidMesh((0.5,1.25,0.5), deg2rad(53.3))
        l1 = SCDI.Lens(uuid4(), s1, x->n)
        SCDI.translate3d!(l1, [-0.25, 0, -0.25])
        SCDI.set_new_origin3d!(s1)
        # Rotate prism to obtain 45° beam input polarization
        SCDI.yrotate3d!(l1, deg2rad(135))
        # Solve system
        system = SCDI.StaticSystem([l1])
        ray = SCDI.PolarizedRay([0, -1, 0], [0, 1., 0], 1000e-9, [0, 0, SCDI.electric_field(1)])
        beam = SCDI.Beam(ray)
        SCDI.solve_system!(system, beam)
        # Assumes propagation along the y-axis after rhomb, calculate polarization state
        Ex = getindex.(SCDI.polarization.(beam.rays), 1)
        Ey = getindex.(SCDI.polarization.(beam.rays), 2)
        Ez = getindex.(SCDI.polarization.(beam.rays), 3)
        # Test for circular polarization and Ey error
        phi = angle(last(Ez)) - angle(last(Ex))
        @test phi ≈ π/2
        @test abs(last(Ey)) < 2e-14
    end

    @testset "Mach-Zehnder Interferometer" begin
        # setup MZI
        m1 = SCDI.RectangularPlanoMirror2D(SCDI.inch)
        m2 = SCDI.RectangularPlanoMirror2D(SCDI.inch)
        b1 = SCDI.ThinBeamSplitter(SCDI.inch, 0.5)
        b2 = SCDI.ThinBeamSplitter(SCDI.inch, 0.5)
        
        system = SCDI.StaticSystem([m1, m2, b1, b2])
        
        SCDI.translate3d!(b1, [0*SCDI.inch, 0*SCDI.inch, 0])
        SCDI.translate3d!(b2, [2*SCDI.inch, 2*SCDI.inch, 0])
        SCDI.translate3d!(m1, [0*SCDI.inch, 2*SCDI.inch, 0])
        SCDI.translate3d!(m2, [2*SCDI.inch, 0*SCDI.inch, 0])
        
        # Rotate with consideration to mirror/bs normal
        SCDI.zrotate3d!(b1, deg2rad(360-135))
        SCDI.zrotate3d!(b2, deg2rad(45))
        SCDI.zrotate3d!(m1, deg2rad(360-135))
        SCDI.zrotate3d!(m2, deg2rad(45))
    
        ray = SCDI.PolarizedRay([0, -0.1, 0], [0., 1., 0], 1000e-9, [0, 0, 1])
        beam = SCDI.Beam(ray)
        
        @testset "z-polarized ray along y-axis" begin
            # Solve with z-polarized ray along y-axis
            SCDI.polarization!(ray, [0, 0, 1])
            SCDI.solve_system!(system, beam)
    
            # Extract E0s: t - transmitted, r - reflected
            t = beam.children[1].rays[1].E0
            r = beam.children[2].rays[1].E0
            tr = beam.children[1].rays[2].E0
            rr = beam.children[2].rays[2].E0
            trt = beam.children[1].children[1].rays[1].E0
            trr = beam.children[1].children[2].rays[1].E0
            rrt = beam.children[2].children[1].rays[1].E0
            rrr = beam.children[2].children[2].rays[1].E0
    
            # Test phase flips
            @test t[3] ≈ sqrt(2)/2
            @test r[3] ≈ -sqrt(2)/2
            @test tr[3] ≈ -sqrt(2)/2
            @test rr[3] ≈ sqrt(2)/2
            @test trt ≈ rrr
            @test trr ≈ rrt
        end
    
        @testset "x-polarized ray along y-axis" begin
            # Test num. of leaves before retracing
            @test length(collect(Leaves(beam))) == 4
            # Retrace with z-polarized ray along y-axis
            SCDI.polarization!(ray, [1, 0, 0])
            SCDI.solve_system!(system, beam)
    
            # Extract E0s: t - transmitted, r - reflected
            t = beam.children[1].rays[1].E0
            r = beam.children[2].rays[1].E0
            tr = beam.children[1].rays[2].E0
            rr = beam.children[2].rays[2].E0
            trt = beam.children[1].children[1].rays[1].E0
            trr = beam.children[1].children[2].rays[1].E0
            rrt = beam.children[2].children[1].rays[1].E0
            rrr = beam.children[2].children[2].rays[1].E0
    
            # Test phase flips
            @test t[1] ≈ sqrt(2)/2
            @test r[2] ≈ -sqrt(2)/2
            @test tr ≈ r
            @test rr ≈ t
            @test trt ≈ rrr
            @test trr ≈ rrt
        end
    end
end