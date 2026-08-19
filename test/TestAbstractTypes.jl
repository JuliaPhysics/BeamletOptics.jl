module TestAbstractTypes

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics
using AbstractTrees

const BMO = BeamletOptics

@testset "Types" begin
    @test isdefined(BMO, :AbstractShape)
    @test isdefined(BMO, :AbstractObject)
    @test isdefined(BMO, :AbstractObjectGroup)
    @test isdefined(BMO, :AbstractRay)
    @test isdefined(BMO, :AbstractBeam)
    @test isdefined(BMO, :AbstractSystem)
    @test isdefined(BMO, :Intersection)
    @test isdefined(BMO, :AbstractIntersection)
    @test isdefined(BMO, :ShapeIntersection)
    @test isdefined(BMO, :ObjectIntersection)
    @test isdefined(BMO, :Hint)
    @test isdefined(BMO, :AbstractInteraction)

    # Generate test structs
    struct TestSystem <: BMO.AbstractSystem end

    @testset "AbstractSystem" begin
        # no tests
        sys = TestSystem()
    end

    mutable struct TestRay{T} <: BMO.AbstractRay{T}
        pos::Vector{T}
        dir::Vector{T}
        λ::T
        n::T
    end

    TestRay(pos::AbstractArray{T}, dir::AbstractArray{T}) where {T} = TestRay(
        pos, dir, T(1000e-9), T(1))

    Base.length(::TestRay) = π

    @testset "AbstractRay" begin
        r = TestRay([0.0, 0, 0], [1.0, 0, 0])
        # Test getters
        @test position(r) == r.pos
        @test BMO.direction(r) == r.dir
        @test BMO.wavelength(r) == r.λ
        @test BMO.refractive_index(r) == r.n
        # Test setters
        n_pos = [1, 1, 1]
        n_dir = [1.0, 1, 0]
        n_lam = 532e-9
        n_rfi = 1.5
        BMO.position!(r, n_pos)
        BMO.direction!(r, n_dir)
        BMO.wavelength!(r, n_lam)
        BMO.refractive_index!(r, n_rfi)
        @test position(r) == n_pos
        @test BMO.direction(r) ≈ n_dir .* (sqrt(2) / 2)
        @test BMO.wavelength(r) == n_lam
        @test BMO.refractive_index(r) == n_rfi
        @test length(r) == π
        # Test ray-plane intersection
        plane_pos = [1, 0, -1]
        plane_nml_1 = [-1, 0, 1]
        plane_nml_2 = [1, 0, 0]
        plane_nml_3 = [0, 0, 1]
        ray = TestRay([0.0, 0, 0], [1.0, 0, 0])
        is_1 = BMO.intersect3d(plane_pos, plane_nml_1, ray)
        is_2 = BMO.intersect3d(plane_pos, plane_nml_2, ray)
        is_3 = BMO.intersect3d(plane_pos, plane_nml_3, ray)
        @test length(is_1) == 2
        @test length(is_2) == 1
        @test isnothing(is_3)
    end

    mutable struct TestBeam{T} <: BMO.AbstractBeam{T, TestRay{T}}
        parent::BMO.Nullable{TestBeam}
        children::Vector{TestBeam}
    end

    TestBeam() = TestBeam{Float64}(nothing, Vector{TestBeam{Float64}}())

    @testset "AbstractBeam" begin
        # Create beam tree
        root = TestBeam()
        cb1 = TestBeam()
        cb2 = TestBeam()
        group = [cb1, cb2]
        cb3 = TestBeam()
        # Add children to root, cb2
        BMO.children!(root, group)
        BMO.children!(cb2, cb3)
        # Test tree structure
        @test treeheight(root) == 2
        @test treebreadth(root) == 2
        @test BMO.children(root) == group
        @test first(BMO.children(cb2)) === cb3
        # Test parent connection
        @test AbstractTrees.parent(root) === nothing
        @test AbstractTrees.parent(cb1) === root
        @test AbstractTrees.parent(cb2) === root
        @test AbstractTrees.parent(cb3) === cb2
        # Replace bottom child
        cbr = TestBeam()
        @test_throws "_modify_beam_head not implemented for $(typeof(cb2))" BMO.children!(
            cb2,
            cbr)
        # Test child removal
        BMO._drop_beams!(cb2)
        @test isempty(BMO.children(cb2))
        # Stuff
        @test_throws "_last_beam_intersection not implemented for $(typeof(cb2))" BMO._last_beam_intersection(cb2)
    end

    mutable struct TestShapeless{T} <: BMO.AbstractShape{T}
        pos::Vector{T}
        dir::Matrix{T}
    end

    TestShapeless() = TestShapeless{Float64}(zeros(3), Matrix{Float64}(I, 3, 3))

    @testset "AbstractShape" begin
        pos = zeros(3)
        dir = Matrix{Float64}(I, 3, 3)
        shape = TestShapeless(pos, dir)
        # Test get/set
        @test BMO.position(shape) == pos
        @test orientation(shape) == dir
        n_pos = [1, 1, 1]
        n_dir = BMO.rotate3d([0, 0, 1], π / 4)
        BMO.position!(shape, n_pos)
        BMO.orientation!(shape, n_dir)
        @test position(shape) == n_pos
        @test orientation(shape) == n_dir
        # Test translation
        translate3d!(shape, n_pos)
        @test position(shape) == 2 * n_pos
        reset_translation3d!(shape)
        @test position(shape) == zeros(3)
        # Test rotation for counter-clockwise in right-hand coord. system
        dir = Matrix{Float64}(I, 3, 3)
        BMO.orientation!(shape, dir)
        rotate3d!(shape, [0, 0, 1], deg2rad(45))
        @test all(orientation(shape)[[1, 2, 5]] .≈ sqrt(2) / 2)
        @test orientation(shape)[4] ≈ -sqrt(2) / 2
        rotate3d!(shape, [0, 0, 1], deg2rad(135))
        @test orientation(shape)[1:4:9] == [-1, -1, 1]
        xrotate3d!(shape, π)
        yrotate3d!(shape, π)
        zrotate3d!(shape, π)
        reset_rotation3d!(shape)
        @test orientation(shape) == dir
        # Test align3d
        target_vec = normalize([1, 1, 1])
        align3d!(shape, target_vec)
        @test orientation(shape)[:, 2] ≈ target_vec
        reset_rotation3d!(shape)
        # The following test are expected to do nothing but not throw exceptions
        ray = TestRay([0.0, 0, 0], [1.0, 0, 0])
        @test_logs (:warn, "No intersect3d method defined for:") BMO.intersect3d(
            shape,
            ray)

        @testset "Testing AbstractRay - AbstractShape" begin
            shape = TestShapeless([1, 0, 0], Matrix{Int}(I, 3, 3))
            ray = TestRay([0.0, 0, 0], [1.0, 0, 0])
            @test BMO.isinfrontof(shape, ray) == true
            BMO.direction!(ray, -[1, 0, 0])
            @test BMO.isinfrontof(shape, ray) == false
            BMO.direction!(ray, [0, 1, 0])
            @test BMO.isinfrontof(shape, ray) == false
            BMO.direction!(ray, [1.0, 1, 0])
            @test BMO.isinfrontof(shape, ray) == true
            @test norm(BMO.direction(ray)) ≈ 1
        end
    end

    struct TestObject{T, S <: BMO.AbstractShape{T}} <: BMO.AbstractObject{T}
        shape::S
    end

    TestObject() = TestObject(TestShapeless())

    @testset "AbstractObject" begin
        object = TestObject()
        @test isa(BMO.shape(object), TestShapeless)
        # Test forwarding of kin. API to object shape
        @test position(object) ==
              position(BMO.shape(object))
        @test position(object) ==
              position(BMO.shape(object))
        translate3d!(object, ones(3))
        rotate3d!(object, [0, 0, 1], π)
        @test position(object) == ones(3)
        @test orientation(object)[1:4:9] == [-1, -1, 1]
        reset_translation3d!(object)
        reset_rotation3d!(object)
        @test position(object) == zeros(3)
        @test orientation(object)[1:4:9] == ones(3)
        # Test translate_to3d
        target_pos = [1, 3, 9]
        translate_to3d!(object, target_pos)
        @test position(object) == target_pos

        @testset "Testing interact3d" begin
            sys = TestSystem()
            obj = TestObject()
            ray = TestRay(zeros(3), ones(3))
            beam = TestBeam()
            @test_logs (:warn, "No interact3d method defined for:") BMO.interact3d(
                sys,
                obj,
                beam,
                ray)===nothing
        end
    end

    @testset "AbstractIntersection" begin
        shape = TestShapeless()
        object = TestObject(shape)
        t = 2.0
        n = [0.0, 0.0, 1.0]

        @testset "ShapeIntersection" begin
            si = BMO.ShapeIntersection(shape, t, Point3(n))
            @test si isa BMO.AbstractIntersection{Float64}
            @test BMO.shape(si) === shape
            @test length(si) == t
            @test BMO.normal3d(si) == n

            # Test convenience constructor
            si2 = BMO.ShapeIntersection(t, n, shape)
            @test BMO.shape(si2) === shape
            @test length(si2) == t
            @test BMO.normal3d(si2) == n
        end

        @testset "ObjectIntersection" begin
            si = BMO.ShapeIntersection(shape, t, Point3(n))
            oi = BMO.ObjectIntersection(object, si)
            @test oi isa BMO.AbstractIntersection{Float64}
            # Test forwarding to the underlying ShapeIntersection
            @test BMO.object(oi) === object
            @test BMO.shape(oi) === shape
            @test length(oi) == t
            @test BMO.normal3d(oi) == n
        end
    end
end

end # MODULE