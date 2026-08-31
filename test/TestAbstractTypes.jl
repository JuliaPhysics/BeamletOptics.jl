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
    @test isdefined(BMO, :AbstractSegment)
    @test isdefined(BMO, :OpenSegment)
    @test isdefined(BMO, :LineSegment)
    @test isdefined(BMO, :AbstractRay)
    @test isdefined(BMO, :AbstractBeam)
    @test isdefined(BMO, :AbstractSystem)
    @test isdefined(BMO, :AbstractIntersection)
    @test isdefined(BMO, :Intersection)
    @test isdefined(BMO, :MultiIntersection)
    @test isdefined(BMO, :Hint)
    @test isdefined(BMO, :AbstractInteraction)
    @test isdefined(BMO, :AbstractVolumeModel)
    @test isdefined(BMO, :VolumeInteraction)

    # Generate test structs
    struct TestSystem <: BMO.AbstractSystem end

    @testset "AbstractSystem" begin
        # no tests
        sys = TestSystem()
    end

    mutable struct TestRay{T, S <: BMO.AbstractSegment{T}} <: BMO.AbstractRay{T, S}
        segment::S
        λ::T
        n::T
    end

    TestRay(pos::AbstractArray{T}, dir::AbstractArray{T}) where {T} = TestRay(
        BMO.OpenSegment(pos, dir, zero(T)), T(1000e-9), T(1))

    Base.length(::TestRay) = π

    @testset "AbstractRay" begin
        r = TestRay([0.0, 0, 0], [1.0, 0, 0])
        # Test getters
        @test position(r) == [0.0, 0, 0]
        @test BMO.direction(r) == [1.0, 0, 0]
        @test BMO.wavelength(r) == r.λ
        @test BMO.refractive_index(r) == r.n
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
        @test is_1 isa BMO.Intersection{Float64}
        @test length(is_1) == 2
        @test length(is_2) == 1
        @test isnothing(is_3)
    end

    mutable struct TestBeam{T, R <: BMO.AbstractRay{T}} <: BMO.AbstractBeam{T, R}
        parent::BMO.Nullable{TestBeam{T, R}}
        children::Vector{TestBeam{T, R}}
    end

    TestBeam() = TestBeam{Float64, TestRay{Float64, BMO.OpenSegment{Float64}}}(nothing, Vector{TestBeam{Float64, TestRay{Float64, BMO.OpenSegment{Float64}}}}())

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
        # Adding another child to a beam that already has one must fail
        cbr = TestBeam()
        @test_throws "Adding child to beam failed" BMO.children!(
            cb2,
            cbr)
        # Test child removal
        BMO._drop_beams!(cb2)
        @test isempty(BMO.children(cb2))
        # Stuff
        @test_throws "_reset_beam! not implemented for $(typeof(cb2))" BMO._reset_beam!(cb2)
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
        # normal3d has no generic AbstractShape implementation
        @test_throws "normal3d not defined for" BMO.normal3d(shape)
        @test_throws "normal3d not defined for" BMO.normal3d(shape, 1)

        @testset "Testing AbstractRay - AbstractShape" begin
            shape = TestShapeless([1, 0, 0], Matrix{Int}(I, 3, 3))
            ray = TestRay([0.0, 0, 0], [1.0, 0, 0])
            @test BMO.isinfrontof(shape, ray) == true
            ray = TestRay([0.0, 0, 0], -[1.0, 0, 0])
            @test BMO.isinfrontof(shape, ray) == false
            ray = TestRay([0.0, 0, 0], [0.0, 1, 0])
            @test BMO.isinfrontof(shape, ray) == false
            ray = TestRay([0.0, 0, 0], [1.0, 1, 0])
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
        p = [0.0, 2.0, 0.0]
        n = [0.0, 0.0, 1.0]

        @testset "Intersection(t, p, n)" begin
            int = Intersection(t, p, n)
            @test int isa BMO.AbstractIntersection{Float64}
            @test length(int) == t
            @test position(int) == Point3(p)
            @test normal3d(int) == Point3(n)
        end

        @testset "MultiIntersection" begin
            hit = Intersection(t, p, n)
            other_object = TestObject(TestShapeless())

            mi = BMO.MultiIntersection(hit; exiting = other_object, entering = object)
            @test mi isa BMO.AbstractIntersection{Float64}
            @test length(mi) == t
            @test position(mi) == Point3(p)
            @test normal3d(mi) == Point3(n)
            @test BMO.exiting(mi) === other_object
            @test BMO.entering(mi) === object
        end
    end

        @testset "AbstractMedium" begin
            amb = Ambient()
            @test BMO.refractive_index(amb, 532e-9) == 1.0
            @test BMO.complex_refractive_index(amb, 532e-9) == 1.0 + 0.0im
            @test BMO.extinction_coefficient(amb, 532e-9) == 0.0
            @test BMO.absorption_coefficient(amb, 532e-9) == 0.0

            glass = IsotropicMedium(:N_BK7, λ -> 1.5168)
            @test BMO.refractive_index(glass, 532e-9) ≈ 1.5168
            @test BMO.complex_refractive_index(glass, 532e-9) ≈ 1.5168 + 0.0im

            const_medium = IsotropicMedium(:Fixed, 1.45)
            @test BMO.refractive_index(const_medium) == 1.45

            # Complex refractive index (e.g. absorbing metal or gain medium)
            metal = IsotropicMedium(:Gold, 0.2 + 3.0im)
            @test BMO.refractive_index(metal) == 0.2
            @test BMO.complex_refractive_index(metal) == 0.2 + 3.0im
            @test BMO.extinction_coefficient(metal) == 3.0

            # Attenuation / linear absorption coefficient
            abs_glass = IsotropicMedium(:AbsGlass, 1.5, 100.0) # α = 100 1/m
            @test BMO.refractive_index(abs_glass, 500e-9) == 1.5
            @test BMO.absorption_coefficient(abs_glass, 500e-9) ≈ 100.0
            @test isapprox(BMO.extinction_coefficient(abs_glass, 500e-9), 100.0 * 500e-9 / (4 * π), atol=1e-10)

            # Temperature and pressure dependent thermo-optic response
            air_tp = IsotropicMedium(:AirTP, (λ; T=293.15, P=101325.0) -> 1.0 + 2.7e-4 * (P / 101325.0) * (293.15 / T))
            @test BMO.refractive_index(air_tp, 500e-9) ≈ 1.00027
            @test isapprox(BMO.refractive_index(air_tp, 500e-9; T=300.0), 1.0 + 2.7e-4 * (293.15 / 300.0), atol=1e-8)
            @test isapprox(BMO.refractive_index(air_tp, 500e-9; P=202650.0), 1.0 + 2 * 2.7e-4, atol=1e-8)
        end

        @testset "AbstractSurfaceModel" begin
            @test FresnelSurface() isa AbstractSurfaceModel
            @test IdealMirrorSurface() isa AbstractSurfaceModel
            @test AbsorbingSurface() isa AbstractSurfaceModel
            @test DetectorSurface() isa AbstractSurfaceModel
            @test CoatedSurface(:AR) isa AbstractSurfaceModel
            @test GratingSurface(600.0, 1) isa AbstractSurfaceModel
        end

        @testset "AbstractVolumeModel" begin
            @test VolumeInteraction() isa AbstractVolumeModel
            struct TestVolumeModel <: AbstractVolumeModel end
            @test TestVolumeModel() isa AbstractVolumeModel
            vm = VolumeInteraction()
            med = IsotropicMedium(:TestMed, 1.5)
            r = Ray([0.0, 0, 0], [0.0, 1, 0], 532e-9, 1.0)
            @test_throws ErrorException interact3d(vm, med, r)
        end

        @testset "Coating" begin
            shape = TestShapeless()
            c = Coating(shape, FresnelSurface())
            @test c isa AbstractCoating
            @test BMO.shape(c) === shape
            @test BMO.surface_model(c) === FresnelSurface()
        end

        @testset "Transition" begin
            trans = resolve_transition(IsotropicMedium(:Glass, 1.5), Ambient(), Ray([0.0, 0, 0], [0.0, 1, 0]), Point3(0.0, -1, 0))
            @test trans isa Transition
            @test trans.is_entering == true

            ray = Ray([0.0, 0, 0], [0.0, 1, 0], 532e-9, 1.0)
            int = Intersection(1.0, [0.0, 1.0, 0.0], [0.0, -1.0, 0.0])
            out_ray = interact3d(FresnelSurface(), trans, int, ray)
            @test out_ray isa Ray
            @test isapprox(norm(BMO.direction(out_ray)), 1.0)

            mirror_ray = interact3d(IdealMirrorSurface(), trans, int, ray)
            @test mirror_ray isa Ray
            @test BMO.direction(mirror_ray)[2] ≈ -1.0

            absorb_ray = interact3d(AbsorbingSurface(), trans, int, ray)
            @test isnothing(absorb_ray)

            # Grating diffraction test: 600 lines/mm, grooves along z
            g = GratingSurface(600e3, 1, Point3(0.0, 0.0, 1.0))
            g_ray = interact3d(g, trans, int, Ray([0.0, 0, 0], [0.0, 1, 0], 500e-9, 1.0))
            @test g_ray isa Ray
            # sin(θm) = 500e-9 * 600e3 = 0.30 -> dx = 0.30, dy = -sqrt(1 - 0.3^2)
            @test isapprox(BMO.direction(g_ray)[1], 0.30, atol=1e-6)
            @test isapprox(BMO.direction(g_ray)[2], -sqrt(1 - 0.3^2), atol=1e-6)

            # Auto groove axis grating
            g_auto = GratingSurface(600e3, 1)
            g_auto_ray = interact3d(g_auto, trans, int, Ray([0.0, 0, 0], [0.0, 1, 0], 500e-9, 1.0))
            @test g_auto_ray isa Ray
            @test isapprox(norm(BMO.direction(g_auto_ray)), 1.0)

            # CoatedSurface tests
            cs_mirror = CoatedSurface(IdealMirrorSurface())
            @test interact3d(cs_mirror, trans, int, ray) isa Ray
            @test BMO.direction(interact3d(cs_mirror, trans, int, ray))[2] ≈ -1.0

            cs_func = CoatedSurface((t, i, r) -> Ray(BMO.position(i), Point3(1.0, 0.0, 0.0), BMO.wavelength(r), 1.0))
            custom_ray = interact3d(cs_func, trans, int, ray)
            @test custom_ray isa Ray
            @test BMO.direction(custom_ray) == Point3(1.0, 0.0, 0.0)

            # PolarizedRay on IdealMirrorSurface
            pol_ray = PolarizedRay([0.0, 0, 0], [0.0, 1, 0], 532e-9, [1.0, 0.0, 0.0])
            pol_mirror_ray = interact3d(IdealMirrorSurface(), trans, int, pol_ray)
            @test pol_mirror_ray isa PolarizedRay
            @test BMO.direction(pol_mirror_ray)[2] ≈ -1.0

            # surface_model fallback tests
            @test surface_model(nothing) isa FresnelSurface
            @test surface_model(IdealMirrorSurface()) isa IdealMirrorSurface

            # Standalone Coating object in a System
            coat_shape = deepcopy(BMO.shape(BMO.ThinBeamsplitter(0.05)))
            translate3d!(coat_shape, [0.0, 0.05, 0.0])
            coat_obj = Coating(coat_shape, IdealMirrorSurface())
            coat_sys = System([coat_obj])
            coat_beam = Beam([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 532e-9)
            solve_system!(coat_sys, coat_beam)
            @test length(BMO.rays(coat_beam)) == 2
            @test BMO.direction(last(BMO.rays(coat_beam)))[2] ≈ -1.0
        end
    end

end # MODULE