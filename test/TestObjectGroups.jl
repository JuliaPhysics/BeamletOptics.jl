module TestObjectGroups

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics
using AbstractTrees

const BMO = BeamletOptics

@testset "Object groups" begin
    mutable struct TestPoint{T} <: BMO.AbstractShape{T}
        pos::Point3{T}
        dir::Matrix{T}
    end

    TestPoint(position::AbstractArray{T}) where {T <: Real} = TestPoint{T}(
        Point3{T}(position),
        Matrix{T}(I, 3, 3))

    struct GroupTestObject{T <: Real, S <: BMO.AbstractShape{T}} <: BMO.AbstractObject{T}
        shape::S
    end

    GroupTestObject(position::AbstractArray) = GroupTestObject(TestPoint(position))

    n = 8
    xs = [cos(x) for x in LinRange(0, 2pi * (n - 1) / n, n)]
    ys = [sin(x) for x in LinRange(0, 2pi * (n - 1) / n, n)]

    # Test center with Float32, rest with Float64
    center = GroupTestObject(zeros(Float32, 3))
    circle = ObjectGroup([GroupTestObject([xs[i], ys[i], 0]) for i in eachindex(xs)])

    objects = ObjectGroup([center, circle])

    # Translation test
    target = [3, 0, 0]
    translate_to3d!(objects, target)

    @testset "translate3d" begin
        # Test if all objects/subgroups have been translated
        @test position(objects) == target
        @test position(center) == target
        @test position(circle) == target
        for (i, obj) in enumerate(BMO.objects(circle))
            @test position(obj) == [xs[i], ys[i], 0] + target
        end
    end

    # Rotation test
    angle = 2π / n
    rotate3d!(objects, [0, 0, 1], angle)

    @testset "rotate3d" begin
        # Test if all objects/subgroups have been rotated relative to the origin
        Rt = BMO.rotate3d([0, 0, 1], angle)
        xt = circshift(xs, -1)
        yt = circshift(ys, -1)
        @test orientation(objects) == Rt
        @test orientation(center) ≈ Rt
        @test orientation(circle) == Rt
        for (i, obj) in enumerate(BMO.objects(circle))
            @test orientation(obj) == Rt
            @test position(obj) ≈ [xt[i], yt[i], 0] + target
        end
    end

    # Reset test
    reset_translation3d!(objects)
    reset_rotation3d!(objects)

    @testset "reset functions" begin
        Ri = Matrix{Float64}(I, 3, 3)
        # Test if objects are reset correctly to initial positioning
        @test position(objects) == zeros(3)
        @test position(center) == zeros(3)
        @test position(circle) == zeros(3)
        @test orientation(objects) == Ri
        @test orientation(center) ≈ Ri
        @test orientation(circle) ≈ Ri
        for (i, obj) in enumerate(Leaves(BMO.objects(circle)))
            @test isapprox(position(obj)[1], xs[i], atol = 5e-16)
            @test isapprox(position(obj)[2], ys[i], atol = 5e-16)
        end
    end

    @testset "System compatibility" begin
        # Test if objects in ObjectGroup are exposed correctly when iterating
        system = System(objects)
        ctr = 0
        # Only the objects within the groups should be exposed
        for obj in BMO.objects(system)
            @test isa(obj, GroupTestObject)
            ctr += 1
        end
        @test ctr == n + 1
    end
end

end # MODULE