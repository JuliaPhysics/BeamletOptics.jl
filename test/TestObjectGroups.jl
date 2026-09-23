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

    @testset "reset_rotation3d! at and near θ = π" begin
        for θ in (π, π - 1e-9, π / 2)
            rotate3d!(objects, [0, 1, 1], θ)
            reset_rotation3d!(objects)
            @test orientation(objects) == Matrix{Float64}(I, 3, 3)
            # center is a Float32 object
            @test orientation(center) ≈ I atol = sqrt(eps(Float32))
            for (i, obj) in enumerate(Leaves(BMO.objects(circle)))
                @test position(obj) ≈ [xs[i], ys[i], 0] atol = 1e-12
                @test orientation(obj) ≈ I atol = 1e-12
            end
        end
    end

    @testset "set_pivot3d!" begin
        # Fresh group built the same way as `objects`/`circle` above
        center2 = GroupTestObject(zeros(Float32, 3))
        circle2 = ObjectGroup([GroupTestObject([xs[i], ys[i], 0]) for i in eachindex(xs)])
        group2 = ObjectGroup([center2, circle2])

        p_before = [position(o) for o in BMO.objects(circle2)]
        o_before = [orientation(o) for o in BMO.objects(circle2)]
        center_p_before = position(center2)
        center_o_before = orientation(center2)

        p = [2.0, -1.0, 0.5]
        @test set_pivot3d!(group2, p) === nothing
        @test position(group2) == p
        @test orientation(group2) == Matrix{Float64}(I, 3, 3)
        # members are untouched (exact)
        for (i, o) in enumerate(BMO.objects(circle2))
            @test position(o) == p_before[i]
            @test orientation(o) == o_before[i]
        end
        @test position(center2) == center_p_before
        @test orientation(center2) == center_o_before

        # rotate3d! now rotates every member about the new pivot p
        Rp = BMO.rotate3d(normalize([1.0, 2.0, 3.0]), 0.7)
        rotate3d!(group2, Rp)
        for (i, o) in enumerate(BMO.objects(circle2))
            @test norm(position(o) - (p + Rp * (p_before[i] - p))) < 1e-12
        end
        # center2 is a Float32 object, use a looser tolerance
        @test norm(position(center2) - (p + Rp * (center_p_before - p))) < 1e-6

        # reset_translation3d! moves p back to the origin, members shifted by -p
        reset_translation3d!(group2)
        @test norm(position(group2)) < 1e-12
        for (i, o) in enumerate(BMO.objects(circle2))
            @test norm(position(o) - Rp * (p_before[i] - p)) < 1e-12
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

@testset "rotate3d! about an arbitrary pivot" begin
    # rotate3d!(x, R, pivot) and rotate3d!(x, axis, θ, pivot) for all shape and object types
    axis = [1.0, -2.0, 0.5]
    θ = 0.8
    pivot = [0.3, -0.7, 1.1]
    Rp = BMO.rotate3d(axis, θ)
    subjects = Any[
        BMO.CylinderSDF(1e-3, 2e-3),                                # shape
        BMO.CylinderSDF(1e-3, 2e-3) + BMO.CylinderSDF(2e-3, 1e-3),  # composite SDF
        BMO.CubeMesh(1e-3),                                         # mesh
        RoundPlanoMirror(10e-3, 2e-3),                              # object (SingleShape)
        CubeBeamsplitter(5e-3, λ -> 1.5),                           # object (MultiShape)
        ObjectGroup([RoundPlanoMirror(10e-3, 2e-3), CubeBeamsplitter(5e-3, λ -> 1.5)]),
    ]
    for x in subjects
        translate3d!(x, [1e-3, 2e-3, 3e-3])
        p0, O0 = position(x), orientation(x)
        a = deepcopy(x)
        b = deepcopy(x)
        @test rotate3d!(a, Rp, pivot) === nothing
        @test rotate3d!(b, 3.7 * axis, θ, Point3(pivot...)) === nothing
        @test position(a) ≈ pivot + Rp * (p0 - pivot) atol = 1e-12
        @test orientation(a) ≈ Rp * O0 atol = 1e-12
        @test position(b) ≈ position(a) atol = 1e-12
        @test orientation(b) ≈ orientation(a) atol = 1e-12
    end
end

end # MODULE