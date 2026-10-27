module TestObjectGroups

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics
using AbstractTrees

const BMO = BeamletOptics

"Static test object: opts out of the kinematic API via its trait."
struct StaticTestObject{T, S <: BMO.AbstractShape{T}} <: BMO.AbstractObject{T}
    shape::S
end
StaticTestObject() = StaticTestObject(BMO.CylinderSDF(1e-3, 2e-3))
BMO.kinematic_trait_of(::StaticTestObject) = BMO.Static()

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

    @testset "Type and members" begin
        @test !(BMO.AbstractObjectGroup <: BMO.AbstractObject)
        @test !(objects isa BMO.AbstractObject)
        @test objects isa BMO.ObjectOrGroup && center isa BMO.ObjectOrGroup
        @test !hasmethod(BMO.shape_trait_of, Tuple{typeof(objects)})
        @test !hasmethod(BMO.shape, Tuple{typeof(objects)})
        @test BMO.objects(objects) isa Vector{BMO.ObjectOrGroup}
        @test BMO.objects(objects)[1] === center && BMO.objects(objects)[2] === circle
        @test collect(AbstractTrees.children(objects)) == [center, circle]
        # only objects and groups are members, e.g. no shapes
        @test_throws TypeError ObjectGroup([center, BMO.CylinderSDF(1e-3, 2e-3)])
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

@testset "push!, delete!, pop! and popat! of an ObjectGroup" begin
    pose(x) = (copy(position(x)), copy(orientation(x)))
    same(xs, ys) = length(xs) == length(ys) && all(xs[i] === ys[i] for i in eachindex(xs))
    m1, m2, m3, m4 = (RoundPlanoMirror(10e-3, 2e-3) for _ in 1:4)
    foreach(((i, m),) -> translate3d!(m, [0, 0.05i, 0]), enumerate((m1, m2, m3, m4)))
    inner = ObjectGroup([m2])
    g = ObjectGroup([m1, inner])
    set_pivot3d!(g, [0.01, 0.02, 0.03])
    zrotate3d!(g, deg2rad(20))
    poses = IdDict{Any, Any}(o => pose(o) for o in (m1, m2, m3, m4, inner))
    g_pose = pose(g)

    # the constructor copies its argument
    v = [m3]
    h = ObjectGroup(v)
    push!(h, m4)
    @test length(v) == 1

    @test push!(g, m3, m4) === g
    @test same(BMO.objects(g), [m1, inner, m3, m4])
    @test pose(g) == g_pose
    @test all(pose(o) == poses[o] for o in (m1, m2, m3, m4, inner))
    @test delete!(g, m3) === g
    @test same(BMO.objects(g), [m1, inner, m4])
    @test pose(m3) == poses[m3]
    @test pose(g) == g_pose
    # nothing happens for an object that is not part of the group
    @test same(BMO.objects(delete!(g, m3)), [m1, inner, m4])

    # refused changes leave the group as it was
    before = copy(BMO.objects(g))
    @test_throws "already part of the group" push!(g, m1)
    @test_throws "already part of the group" push!(g, m2)  # within a subgroup
    @test_throws "already part of the group" push!(g, m3, m3)
    @test_throws "already part of the group" push!(g, g)
    @test_throws "already part of the group" push!(inner, g)  # a group can not hold itself
    @test_throws "already part of the group" push!(g, m3, ObjectGroup([m1]))
    @test_throws "remove it from its subgroup" delete!(g, m2)
    @test same(BMO.objects(g), before)
    @test same(BMO.objects(inner), [m2])

    # static and movable objects are not mixed
    @test_throws ArgumentError push!(g, StaticTestObject())
    @test same(BMO.objects(g), before)
    sg = ObjectGroup([])
    push!(sg, StaticTestObject())
    @test BMO.is_static(sg)
    @test_throws ArgumentError push!(sg, m3)
    @test length(BMO.objects(sg)) == 1

    # the last object can be removed
    e = ObjectGroup([m3])
    delete!(e, m3)
    @test isempty(BMO.objects(e))
    translate3d!(e, [1.0, 0, 0])
    @test position(e) == [1.0, 0, 0]
    @test pose(m3) == poses[m3]

    # pop! and popat! remove at the top level, a subgroup counts as one
    p = ObjectGroup([m3, inner, m4])
    p_pose = pose(p)
    @test pop!(p) === m4
    @test popat!(p, 2) === inner
    @test same(BMO.objects(p), [m3])
    @test same(BMO.objects(inner), [m2])
    @test_throws BoundsError popat!(p, 2)
    @test pop!(p) === m3
    @test_throws ArgumentError pop!(p)
    @test pose(p) == p_pose
    @test all(pose(o) == poses[o] for o in (m2, m3, m4, inner))

    # tracing follows the objects of a group in a system
    mirror = RoundPlanoMirror(10e-3, 2e-3)
    translate3d!(mirror, [0, 0.1, 0])
    mg = ObjectGroup([])
    system = System([mg])
    beam = Beam([0.0, 0, 0], [0.0, 1, 0], 1e-6)
    solve_system!(system, beam)
    @test isnothing(BMO.intersection(first(BMO.rays(beam))))
    push!(mg, mirror)
    solve_system!(system, beam)
    @test BMO.object(BMO.intersection(first(BMO.rays(beam)))) === mirror
    delete!(mg, mirror)
    solve_system!(system, beam)
    @test isnothing(BMO.intersection(first(BMO.rays(beam))))

    # an empty group or subgroup is no object of its system, also after its last object was removed
    a, b, c = (RoundPlanoMirror(10e-3, 2e-3) for _ in 1:3)
    sub = ObjectGroup([a])
    outer = ObjectGroup([sub, b])
    es = System([outer])
    delete!(sub, a)
    @test same(collect(BMO.objects(es)), [b])
    @test push!(es, a) === es
    @test same(collect(BMO.objects(es)), [b, a])
    @test pop!(outer) === b
    @test same(collect(BMO.objects(es)), [a])
    @test_throws "already an object of the system" push!(es, ObjectGroup([ObjectGroup([]), a]))
    push!(es, b)
    @test same(collect(BMO.objects(es)), [a, b])
    @test isempty(collect(BMO.objects(System([ObjectGroup([ObjectGroup([])])]))))
    @test same(collect(BMO.objects(push!(System([ObjectGroup([])]), c))), [c])

    # a StaticSystem keeps the objects that its groups had when it was built, without empty groups
    d = RoundPlanoMirror(10e-3, 2e-3)
    cg = ObjectGroup([c])
    ss = StaticSystem([ObjectGroup([]), cg])
    @test same(BMO.objects(ss), [c])
    push!(cg, d)
    @test same(BMO.objects(ss), [c])
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