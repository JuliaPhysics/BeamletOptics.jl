module TestBoundingSphere

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics
using Random

const BMO = BeamletOptics

# Wraps an sdf: counts its evaluations and, with a `radius`, has a bounding sphere around the local
# `center`, whose computations are counted too
mutable struct WrappedSDF{T, S <: BMO.AbstractSDF{T}, R} <: BMO.AbstractSDF{T}
    shape::S
    center::Point3{T}
    radius::R
    count::Int
    spheres::Int
end

WrappedSDF(shape::BMO.AbstractSDF{T}, radius = nothing; center = Point3{T}(0)) where {T} =
    WrappedSDF(shape, Point3{T}(center), radius, 0, 0)

BMO.position(w::WrappedSDF) = position(w.shape)
BMO.position!(w::WrappedSDF, pos) = BMO.position!(w.shape, pos)
BMO.orientation(w::WrappedSDF) = BMO.orientation(w.shape)
BMO.orientation!(w::WrappedSDF, dir) = BMO.orientation!(w.shape, dir)
BMO.sdf(w::WrappedSDF, point) = (w.count += 1; BMO.sdf(w.shape, point))

function BMO.bounding_sphere_of(w::WrappedSDF{<:Any, <:Any, <:Real})
    w.spheres += 1
    return BMO.SingleBoundingSphere(w, w.center, w.radius)
end

# A shape whose bounding sphere must never be asked for
struct ThrowingShape <: BMO.AbstractShape{Float64} end
BMO.bounding_sphere_of(::ThrowingShape) = error("the bounding sphere of this shape must not be computed")

# The test with the bounding sphere of the shape, outside of a solve
bounded(shape, ray) = BMO.intersect3d(BMO.bounding_sphere_of(shape), shape, ray)

function unit_vector(rng)
    v = randn(rng, 3)
    return v / norm(v)
end

same(::Nothing, ::Nothing) = true
same(a, b) = !isnothing(a) && !isnothing(b) && length(a) == length(b) && BMO.object(a) === BMO.object(b)

# Whether the sphere `outer` encloses the sphere `inner`
encloses(outer, inner) = norm(inner.pos - outer.pos) + inner.radius ≤ outer.radius + 1e-12

sphere_prism(radius, bound = radius) = BMO.Prism(WrappedSDF(BMO.SphereSDF(radius), bound), λ -> 1.5)

@testset "Bounding sphere" begin
    radius = 10e-3
    n = λ -> 1.5

    @testset "Shapes" begin
        plain = WrappedSDF(BMO.SphereSDF(radius))
        @test BMO.bounding_sphere_of(plain) === BMO.NoBoundingSphere()
        # the center is given in the local frame and converted for the current pose
        thickness, diameter = 4e-3, 20e-3
        plate = BMO.PlanoSurfaceSDF(thickness, diameter)
        r = sqrt((thickness / 2)^2 + (diameter / 2)^2)
        shape = WrappedSDF(plate, r; center = Point3(0, thickness / 2, 0))
        translate3d!(plate, [0.1, -0.2, 0.3])
        rotate3d!(plate, normalize([1.0, 2, 3]), 0.7)
        sphere = @inferred BMO.bounding_sphere_of(shape)
        @test sphere isa BMO.SingleBoundingSphere{Float64}
        @test sphere.radius == r
        @test sphere.pos ≈ position(plate) + BMO.orientation(plate) * [0, thickness / 2, 0]
        @test position(sphere) == sphere.pos
        # the center of the plate lies half its thickness below the surface
        @test BMO.sdf(plate, sphere.pos) ≈ -thickness / 2
        # mixed number types are promoted
        @test BMO.SingleBoundingSphere(Point3(0.0f0, 0, 0), 1.0) isa BMO.SingleBoundingSphere{Float64}
    end

    @testset "MultiBoundingSphere" begin
        none = BMO.NoBoundingSphere()
        a = BMO.SingleBoundingSphere([0.0, 0, 0], 1.0)
        b = BMO.SingleBoundingSphere([4.0, 0, 0], 2.0)
        # the smallest sphere around two spheres, which touch it from within
        ab = @inferred BMO.MultiBoundingSphere(a, b)
        @test ab isa BMO.MultiBoundingSphere{Float64}
        @test ab.pos ≈ [2.5, 0, 0]
        @test ab.radius ≈ 3.5
        @test position(ab) == ab.pos
        @test BMO.radius(ab) == ab.radius
        @test BMO.radius(a) == a.radius
        @test BMO.MultiBoundingSphere(b, a).pos ≈ ab.pos
        @test BMO.MultiBoundingSphere(b, a).radius ≈ ab.radius
        @test encloses(ab, a) && encloses(ab, b)
        # one sphere within the other, also the same sphere twice
        big = BMO.SingleBoundingSphere([0.5, 0, 0], 3.0)
        for (x, y) in ((a, big), (big, a), (big, big))
            inner = BMO.MultiBoundingSphere(x, y)
            @test inner isa BMO.MultiBoundingSphere{Float64}
            @test inner.pos == big.pos
            @test inner.radius == big.radius
        end
        # a sphere around several spheres is a part like any other
        c = BMO.SingleBoundingSphere([0.0, 10, 0], 1.0)
        abc = BMO.MultiBoundingSphere(ab, c)
        @test abc isa BMO.MultiBoundingSphere{Float64}
        @test all(s -> encloses(abc, s), (a, b, c, ab))
        # all spheres of a tuple
        all3 = @inferred BMO.MultiBoundingSphere((a, b, c))
        @test all3.pos ≈ abc.pos
        @test all3.radius ≈ abc.radius
        @test BMO.MultiBoundingSphere((a,)) == BMO.MultiBoundingSphere(a.pos, a.radius)
        @test BMO.MultiBoundingSphere((ab,)) === ab
        # no sphere if one of them is none, and for nothing at all
        @test BMO.MultiBoundingSphere(a, none) === none
        @test BMO.MultiBoundingSphere(none, a) === none
        @test BMO.MultiBoundingSphere(none, none) === none
        @test BMO.MultiBoundingSphere((a, none, b)) === none
        @test BMO.MultiBoundingSphere((none,)) === none
        @test BMO.MultiBoundingSphere(()) === none
        # the same sphere as the other kind
        @test BMO.SingleBoundingSphere(ab) == BMO.SingleBoundingSphere(ab.pos, ab.radius)
        @test BMO.MultiBoundingSphere(a) == BMO.MultiBoundingSphere(a.pos, a.radius)
        @test BMO.SingleBoundingSphere(a) === a
        @test BMO.SingleBoundingSphere(none) === none
        @test BMO.MultiBoundingSphere(none) === none
        # no sphere has no center and no radius, with an error that says so
        @test_throws ArgumentError position(none)
        @test_throws ArgumentError BMO.radius(none)
        @test_throws "NoBoundingSphere has no center" position(none)
        @test_throws "NoBoundingSphere has no radius" BMO.radius(none)
        @test_throws ArgumentError BMO.bounding_box(none)
        @test_throws ArgumentError BMO._sphere_exit(none, Ray([0.0, 0, 0], [0.0, 1, 0]))
        # number types are promoted
        @test BMO.MultiBoundingSphere(Point3(0.0f0, 0, 0), 1.0) isa BMO.MultiBoundingSphere{Float64}
        f32 = BMO.SingleBoundingSphere(Point3(1.0f0, 0, 0), 1.0f0)
        @test BMO.MultiBoundingSphere(f32, f32) isa BMO.MultiBoundingSphere{Float32}
        @test BMO.MultiBoundingSphere(a, f32) isa BMO.MultiBoundingSphere{Float64}
        # tested against a ray like the sphere of a single shape
        through = Ray([2.5, -10, 0], [0.0, 1, 0])
        passing = Ray([7.0, -10, 0], [0.0, 1, 0])
        @test BMO._sphere_exit(ab, through) ≈ BMO._sphere_exit(BMO.SingleBoundingSphere(ab), through)
        @test isnothing(BMO._sphere_exit(ab, passing))
        @test BMO.bounding_box(ab) == BMO.bounding_box(BMO.SingleBoundingSphere(ab))
    end

    @testset "Meshes" begin
        mesh = BMO.CubeMesh(20e-3)
        translate3d!(mesh, [0.1, 0.2, -0.3])
        rotate3d!(mesh, normalize([1.0, -1, 2]), 0.4)
        sphere = BMO.bounding_sphere_of(mesh)
        v = BMO.vertices(mesh)
        distances = [norm(v[i, :] - sphere.pos) for i in axes(v, 1)]
        @test maximum(distances) ≈ sphere.radius
        @test sphere.radius ≈ 10e-3 * sqrt(3)
        # an object that is never hit has none, and its shape is not asked
        dummy = BMO.NonInteractableObject(ThrowingShape())
        @test BMO.bounding_sphere_of(dummy) === BMO.NoBoundingSphere()
        @test BMO.bounding_sphere_of(BMO.NonInteractableObject(mesh)) === BMO.NoBoundingSphere()
    end

    @testset "Objects" begin
        prism = sphere_prism(radius)
        @test BMO.bounding_sphere_of(prism) == BMO.bounding_sphere_of(BMO.shape(prism))
        @test BMO.bounding_sphere_of(BMO.Prism(WrappedSDF(BMO.SphereSDF(radius)), n)) === BMO.NoBoundingSphere()
        # objects with several shapes and groups: the sphere around the spheres of all parts
        doublet = SphericalDoubletLens(100e-3, -60e-3, -200e-3, 6e-3, 3e-3, 25.4e-3, n, n)
        cube = CubeBeamsplitter(20e-3, n)
        translate3d!(cube, [0, 0.1, 0])
        zrotate3d!(cube, 0.3)
        inner = ObjectGroup([doublet, cube])
        other = sphere_prism(radius)
        translate3d!(other, [0.2, -0.1, 0.05])
        outer = ObjectGroup([inner, other])
        for object in (doublet, cube, inner, outer)
            main = BMO.bounding_sphere_of(object)
            @test main isa BMO.MultiBoundingSphere{Float64}
            @test all(part -> encloses(main, BMO.bounding_sphere_of(part)), BMO.shape(object))
        end
        # tight for two spheres: they touch the sphere around them from within
        a, b = sphere_prism(radius), sphere_prism(2radius)
        translate3d!(b, [0.1, 0, 0])
        main = BMO.bounding_sphere_of(ObjectGroup([a, b]))
        @test main.radius ≈ (0.1 + 3radius) / 2
        @test main.pos ≈ [(0.1 + radius) / 2, 0, 0]
        # one sphere within the other
        @test BMO.bounding_sphere_of(ObjectGroup([sphere_prism(radius), sphere_prism(3radius)])).radius ≈ 3radius
        # no sphere if a part has none
        none = BMO.Prism(WrappedSDF(BMO.SphereSDF(radius)), n)
        @test BMO.bounding_sphere_of(ObjectGroup([sphere_prism(radius), none])) === BMO.NoBoundingSphere()
    end

    @testset "Table of a solve" begin
        prism = sphere_prism(radius)
        doublet = SphericalDoubletLens(100e-3, -60e-3, -200e-3, 6e-3, 3e-3, 25.4e-3, n, n)
        group = ObjectGroup([sphere_prism(radius), sphere_prism(radius)])
        dummy = BMO.NonInteractableObject(ThrowingShape())
        none = BMO.Prism(WrappedSDF(BMO.SphereSDF(radius)), n)
        system = System([prism, doublet, group, dummy, none])
        table = BMO.bounding_spheres(system)
        @test table isa BMO.BoundingSphereTable
        # the shapes of the single-shape objects, and the multi-shape object and the group themselves
        shapes = [BMO.shape(prism), BMO.shape(doublet.front), BMO.shape(doublet.back),
            BMO.shape.(BMO.shape(group))...]
        @test all(s -> haskey(table, s), shapes)
        @test haskey(table, doublet) && haskey(table, group)
        @test length(table) == length(shapes) + 2
        @test !haskey(table, BMO.shape(dummy)) && !haskey(table, BMO.shape(none))
        # stored with a margin for the rounding of the shape
        computed = BMO.bounding_sphere_of(prism)
        @test table[BMO.shape(prism)].pos == computed.pos
        @test computed.radius < table[BMO.shape(prism)].radius < computed.radius * (1 + 1e-6)
        @test encloses(table[group], table[BMO.shape(first(BMO.shape(group)))])
        # the kind of the sphere is kept
        @test all(s -> table[s] isa BMO.SingleBoundingSphere{Float64}, shapes)
        @test table[doublet] isa BMO.MultiBoundingSphere{Float64}
        @test table[group] isa BMO.MultiBoundingSphere{Float64}
        @test isempty(BMO.bounding_spheres(System()))
        @test length(BMO.bounding_spheres(StaticSystem([prism, group]))) == 3

        # the lookup in a table, or in none
        shape = BMO.shape(prism)
        @test BMO.bounding_sphere_of(table, shape) === table[shape]
        @test BMO.bounding_sphere_of(table, group) === table[group]
        @test BMO.bounding_sphere_of(table, BMO.shape(none)) === BMO.NoBoundingSphere()
        @test BMO.bounding_sphere_of(nothing, shape) === BMO.NoBoundingSphere()
        @test BMO.bounding_sphere_of(nothing, group) === BMO.NoBoundingSphere()

        # the table of the current task is only set in a running solve
        @test isnothing(BMO.current_bounding_spheres())
        @test BMO.bounding_sphere_of(BMO.current_bounding_spheres(), shape) === BMO.NoBoundingSphere()
        stranger = sphere_prism(radius)
        BMO.with_bounding_spheres(system) do
            @test BMO.bounding_sphere_of(BMO.current_bounding_spheres(), shape) === table[shape]
            @test BMO.bounding_sphere_of(BMO.current_bounding_spheres(), group) === table[group]
            @test BMO.bounding_sphere_of(BMO.current_bounding_spheres(), BMO.shape(none)) === BMO.NoBoundingSphere()
            @test BMO.bounding_sphere_of(BMO.current_bounding_spheres(), BMO.shape(stranger)) === BMO.NoBoundingSphere()
            # a table that is set is kept
            BMO.with_bounding_spheres(System([stranger])) do
                @test BMO.bounding_sphere_of(BMO.current_bounding_spheres(), shape) === table[shape]
                @test BMO.bounding_sphere_of(BMO.current_bounding_spheres(), BMO.shape(stranger)) === BMO.NoBoundingSphere()
            end
        end
        @test BMO.bounding_sphere_of(BMO.current_bounding_spheres(), shape) === BMO.NoBoundingSphere()

        # the table belongs to the task: a spawned task sees none, unless it is passed on
        BMO.with_bounding_spheres(system) do
            current = BMO.current_bounding_spheres()
            @test current isa BMO.BoundingSphereTable
            @test isnothing(fetch(Threads.@spawn BMO.current_bounding_spheres()))
            passed = fetch(Threads.@spawn BMO.with_bounding_spheres(BMO.current_bounding_spheres, current))
            @test passed === current
            # a given table replaces the one that is set, and the latter is restored
            other = BMO.bounding_spheres(System([stranger]))
            @test BMO.with_bounding_spheres(BMO.current_bounding_spheres, other) === other
            @test BMO.current_bounding_spheres() === current
        end
        # no table is left if the function throws
        @test_throws ErrorException BMO.with_bounding_spheres(() -> error("failed"), system)
        @test isnothing(BMO.current_bounding_spheres())
    end

    @testset "Sphere before a shape or parts" begin
        a = BMO.SphereSDF(radius)
        b = BMO.SphereSDF(radius)
        translate3d!(b, [0, 4radius, 0])
        around = BMO.SingleBoundingSphere([0.0, 2radius, 0], 3radius)
        beside = BMO.SingleBoundingSphere([1.0, 0, 0], radius)
        hit = Ray([0.0, -3radius, 0], [0.0, 1, 0])
        passing = Ray([4radius, -3radius, 0], [0.0, 1, 0])
        # a shape, and the parts of an object or group as a tuple or a vector
        for x in (a, (b, a), [b, a])
            ref = BMO.intersect3d(x, hit)
            @test length(ref) ≈ 2radius
            @test length(BMO.intersect3d(BMO.NoBoundingSphere(), x, hit)) == length(ref)
            @test length(BMO.intersect3d(around, x, hit)) ≈ length(ref)
            @test isnothing(BMO.intersect3d(BMO.NoBoundingSphere(), x, passing))
            @test isnothing(BMO.intersect3d(around, x, passing))
            # the sphere decides: `x` is not tested if the ray misses the sphere
            @test isnothing(BMO.intersect3d(beside, x, hit))
        end
    end

    @testset "Evaluations of the shape" begin
        shape = WrappedSDF(BMO.SphereSDF(radius), 2radius)
        prism = BMO.Prism(shape, n)
        system = System([prism])
        passing = Ray([1.0, -1, 0], [0.0, 1, 0])
        behind = Ray([3radius, 0, 0], [1.0, 0, 0])
        hit = Ray([1.0, 0, 0], [-1.0, 0, 0])
        # outside of a solve the object is tested without its sphere
        @test isnothing(BMO.intersect3d(prism, passing))
        @test shape.count > 10
        @test shape.spheres == 0
        BMO.with_bounding_spheres(system) do
            # the line of the ray passes the sphere, and the sphere lies behind the ray
            shape.count = 0
            @test isnothing(BMO.intersect3d(prism, passing))
            @test isnothing(BMO.intersect3d(prism, behind))
            @test shape.count == 0
            # the ray passes through the sphere but misses the shape: marched up to the exit only
            through = Ray([1.5radius, -1, 0], [0.0, 1, 0])
            @test isnothing(BMO.intersect3d(prism, through))
            limited = shape.count
            shape.count = 0
            @test isnothing(BMO.intersect3d(shape, through))
            @test 0 < limited < shape.count
            # a ray that leaves the shape from its surface
            for dir in ([1.0, 0, 0], normalize([1.0, 1, 0]), normalize([1.0, 1, 1]))
                leaving = Ray([radius, 0, 0], dir)
                shape.count = 0
                @test isnothing(BMO.intersect3d(prism, leaving))
                limited = shape.count
                shape.count = 0
                @test isnothing(BMO.intersect3d(shape, leaving))
                @test limited < 60
                @test limited < shape.count / 2
            end
            # a hit is unchanged
            found = BMO.intersect3d(prism, hit)
            @test BMO.object(found) === prism
            @test length(found) == length(BMO.intersect3d(shape, hit)) ≈ 1 - radius
        end
    end

    @testset "Groups" begin
        members = [sphere_prism(radius) for _ in 1:5]
        foreach((m, i) -> translate3d!(m, [0.1, 0.03i, 0]), members, 1:5)
        target = sphere_prism(radius)
        translate3d!(target, [0, 0.2, 0])
        group = ObjectGroup(members)
        system = System([group, target])
        count() = sum(m -> BMO.shape(m).count, members)
        ray = Ray([0.0, -1, 0], [0.0, 1, 0])
        # the ray misses the sphere of the group: no member is evaluated
        BMO.with_bounding_spheres(system) do
            @test BMO.object(BMO.trace_all(system, ray)) === target
            @test count() == 0
        end
        beam = Beam([0.0, -1, 0], [0.0, 1, 0], 1e-6)
        solve_system!(system, beam)
        @test count() == 0
        @test BMO.object(BMO.intersection(first(BMO.rays(beam)))) === target
        # the hit of a member keeps the member, in a solve and for a direct call
        towards = Ray([0.1, -1, 0], [0.0, 1, 0])
        @test BMO.object(BMO.intersect3d(group, towards)) === members[1]
        BMO.with_bounding_spheres(system) do
            @test BMO.object(BMO.trace_all(system, towards)) === members[1]
        end
        # also for a group within a group
        nested = ObjectGroup([ObjectGroup([sphere_prism(radius)]), sphere_prism(radius)])
        @test BMO.object(BMO.intersect3d(nested, ray)) === only(BMO.shape(first(BMO.shape(nested))))
    end

    @testset "Same result as without spheres" begin
        rng = MersenneTwister(3)
        lens = SphericalLens(50e-3, -50e-3, 8e-3, 25.4e-3, n)
        doublet = SphericalDoubletLens(100e-3, -60e-3, -200e-3, 6e-3, 3e-3, 25.4e-3, n, n)
        translate3d!(doublet, [0, 30e-3, 0])
        plate = RoundPlateBeamsplitter(25.4e-3, 5e-3, n)
        translate3d!(plate, [0, 60e-3, 0])
        zrotate3d!(plate, deg2rad(45))
        cube = CubeBeamsplitter(20e-3, n)
        translate3d!(cube, [0, 90e-3, 0])
        mirror = SquarePlanoMirror2D(25.4e-3)
        translate3d!(mirror, [0, 120e-3, 0])
        dummy = BMO.NonInteractableObject(BMO.CubeMesh(15e-3))
        translate3d!(dummy, [30e-3, 60e-3, 0])
        a, b = SphericalLens(80e-3, Inf, 5e-3, 25.4e-3, n), RoundPlanoMirror(25.4e-3, 5e-3)
        translate3d!(a, [-40e-3, 20e-3, 0])
        translate3d!(b, [-40e-3, 80e-3, 10e-3])
        group = ObjectGroup([ObjectGroup([a]), b])
        system = System([lens, doublet, plate, cube, mirror, dummy, group])
        center = [0, 60e-3, 0]
        rays = map(1:10^4) do _
            start = center + 0.25 * unit_vector(rng)
            target = center + [50e-3, 80e-3, 30e-3] .* (2 .* rand(rng, 3) .- 1)
            Ray(start, target - start)
        end
        exact = [BMO.trace_all(system, ray) for ray in rays]
        with_table = BMO.with_bounding_spheres(system) do
            [BMO.trace_all(system, ray) for ray in rays]
        end
        @test all(same.(exact, with_table))
        hits = filter(!isnothing, exact)
        @test 1000 < length(hits) < 9000
        # every kind of object is among the hits, the dummy never
        hit_objects = unique(BMO.object.(hits))
        @test all(o -> any(h -> h === o, hit_objects), (lens, doublet, plate, cube, mirror, a, b))
        @test !any(h -> h === dummy, hit_objects)
    end

    @testset "Solves" begin
        shape = WrappedSDF(BMO.SphereSDF(radius), radius)
        prism = BMO.Prism(shape, n)
        translate3d!(prism, [0, 0.1, 0])
        dummy = BMO.NonInteractableObject(ThrowingShape())
        system = System([prism, dummy])
        # one table per solve of a beam, and one for all beams of a group
        beam = Beam([0.0, 0, 0], [0.0, 1, 0], 1e-6)
        solve_system!(system, beam)
        @test shape.spheres == 1
        @test length(BMO.rays(beam)) > 1
        source = CollimatedSource([0.0, 0, 0], [0.0, 1, 0], 5e-3, 1e-6; num_rings = 5)
        shape.spheres = 0
        solve_system!(system, source)
        @test shape.spheres == 1
        @test all(b -> length(BMO.rays(b)) > 1, BMO.beams(source))
        # the next solve sees a moved object
        translate3d!(prism, [1.0, 0, 0])
        solve_system!(system, beam)
        @test length(BMO.rays(beam)) == 1
        translate3d!(prism, [-1.0, 0, 0])
        solve_system!(system, beam)
        @test BMO.object(BMO.intersection(first(BMO.rays(beam)))) === prism
        # a hint is tested with the sphere of its shape
        BMO.with_bounding_spheres(system) do
            shape.count = 0
            passing = Ray([1.0, -1, 0], [0.0, 1, 0])
            @test isnothing(BMO.trace_one(system, passing, BMO.Hint(prism)))
            @test shape.count == 0
            hinted = BMO.trace_one(system, Ray([0.0, -1, 0], [0.0, 1, 0]), BMO.Hint(prism))
            @test BMO.object(hinted) === prism
            @test length(hinted) ≈ 1.1 - radius
        end
    end

    @testset "Shape that touches its sphere" begin
        rng = MersenneTwister(1)
        shape = WrappedSDF(BMO.SphereSDF(radius), radius)
        translate3d!(shape.shape, [0.3, -0.1, 0.2])
        c = position(shape)
        hits = 0
        for _ in 1:1000
            start = c + 5radius * unit_vector(rng)
            target = c + 2radius * rand(rng) * unit_vector(rng)
            ray = Ray(start, target - start)
            exact = BMO.intersect3d(shape, ray)
            @test same(exact, bounded(shape, ray))
            hits += !isnothing(exact)
        end
        @test 100 < hits < 900
        # tangent rays
        for _ in 1:1000
            p = unit_vector(rng)
            u = normalize(cross(p, unit_vector(rng)))
            ray = Ray(c + radius * p - 3radius * u, u)
            @test same(BMO.intersect3d(shape, ray), bounded(shape, ray))
        end
    end

    @testset "Float32 shape, Float64 ray" begin
        rng = MersenneTwister(2)
        r = 10.0f-3
        shape = WrappedSDF(BMO.SphereSDF(r), r)
        translate3d!(shape.shape, Float32[1, 2, 3])
        prism = BMO.Prism(shape, λ -> 1.5f0)
        system = System([prism])
        c = Float64.(position(shape))
        @test BMO.bounding_sphere_of(shape) isa BMO.SingleBoundingSphere{Float32}
        stored = BMO.bounding_spheres(system)[shape]
        @test stored isa BMO.SingleBoundingSphere{Float64}
        # the margin covers the rounding of a Float32 pose
        @test stored.radius - r > 1e-4 * norm(c)
        rays = map(1:1000) do _
            start = c + 5r * unit_vector(rng)
            target = c + 2r * rand(rng) * unit_vector(rng)
            Ray(start, target - start)
        end
        exact = [BMO.intersect3d(prism, ray) for ray in rays]
        BMO.with_bounding_spheres(system) do
            @test all(same.(exact, [BMO.intersect3d(prism, ray) for ray in rays]))
            shape.count = 0
            @test isnothing(BMO.intersect3d(prism, Ray(c + [1.0, -1, 0], [0.0, 1, 0])))
            @test shape.count == 0
        end
        @test 100 < count(!isnothing, exact) < 900
    end
end

end # MODULE
