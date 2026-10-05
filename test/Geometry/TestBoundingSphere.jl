module TestBoundingSphere

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics
using Random

const BMO = BeamletOptics

# Wraps an sdf: counts its evaluations and, with a `radius`, has a bounding sphere around `center`
mutable struct WrappedSDF{T, S <: BMO.AbstractSDF{T}, R} <: BMO.AbstractSDF{T}
    shape::S
    center::Point3{T}
    radius::R
    count::Int
end

WrappedSDF(shape::BMO.AbstractSDF{T}, radius = nothing; center = Point3{T}(0)) where {T} =
    WrappedSDF(shape, Point3{T}(center), radius, 0)

BMO.position(w::WrappedSDF) = position(w.shape)
BMO.position!(w::WrappedSDF, pos) = BMO.position!(w.shape, pos)
BMO.orientation(w::WrappedSDF) = BMO.orientation(w.shape)
BMO.orientation!(w::WrappedSDF, dir) = BMO.orientation!(w.shape, dir)
BMO.sdf(w::WrappedSDF, point) = (w.count += 1; BMO.sdf(w.shape, point))

BMO.bounding_sphere(w::WrappedSDF{<:Any, <:Any, <:Real}) = (w.center, w.radius)

# The test with the bounding sphere, as called by the solver
bounded(shape, ray) = BMO.intersect3d(BMO.world_bounding_sphere(shape), shape, ray)

# A point on the unit sphere
function unit_vector(rng)
    v = randn(rng, 3)
    return v / norm(v)
end

same(::Nothing, ::Nothing) = true
same(a, b) = !isnothing(a) && !isnothing(b) && length(a) == length(b)

@testset "Bounding sphere" begin
    radius = 10e-3

    @testset "Shape without a sphere" begin
        mesh = BMO.QuadraticFlatMesh(0.05)
        @test isnothing(BMO.bounding_sphere(mesh))
        @test isnothing(BMO.world_bounding_sphere(mesh))
        hit = Ray([0.0, -1, 0], [0.0, 1, 0])
        miss = Ray([1.0, -1, 0], [0.0, 1, 0])
        @test length(bounded(mesh, hit)) == length(BMO.intersect3d(mesh, hit)) == 1
        @test isnothing(bounded(mesh, miss))
        # an sdf without a sphere is marched as before
        plain = WrappedSDF(BMO.SphereSDF(radius))
        @test isnothing(BMO.world_bounding_sphere(plain))
        @test isnothing(bounded(plain, miss))
        @test plain.count > 10
    end

    @testset "Sphere in world coordinates" begin
        thickness, diameter = 4e-3, 20e-3
        plate = BMO.PlanoSurfaceSDF(thickness, diameter)
        r = sqrt((thickness / 2)^2 + (diameter / 2)^2)
        shape = WrappedSDF(plate, r; center = Point3(0, thickness / 2, 0))
        translate3d!(plate, [0.1, -0.2, 0.3])
        rotate3d!(plate, normalize([1.0, 2, 3]), 0.7)
        center, r_world = @inferred BMO.world_bounding_sphere(shape)
        @test r_world == r
        @test center ≈ position(plate) + BMO.orientation(plate) * [0, thickness / 2, 0]
        # the center of the plate lies half its thickness below the surface
        @test BMO.sdf(plate, center) ≈ -thickness / 2
    end

    @testset "Evaluations of the shape" begin
        shape = WrappedSDF(BMO.SphereSDF(radius), 2radius)
        # the line of the ray passes the sphere
        @test isnothing(bounded(shape, Ray([1.0, -1, 0], [0.0, 1, 0])))
        @test shape.count == 0
        # the sphere lies behind the ray
        @test isnothing(bounded(shape, Ray([3radius, 0, 0], [1.0, 0, 0])))
        @test shape.count == 0
        # the ray passes through the sphere but misses the shape: marched up to the exit only
        through = Ray([1.5radius, -1, 0], [0.0, 1, 0])
        @test isnothing(bounded(shape, through))
        limited = shape.count
        shape.count = 0
        @test isnothing(BMO.intersect3d(shape, through))
        @test 0 < limited < shape.count
        # a ray that leaves the shape from its surface
        for dir in ([1.0, 0, 0], normalize([1.0, 1, 0]), normalize([1.0, 1, 1]))
            leaving = Ray([radius, 0, 0], dir)
            shape.count = 0
            @test isnothing(bounded(shape, leaving))
            limited = shape.count
            shape.count = 0
            @test isnothing(BMO.intersect3d(shape, leaving))
            @test limited < 60
            @test limited < shape.count / 2
        end
        # a hit is unchanged
        hit = Ray([1.0, 0, 0], [-1.0, 0, 0])
        @test length(bounded(shape, hit)) == length(BMO.intersect3d(shape, hit)) ≈ 1 - radius
    end

    @testset "Paths of the solver" begin
        shape = WrappedSDF(BMO.SphereSDF(radius), radius)
        prism = BMO.Prism(shape, λ -> 1.5)
        other = BMO.Prism(WrappedSDF(BMO.SphereSDF(radius), radius), λ -> 1.5)
        translate3d!(other, [0, 0.1, 0])
        count() = shape.count + BMO.shape(other).count
        passing = Ray([1.0, -1, 0], [0.0, 1, 0])
        hitting = Ray([0.0, -1, 0], [0.0, 1, 0])
        # object with a single shape
        @test isnothing(BMO.intersect3d(prism, passing))
        @test count() == 0
        @test BMO.object(BMO.intersect3d(prism, hitting)) === prism
        # object with several shapes
        group = ObjectGroup([prism, other])
        shape.count = 0
        @test isnothing(BMO.intersect3d(group, passing))
        @test count() == 0
        @test length(BMO.intersect3d(group, hitting)) ≈ 1 - radius
        # hint
        system = System([prism, other])
        shape.count = BMO.shape(other).count = 0
        @test isnothing(BMO.trace_one(system, passing, BMO.Hint(other)))
        @test count() == 0
        hinted = BMO.trace_one(system, hitting, BMO.Hint(other))
        @test BMO.object(hinted) === other
        @test length(hinted) ≈ 1.1 - radius
        # the full search
        @test BMO.object(BMO.trace_all(system, hitting)) === prism
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
        c = Float64.(position(shape))
        @test BMO.world_bounding_sphere(shape) isa Tuple{Point3{Float32}, Float32}
        @test isnothing(bounded(shape, Ray(c + [1.0, -1, 0], [0.0, 1, 0])))
        @test shape.count == 0
        for _ in 1:200
            start = c + 5r * unit_vector(rng)
            target = c + 2r * rand(rng) * unit_vector(rng)
            ray = Ray(start, target - start)
            @test same(BMO.intersect3d(shape, ray), bounded(shape, ray))
        end
    end
end

end # MODULE
