module TestAbstractSDF

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics

const BMO = BeamletOptics

# Orientation-less test point sdf
mutable struct TestPointSDF{T} <: BMO.AbstractSDF{T}
    position::Point3{T}
    orientation::Matrix{T}
end

TestPointSDF(p::AbstractArray{T}) where {T} = TestPointSDF{T}(
    Point3{T}(p), Matrix{T}(I, 3, 3))
TestPointSDF(T = Float64) = TestPointSDF{T}(Point3{T}(0), Matrix{T}(I, 3, 3))

BMO.position(tps::TestPointSDF) = tps.position
BMO.position!(tps::TestPointSDF{T}, new::Point3{T}) where {T} = (tps.position = new)

BMO.orientation(tps::TestPointSDF) = tps.orientation
BMO.orientation!(tps::TestPointSDF{T}, new::Matrix{T}) where {T} = (tps.orientation = new)

BMO.transposed_orientation(tps::TestPointSDF) = transpose(tps.orientation)
BMO.transposed_orientation!(::TestPointSDF, ::Any) = nothing

function BMO.sdf(tps::TestPointSDF, point)
    p = BMO._world_to_sdf(tps, point)
    return norm(p)
end

# Counts the evaluations of the sdf of the wrapped shape
mutable struct CountingSDF{T, S <: BMO.AbstractSDF{T}} <: BMO.AbstractSDF{T}
    shape::S
    count::Int
end

CountingSDF(shape) = CountingSDF(shape, 0)

BMO.sdf(c::CountingSDF, point) = (c.count += 1; BMO.sdf(c.shape, point))

@testset "Abstract SDF" begin
    @testset "Testing type definitions" begin
        @test isdefined(BMO, :AbstractSDF)
        @test isdefined(BMO, :SphereSDF)
        @test isdefined(BMO, :CylinderSDF)
        @test isdefined(BMO, :CutSphereSDF)
        @test isdefined(BMO, :ThinLensSDF)
    end

    @testset "Testing kinematics and transforms" begin
        point = TestPointSDF()
        t = 10
        θ = deg2rad(30)
        translate3d!(point, [t, 0, 0])
        rotate3d!(point, [0, 1, 0], θ)
        pt = BMO._world_to_sdf(point, [0, 0, 0])
        @test pt[1] ≈ -t * cos(θ)
        @test pt[2] ≈ 0
        @test pt[3] ≈ -t * sin(θ)
    end

    @testset "Testing intersect3d" begin
        t = 10.0
        point = TestPointSDF(zeros(3))
        translate3d!(point, [t, 0, 0])

        r1 = Ray(zeros(3), [1.0, 0, 0])
        r2 = Ray(zeros(3), [1.0, 1, 0])
        r3 = Ray(zeros(3), [1.0, 0, 1])

        i1 = BMO.intersect3d(point, r1)
        i2 = BMO.intersect3d(point, r2)
        i3 = BMO.intersect3d(point, r3)

        @test length(i1) == t
        @test isnothing(i2)
        @test isnothing(i3)
    end

    @testset "Testing misses" begin
        radius = 10e-3
        counting = CountingSDF(BMO.SphereSDF(radius))
        # a ray that leaves the shape is a miss after few evaluations, whatever its start distance
        for distance in (1e-3, 1.0, 1e3), dir in ([0.0, 1, 0], [1.0, 0, 0], [1.0, 1, 1] / sqrt(3))
            counting.count = 0
            ray = Ray([radius + distance, 0, 0], dir)
            @test isnothing(BMO.intersect3d(counting, ray))
            @test counting.count < 100
        end
        # a hit needs few evaluations as well
        counting.count = 0
        hit = BMO.intersect3d(counting, Ray([1.0, 0, 0], [-1.0, 0, 0]))
        @test length(hit) ≈ 1 - radius
        @test counting.count < 100
        # a shape far away is still hit, one beyond the miss distance is not
        far = BMO.SphereSDF(1.0)
        translate3d!(far, [0, 1e6, 0])
        @test length(BMO.intersect3d(far, Ray(zeros(3), [0.0, 1, 0]))) ≈ 1e6 - 1
        translate3d!(far, [0, 10 * get_sdf_miss_distance(), 0])
        @test isnothing(BMO.intersect3d(far, Ray(zeros(3), [0.0, 1, 0])))
        # the miss distance of a single test
        near = BMO.SphereSDF(1.0)
        translate3d!(near, [0, 10, 0])
        @test length(BMO.intersect3d(near, Ray(zeros(3), [0.0, 1, 0]), 9.5)) ≈ 9
        @test isnothing(BMO.intersect3d(near, Ray(zeros(3), [0.0, 1, 0]), 8.5))
    end

    @testset "thickness on the optical axis" begin
        # the `height` field of a cylinder is half its height
        @test BMO.thickness(BMO.CylinderSDF(0.5, 1.0)) === 2.0
        # no material on the axis of a ring
        @test BMO.thickness(BMO.RingSDF(10.0f-3, 2.0f-3, 5.0f-3)) === 0.0f0
        # the fallback throws instead of returning 0
        @test_throws "thickness not defined for" BMO.thickness(BMO.CutSphereSDF(1.0, 0.5))
    end

    @testset "Testing normal3d" begin
        point = TestPointSDF(zeros(3))
        offset = [5, 0, 0]
        translate3d!(point, offset)
        p1 = [1, 0, 0]
        p2 = [0, 1, 0]
        p3 = [0, 0, 1]
        @test BMO.normal3d(point, p1 + offset) == p1
        @test BMO.normal3d(point, p2 + offset) == p2
        @test BMO.normal3d(point, p3 + offset) == p3
    end

    @testset "Rays inside leave at the first exit" begin
        mm = 1e-3
        # a ray within the wall of a ring (radius 10 to 12 mm) that runs across the hole leaves the
        # wall at the inner face, not at the far side of the ring, also if it starts in the middle
        # of the wall, where the SDF increases along the ray
        ring = BMO.RingSDF(10mm, 2mm, 5mm)
        i = BMO.intersect3d(ring, Ray([-11mm, 0, 0], [1.0, 0, 0]))
        @test i.t ≈ 1mm atol = 1e-9
        @test BMO.normal3d(i) ≈ [1, 0, 0] atol = 1e-6
        # the same from the surface: the ray enters the wall at -12 mm, the exit is at -10 mm
        i = BMO.intersect3d(ring, Ray([-12mm, 0, 0], [1.0, 0, 0]))
        @test i.t ≈ 2mm atol = 1e-9
        # from the inner face into the hole, the next hit is the opposite wall
        i = BMO.intersect3d(ring, Ray([-10mm, 0, 0], [1.0, 0, 0]))
        @test i.t ≈ 20mm atol = 1e-9
        # a lens is a union of parts: the faces between the parts within the lens are no exit
        lens = BMO.shape(SphericalLens(50mm, -50mm, 5mm, 25.4mm, λ -> 1.5))
        front = BMO.intersect3d(lens, Ray([2mm, -50mm, 1mm], [0.0, 1, 0]))
        p = Point3(2mm, -50mm, 1mm) + front.t * Point3(0.0, 1, 0)
        d = normalize(Point3(0.1, 1, 0.05))
        i = BMO.intersect3d(lens, Ray(p, d))
        @test BMO.sdf(lens, p + i.t * d) ≈ 0 atol = 1e-9
        @test BMO.sdf(lens, p + (i.t - 1e-6) * d) < 0
        @test BMO.sdf(lens, p + (i.t + 1e-6) * d) > 0
    end

    @testset "Hits from outside lie on the surface" begin
        # the hit is the root of the SDF, not the first step within the marching tolerance, hence it
        # varies smoothly with the ray
        mm = 1e-3
        R = 10mm
        sphere = BMO.SphereSDF(R)
        for b in (0.0, 0.5R, 0.9R, 0.99R)
            ray = Ray([b, -100mm, 0.0], [0.0, 1.0, 0])
            i = BMO.intersect3d(sphere, ray)
            @test i.t ≈ 100mm - sqrt(R^2 - b^2) atol = 1e-13
            @test BMO.normal3d(i) ≈ normalize([b, -sqrt(R^2 - b^2), 0]) atol = 1e-9
        end
        # the same for a cylinder cap hit under an angle
        cyl = BMO.CylinderSDF(R, 5mm)
        d = normalize([0.1, 1.0, 0.05])
        p0 = Point3(0.0, -50mm, 0.0)
        i = BMO.intersect3d(cyl, Ray(p0, d))
        @test p0[2] + i.t * d[2] ≈ -5mm atol = 1e-13
        @test BMO.normal3d(i) ≈ [0, -1, 0]
        # finite differences of the hit distance are smooth across the aperture: hits at the first
        # probe within the marching tolerance put relative errors of up to 1e-4 into them
        t(x) = BMO.intersect3d(sphere, Ray([x, -100mm, 0.0], [0.0, 1.0, 0])).t
        h = 1e-7
        rel_err(x) = abs((t(x + h) - t(x - h)) / 2h / (x / sqrt(R^2 - x^2)) - 1)
        @test maximum(rel_err, range(1mm, 9mm, 81)) < 1e-6
    end
end

end # MODULE
