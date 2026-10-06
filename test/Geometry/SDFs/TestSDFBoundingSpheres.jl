module TestSDFBoundingSpheres

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics
using Random

const BMO = BeamletOptics

# SDF without a `bounding_sphere_of` method
mutable struct NoSphereSDF <: BMO.AbstractSDF{Float64}
    dir::Matrix{Float64}
    transposed_dir::Matrix{Float64}
    pos::Point3{Float64}
end
NoSphereSDF() = NoSphereSDF(Matrix(1.0I, 3, 3), Matrix(1.0I, 3, 3), Point3(0.0))
BMO.sdf(s::NoSphereSDF, p) = norm(BMO._world_to_sdf(s, p)) - 5e-3

const N_RAYS = 10^4
const ENCLOSURE_TOL = 1e-9  # [m]
const TIGHTNESS = 0.6

"""
Returns `label => shape` pairs that cover the 20 concrete SDF types, both signs of curvature where
a type has them, and the shapes of the public lens, prism and mirror constructors.
"""
function instances()
    n = λ -> 1.5
    asph = [0, 1e3]
    acyl = [0, 1.1926075e-5 * (1e3)^3, -2.9323497e-9 * (1e3)^5]
    # operands that are moved one by one, such that their frames differ from that of the composite
    box = BMO.BoxSDF(20e-3, 10e-3, 30e-3)
    translate3d!(box, [5e-3, 20e-3, -10e-3])
    xrotate3d!(box, deg2rad(25))
    bore = BMO.CylinderSDF(3e-3, 20e-3)
    translate3d!(bore, [5e-3, 20e-3, -10e-3])
    s1, s2, s3 = BMO.SphereSDF(4e-3), BMO.SphereSDF(7e-3), BMO.SphereSDF(2e-3)
    translate3d!(s2, [30e-3, 5e-3, 0])
    translate3d!(s3, [0, -20e-3, 10e-3])
    inner = BMO.SphereSDF(2e-3)
    translate3d!(inner, [1e-3, 0, 0])
    return [
        # PrimitiveSDF.jl
        "BoxSDF" => BMO.BoxSDF(20e-3, 10e-3, 30e-3),
        "CylinderSDF" => BMO.CylinderSDF(12e-3, 4e-3),
        "CutSphereSDF, cap" => BMO.CutSphereSDF(10e-3, 4e-3),
        "CutSphereSDF, more than a hemisphere" => BMO.CutSphereSDF(10e-3, -3e-3),
        "RingSDF" => BMO.RingSDF(10e-3, 4e-3, 6e-3),
        "RightAnglePrismSDF" => BMO.RightAnglePrismSDF(20e-3, 15e-3),
        "PolygonPrismSDF, off-center pentagon" => BMO.PolygonPrismSDF(
            [(0.0, 0.0), (20e-3, 0.0), (25e-3, 10e-3), (8e-3, 18e-3), (-4e-3, 9e-3)], 7e-3),
        "PolygonPrismSDF, right triangle" => BMO.PolygonPrismSDF(
            [(0.0, 0.0), (20e-3, 0.0), (0.0, 20e-3)], 30e-3),
        "EquilateralPrism" => BMO.shape(EquilateralPrism(20e-3, 10e-3, 1.5)),
        "DovePrism" => BMO.shape(DovePrism(50e-3, 10e-3, 10e-3, 1.5)),
        "RightAnglePrism" => BMO.shape(RightAnglePrism(25e-3, 20e-3, 1.5)),
        # SphericalLensSDF.jl
        "PlanoSurfaceSDF" => BMO.PlanoSurfaceSDF(5e-3, 25e-3),
        "SphereSDF" => BMO.SphereSDF(8e-3),
        "ConcaveSphericalSurfaceSDF" => BMO.ConcaveSphericalSurfaceSDF(30e-3, 25e-3),
        "ConcaveSphericalSurfaceSDF, deep" => BMO.ConcaveSphericalSurfaceSDF(13e-3, 25e-3),
        "ConvexSphericalSurfaceSDF" => BMO.ConvexSphericalSurfaceSDF(30e-3, 25e-3),
        "ConvexSphericalSurfaceSDF, deep" => BMO.ConvexSphericalSurfaceSDF(13e-3, 25e-3),
        # CylindricalSDF.jl
        "ConvexCylinderSDF" => BMO.ConvexCylinderSDF(30e-3, 20e-3, 25e-3),
        "ConcaveCylinderSDF, R > 0" => BMO.ConcaveCylinderSDF(30e-3, 20e-3, 25e-3),
        "ConcaveCylinderSDF, R < 0" => BMO.ConcaveCylinderSDF(-30e-3, 20e-3, 25e-3),
        # AsphericalLensSDF.jl
        "ConvexAsphericalSurfaceSDF, R > 0" => BMO.ConvexAsphericalSurfaceSDF(asph, 20e-3, -1.0, 10e-3),
        "ConvexAsphericalSurfaceSDF, R < 0" => BMO.ConvexAsphericalSurfaceSDF(-asph, -20e-3, -1.0, 10e-3),
        "ConcaveAsphericalSurfaceSDF, R < 0" => BMO.ConcaveAsphericalSurfaceSDF(-asph, -20e-3, -1.0, 10e-3),
        "ConcaveAsphericalSurfaceSDF, R > 0" => BMO.ConcaveAsphericalSurfaceSDF(asph, 20e-3, -1.0, 10e-3),
        # AcylindricalSDF.jl
        "AconvexCylinderSDF" => BMO.AconvexCylinderSDF(15.538e-3, 25e-3, 50e-3, -1.0, acyl),
        "AconcaveCylinderSDF, R < 0" => BMO.AconcaveCylinderSDF(-15.538e-3, 25e-3, 50e-3, -1.0, -acyl),
        "AconcaveCylinderSDF, R > 0" => BMO.AconcaveCylinderSDF(15.538e-3, 25e-3, 50e-3, -1.0, acyl),
        # ConicSDF.jl
        "ConicSDF, paraboloid" => BMO.ConicSDF(0.2, -1, 0, 60e-3, 20e-3),
        "ConicSDF, off-axis paraboloid" => BMO.ConicSDF(0.2, -1, 50e-3, 50e-3, 30e-3),
        "ConicSDF, convex sphere" => BMO.ConicSDF(-0.5, 0, 0, 0.2, 50e-3),
        "ConicSDF, off-axis hyperboloid" => BMO.ConicSDF(0.3, -2.5, 40e-3, 50e-3, 25e-3),
        "ConicSDF, off-axis convex ellipsoid" => BMO.ConicSDF(-0.3, -0.5, 40e-3, 50e-3, 25e-3),
        "OffAxisParabolicMirror" => BMO.shape(OffAxisParabolicMirror(50e-3, 25e-3; angle = 90)),
        # MeniscusLensSDF.jl
        "MeniscusLensSDF, R > 0" => BMO.MeniscusLensSDF(50e-3, 30e-3, 3e-3, 25.4e-3),
        "MeniscusLensSDF, R < 0" => BMO.MeniscusLensSDF(-30e-3, -50e-3, 3e-3, 25.4e-3),
        "MeniscusLensSDF, thick" => BMO.MeniscusLensSDF(20e-3, 25e-3, 15e-3, 25.4e-3),
        "MeniscusLensSDF, deep" => BMO.MeniscusLensSDF(13e-3, 14e-3, 3e-3, 25.4e-3),
        "MeniscusLensSDF with ring" => BMO.MeniscusLensSDF(50e-3, 30e-3, 3e-3, 25.4e-3, 30e-3),
        "meniscus Lens, R > 0" => BMO.shape(SphericalLens(50e-3, 30e-3, 3e-3, 25.4e-3, n)),
        "meniscus Lens, R < 0" => BMO.shape(SphericalLens(-30e-3, -50e-3, 3e-3, 25.4e-3, n)),
        # UnionSDF.jl
        "UnionSDF, separate spheres" => s1 + s2 + s3,
        "bi-convex SphericalLens" => BMO.shape(SphericalLens(50e-3, -50e-3, 5e-3, 25.4e-3, n)),
        "plano-convex SphericalLens" => BMO.shape(SphericalLens(50e-3, Inf, 5e-3, 25.4e-3, n)),
        "convex-plano SphericalLens" => BMO.shape(SphericalLens(Inf, -50e-3, 5e-3, 25.4e-3, n)),
        "bi-concave SphericalLens" => BMO.shape(SphericalLens(-50e-3, 50e-3, 3e-3, 25.4e-3, n)),
        "plano-concave SphericalLens" => BMO.shape(SphericalLens(Inf, 50e-3, 3e-3, 25.4e-3, n)),
        "ThinLens" => BMO.shape(ThinLens(50e-3, 80e-3, 25.4e-3, n)),
        "Lens with mechanical diameter" => BMO.shape(Lens(SphericalSurface(-40e-3, 20e-3, 25e-3),
            SphericalSurface(60e-3, 16e-3, 25e-3), 3e-3, n)),
        "aspheric Lens, convex" => BMO.shape(Lens(EvenAsphericalSurface(20e-3, 10e-3, -1.0, asph), 4e-3, n)),
        "aspheric Lens, concave" => BMO.shape(Lens(EvenAsphericalSurface(-20e-3, 10e-3, -1.0, -asph), 4e-3, n)),
        "cylindrical Lens, plano-convex" => BMO.shape(Lens(CylindricalSurface(5.2e-3, 10e-3, 20e-3), 5.9e-3, n)),
        "cylindrical Lens, bi-concave" => BMO.shape(Lens(CylindricalSurface(-20e-3, 10e-3, 20e-3),
            CylindricalSurface(20e-3, 10e-3, 20e-3), 3e-3, n)),
        "acylindrical Lens" => BMO.shape(Lens(BMO.AcylindricalSurface(15.538e-3, 25e-3, 50e-3, -1.0, acyl), 7.5e-3, n)),
        "SphericalMirror" => BMO.shape(SphericalMirror(100e-3, 25e-3, 5e-3)),
        # DifferenceSDF.jl
        "DifferenceSDF, moved base" => box - bore,
        "DifferenceSDF, inner tool" => BMO.SphereSDF(6e-3) - inner,
        "RoundPlanoMirror with hole" => BMO.shape(RoundPlanoMirror(25e-3, 5e-3; hole_diameter = 5e-3)),
        "ParabolicMirror with hole" => BMO.shape(ParabolicMirror(25e-3, 50e-3; hole_diameter = 5e-3)),
    ]
end

allocations(s) = @allocated BMO.bounding_sphere_of(s)

random_direction(rng) = normalize(Point3(randn(rng), randn(rng), randn(rng)))

# The test with the bounding sphere, as called by the solver
bounded(shape, ray) = BMO.intersect3d(BMO.bounding_sphere_of(shape), shape, ray)

same(::Nothing, ::Nothing) = true
same(a, b) = !isnothing(a) && !isnothing(b) && a.t == b.t

"""
Returns the number of hits, the largest distance of a hit point from `center` and the number of
rays for which the test with the bounding sphere differs from the exact one, for `N_RAYS` rays that
start outside of the sphere `(center, r)` and aim at random points within and around it, and for
one ray from each hit point into a random direction, i.e. leaving or entering the shape. The hits
stem from the exact test `intersect3d(shape, ray)`, which does not use the bounding sphere.
"""
function farthest_hit(rng, shape, center, r)
    hits, dmax, differing = 0, 0.0, 0
    for _ in 1:N_RAYS
        pos = center + 3r * random_direction(rng)
        target = center + 1.3r * cbrt(rand(rng)) * random_direction(rng)
        ray = Ray(pos, target - pos)
        intersection = BMO.intersect3d(shape, ray)
        differing += !same(intersection, bounded(shape, ray))
        isnothing(intersection) && continue
        p = BMO.position(ray) + intersection.t * BMO.direction(ray)
        hits += 1
        dmax = max(dmax, norm(p - center))
        next = Ray(p, random_direction(rng))
        differing += !same(BMO.intersect3d(shape, next), bounded(shape, next))
    end
    return hits, dmax, differing
end

@testset "SDF bounding spheres" begin
    list = instances()

    @testset "all concrete SDF types are covered" begin
        covered = Set(nameof(typeof(s)) for (_, s) in list)
        expected = [:BoxSDF, :CylinderSDF, :CutSphereSDF, :RingSDF, :RightAnglePrismSDF,
            :PolygonPrismSDF, :PlanoSurfaceSDF, :SphereSDF, :ConcaveSphericalSurfaceSDF,
            :ConvexSphericalSurfaceSDF, :ConvexCylinderSDF, :ConcaveCylinderSDF,
            :ConvexAsphericalSurfaceSDF, :ConcaveAsphericalSurfaceSDF, :AconvexCylinderSDF,
            :AconcaveCylinderSDF, :ConicSDF, :MeniscusLensSDF, :UnionSDF, :DifferenceSDF]
        @test length(expected) == 20
        @test issubset(expected, covered)
    end

    @testset "$label" for (i, (label, s)) in collect(enumerate(list))
        T = eltype(position(s))
        sphere = @inferred BMO.bounding_sphere_of(s)
        @test sphere isa BMO.SingleBoundingSphere{T}
        @test isfinite(sphere.radius) && sphere.radius > 0
        @test allocations(s) == 0

        # the sphere is fixed in the local frame of the shape
        local_center = BMO.transposed_orientation(s) * (sphere.pos - position(s))
        translate3d!(s, [0.3, -0.2, 0.5])
        rotate3d!(s, normalize([1.0, 2.0, 3.0]), 0.7)
        moved = @inferred BMO.bounding_sphere_of(s)
        center, r = moved.pos, moved.radius
        @test center ≈ position(s) + BMO.orientation(s) * local_center
        @test r ≈ sphere.radius

        rng = MersenneTwister(1000 + i)
        hits, dmax, differing = farthest_hit(rng, s, center, r)
        @test hits > N_RAYS / 100
        @test differing == 0                  # same hits and lengths as without the sphere
        @test dmax ≤ r + ENCLOSURE_TOL        # enclosure
        @test dmax ≥ TIGHTNESS * r            # tightness
    end

    @testset "composites" begin
        # an operand without a sphere
        @test BMO.bounding_sphere_of(BMO.SphereSDF(1e-2) + NoSphereSDF()) === BMO.NoBoundingSphere()
        @test BMO.bounding_sphere_of(NoSphereSDF() + BMO.SphereSDF(1e-2) + BMO.SphereSDF(1e-2)) === BMO.NoBoundingSphere()
        @test BMO.bounding_sphere_of(NoSphereSDF() - BMO.SphereSDF(1e-3)) === BMO.NoBoundingSphere()
        # a difference is bounded by its base, a tool without a sphere does not matter
        @test BMO.bounding_sphere_of(BMO.SphereSDF(1e-2) - NoSphereSDF()) == BMO.SingleBoundingSphere(Point3(0.0), 1e-2)

        # union of two spheres: the smallest sphere around both
        a, b = BMO.SphereSDF(1e-2), BMO.SphereSDF(2e-2)
        translate3d!(b, [0, 0.1, 0])
        u = a + b
        sphere = BMO.bounding_sphere_of(u)
        center, r = sphere.pos, sphere.radius
        @test r ≈ (0.1 + 1e-2 + 2e-2) / 2
        @test center ≈ Point3(0, 0.055, 0)
        # one sphere within the other
        c = BMO.SphereSDF(5e-2)
        translate3d!(c, [0, 0.01, 0])
        @test BMO.bounding_sphere_of(a + c) == BMO.SingleBoundingSphere(Point3(0, 0.01, 0), 5e-2)
        @test BMO.bounding_sphere_of(c + a) == BMO.SingleBoundingSphere(Point3(0, 0.01, 0), 5e-2)

        # an operand that is moved on its own: the sphere follows, in the frame of the composite
        zrotate3d!(u, π / 2)
        translate3d!(u, [1.0, 2.0, 3.0])
        translate3d!(b, [-0.1, 0, 0])
        sphere = BMO.bounding_sphere_of(u)
        center, r = sphere.pos, sphere.radius
        @test r ≈ (0.2 + 1e-2 + 2e-2) / 2
        @test center ≈ Point3(1.0 - 0.105, 2.0, 3.0)
        @test BMO.transposed_orientation(u) * (center - position(u)) ≈ Point3(0, 0.105, 0)

        base = BMO.SphereSDF(1e-2)
        translate3d!(base, [0, 0.1, 0])
        d = base - BMO.SphereSDF(1e-3)
        @test BMO.bounding_sphere_of(d) == BMO.SingleBoundingSphere(Point3(0, 0.1, 0), 1e-2)
        yrotate3d!(d, π / 3)
        translate3d!(d, [0.5, 0, 0])
        @test BMO.bounding_sphere_of(d).pos ≈ position(base)
        translate3d!(base, [0, 0, 0.2])
        @test BMO.bounding_sphere_of(d).pos ≈ position(base)
        @test BMO.bounding_sphere_of(d).radius == 1e-2
    end
end

end
