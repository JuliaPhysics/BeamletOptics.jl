module TestPrisms

using BeamletOptics
using LinearAlgebra
using Random
using Test

const BMO = BeamletOptics
const mm = 1e-3

@testset "Prisms" begin
    @testset "PolygonPrismSDF" begin
        # unit square, 2 high: exact distances inside, outside, on the faces and at the edges/corners
        sq = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
        for verts in (sq, reverse(sq))
            s = BMO.PolygonPrismSDF(verts, 2.0)
            @test s.vertices == BMO.Point2.(sq) # normalized to counter-clockwise
            @test BMO.sdf(s, [0, 0, 0.0]) ≈ -1
            @test BMO.sdf(s, [0.5, 0, 0.0]) ≈ -0.5
            @test BMO.sdf(s, [0, 0, 0.75]) ≈ -0.25
            @test abs(BMO.sdf(s, [1, 0.3, 0.2])) < 1e-12
            @test abs(BMO.sdf(s, [0.2, -1, 1.0])) < 1e-12
            @test BMO.sdf(s, [3, 0, 0.0]) ≈ 2
            @test BMO.sdf(s, [2, 2, 0.0]) ≈ sqrt(2)
            @test BMO.sdf(s, [2, 2, 3.0]) ≈ sqrt(1 + 1 + 4)
            @test BMO.sdf(s, [0, 0, 3.0]) ≈ 2
            @test BMO.normal3d(s, [1, 0.3, 0.2]) ≈ [1, 0, 0]
            @test BMO.normal3d(s, [0.1, 0.2, 1.0]) ≈ [0, 0, 1]
        end
        # triangle: the distance is exact (compare with the brute-force distance to the faces)
        tri = BMO.PolygonPrismSDF([(0.0, 0.0), (4.0, 0.0), (1.0, 3.0)], 1.0)
        @test BMO.sdf(tri, [-1, -1, 0.0]) ≈ sqrt(2)
        @test BMO.sdf(tri, [2.5, 1.5, 0.0]) ≈ 0 atol = 1e-12 # on the hypotenuse
        @test BMO.sdf(tri, [3, 3, 0.0]) ≈ norm([3, 3] - [2.5, 1.5] + [0, 0]) atol = 1.5 # outside, finite
        @test BMO.sdf(tri, [3, 3, 0.0]) > 0
        @test BMO.sdf(tri, [1.5, 1, 0.0]) < 0
        # exterior distance equals the brute-force distance to the nearest edge segment (random convex polygons)
        segdist(q, a, b) = (e = b - a; norm(q - (a + clamp(dot(q - a, e) / dot(e, e), 0, 1) * e)))
        rng = MersenneTwister(1)
        let nchecked = 0
            for _ in 1:100
                θ = sort(2π * rand(rng, rand(rng, 3:9)))
                V = [(cos(t), sin(t)) .* (0.2 + 3rand(rng)) for t in θ]
                s = try BMO.PolygonPrismSDF(V, 1.0) catch; continue end
                P = s.vertices
                n = length(P)
                for _ in 1:100
                    q = BMO.Point2(4randn(rng), 4randn(rng))
                    ref = minimum(i -> segdist(q, P[i], P[mod1(i + 1, n)]), 1:n)
                    inside = all(i -> (e = P[mod1(i + 1, n)] - P[i]; w = q - P[i]; e[1] * w[2] - e[2] * w[1] >= 0), 1:n)
                    inside && continue
                    nchecked += 1
                    @test BMO.sdf(s, [q[1], q[2], 0.0]) ≈ ref atol = 1e-12
                end
            end
            @test nchecked > 1000
        end
        # a polygon with a reflex corner (where a vertex-only fallback would be wrong) is rejected
        @test_throws ArgumentError BMO.PolygonPrismSDF([(0.0, 0.0), (5.0, 0.0), (4.0, 1.0), (3.0, 5.0), (0.0, 4.0)], 1.0)
        # the degeneracy test is translation invariant, also in Float32
        for off in (0f0, 1f0, 100f0)
            @test BMO.PolygonPrismSDF([(off, off), (off + 5f-4, off), (off + 5f-4, off + 5f-4), (off, off + 5f-4)], 1f0) isa BMO.PolygonPrismSDF{Float32}
        end
        # kinematics: the SDF follows the pose
        BMO.translate3d!(tri, [1, 1, 1.0])
        BMO.zrotate3d!(tri, deg2rad(30))
        @test BMO.sdf(tri, [1, 1, 1.0] + BMO.orientation(tri) * [1.5, 1, 0]) < 0
        @test BMO.thickness(BMO.PolygonPrismSDF(sq, 1.0)) ≈ 2

        # validation
        @test_throws ArgumentError BMO.PolygonPrismSDF([(0.0, 0.0), (1.0, 0.0)], 1.0)
        @test_throws ArgumentError BMO.PolygonPrismSDF(sq, -1.0)
        @test_throws ArgumentError BMO.PolygonPrismSDF([(0.0, 0.0), (2.0, 0.0), (1.0, 0.2), (2.0, 2.0), (0.0, 2.0)], 1.0) # reflex
        @test_throws ArgumentError BMO.PolygonPrismSDF([(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (1.0, 1.0)], 1.0) # collinear
        @test_throws ArgumentError BMO.PolygonPrismSDF([(0.0, 0.0), (1.0, 1.0), (1.0, 0.0), (0.0, 1.0)], 1.0) # bow tie
        @test_throws ArgumentError BMO.PolygonPrismSDF([(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)], 1.0) # zero area
    end

    @testset "Constructors" begin
        @test Prism([(0, 0), (1, 0), (0, 1)], 1, 1.5) isa Prism
        eq = EquilateralPrism(20mm, 10mm, 1.5)
        @test eq isa Prism
        V = eq.shape.vertices
        @test sum(V) / 3 ≈ zeros(2) atol = 1e-15 # centroid at the origin
        @test V[3][2] > 0 && V[3][1] == 0 # apex along +y
        @test all(norm(V[i] - V[mod1(i + 1, 3)]) ≈ 20mm for i in 1:3)
        dove = DovePrism(50mm, 10mm, 10mm, 1.5)
        @test dove isa Prism
        @test BMO.thickness(dove) ≈ 50mm
        @test_throws ArgumentError DovePrism(10mm, 10mm, 10mm, 1.5)
        # a Dove prism does not deviate a beam on its axis (n = 1.5, parallel to the base)
        system = System([dove])
        beam = Beam([-1mm, -40mm, 0], [0, 1, 0], 1e-6)
        solve_system!(system, beam)
        @test length(beam.rays) == 4 # refraction, total internal reflection, refraction
        @test BMO.direction(beam.rays[end]) ≈ [0, 1, 0] atol = 1e-9
    end

    @testset "Equilateral prism at minimum deviation" begin
        # N-SF11
        nsf11 = SellmeierEquation(1.73759695, 0.313747346, 1.89878101,
            0.013188707, 0.0623068142, 155.23629)
        for λ in (450e-9, 589e-9, 800e-9)
            n = nsf11(λ)
            α = deg2rad(60)
            δ_min = 2asin(n * sin(α / 2)) - α
            θ = asin(n * sin(α / 2)) - α / 2 # incidence offset from the base-parallel inner ray
            side = 40mm
            prism = EquilateralPrism(side, 20mm, nsf11)
            r = side * sqrt(3) / 6
            # inner ray is parallel to the base, entering at the midpoint of the left face
            d_in = [cos(θ), sin(θ), 0]
            hit = [-side / 4, r / 2, 0]
            beam = Beam(hit - 30mm * d_in, d_in, λ)
            solve_system!(System([prism]), beam)
            @test length(beam.rays) == 3
            @test isempty(beam.children)
            r1, r2, r3 = beam.rays
            # rays enter and leave exactly on the faces
            p1, p2 = BMO.position(r2), BMO.position(r3)
            @test BMO.sdf(prism.shape, p1) ≈ 0 atol = 1e-9
            @test BMO.sdf(prism.shape, p2) ≈ 0 atol = 1e-9
            @test p1 ≈ hit atol = 1e-9
            @test p2 ≈ [side / 4, r / 2, 0] atol = 1e-9
            @test BMO.refractive_index(r1) == 1
            @test BMO.refractive_index(r2) ≈ n
            @test BMO.refractive_index(r3) == 1
            @test BMO.direction(r2) ≈ [1, 0, 0] atol = 1e-9
            # deviation matches the analytic minimum deviation
            δ = acos(clamp(dot(BMO.direction(r1), BMO.direction(r3)), -1, 1))
            @test δ ≈ δ_min atol = 1e-8
            @test BMO.direction(r3)[2] < 0 # deflected towards the base
        end
    end
end

end # module
