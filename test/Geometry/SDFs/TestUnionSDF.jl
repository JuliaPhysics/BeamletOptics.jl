module TestUnionSDF

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics

const BMO = BeamletOptics

"""
    NormalProbeSDF

Sphere-shaped probe with a specialized `normal3d` to prove that `normal3d(::UnionSDF, pos)` dispatches to the
operand's own method rather than falling back to the finite-difference normal.
"""
mutable struct NormalProbeSDF{T} <: BMO.AbstractSDF{T}
    pos::Point3{T}
    dir::BMO.SMatrix{3, 3, T, 9}
    transposed_dir::BMO.SMatrix{3, 3, T, 9}
    radius::T
end

NormalProbeSDF(r::T) where {T} = NormalProbeSDF{T}(
    Point3{T}(0), BMO.SMatrix{3, 3, T, 9}(I), BMO.SMatrix{3, 3, T, 9}(I), r)

const SENTINEL = Point3(0.0, 0.0, 1.0)

BMO.sdf(s::NormalProbeSDF, point) = norm(BMO._world_to_sdf(s, point)) - s.radius
BMO.normal3d(::NormalProbeSDF, ::Any) = SENTINEL

@testset "Union SDFs" begin
    @testset "Operator and flattening" begin
        a = BMO.SphereSDF(1.0)
        b = BMO.SphereSDF(0.5)
        c = BMO.SphereSDF(0.3)
        d = BMO.SphereSDF(0.2)

        # sdf + sdf
        ab = a + b
        @test ab isa BMO.UnionSDF
        @test BMO.operands(ab) === (a, b)

        # union + sdf appends, sdf + union prepends, union + union concatenates;
        # all of them stay flat instead of nesting unions
        @test BMO.operands(ab + c) === (a, b, c)
        @test BMO.operands(c + ab) === (c, a, b)
        @test BMO.operands(ab + (c + d)) === (a, b, c, d)

        # `shape += tool` accumulator idiom
        shape = a + b
        shape += c
        @test BMO.operands(shape) === (a, b, c)

        # mixed precision is rejected rather than silently promoted
        @test_throws MethodError BMO.SphereSDF(1.0f0) + BMO.SphereSDF(1.0)
    end

    @testset "Explicit constructor" begin
        a = BMO.SphereSDF(1.0)
        b = BMO.SphereSDF(0.5)
        c = BMO.CylinderSDF(0.5, 1.0)
        u = BMO.UnionSDF{Float64}(a, b, c)
        @test BMO.operands(u) === (a, b, c)
        @test position(u) == Point3(0.0, 0.0, 0.0)
        @test orientation(u) == I
        @test BMO.transposed_orientation(u) == I
    end

    @testset "sdf" begin
        # Two disjoint unit spheres: the union of exact SDFs is exact, i.e. the minimum of
        # the analytic distances to both spheres, inside and outside.
        c1 = Point3(-3.0, 0.0, 0.0)
        c2 = Point3(3.0, 0.0, 0.0)
        s1 = BMO.SphereSDF(1.0)
        s2 = BMO.SphereSDF(1.0)
        translate3d!(s1, c1)
        translate3d!(s2, c2)
        u = s1 + s2

        for x in range(-5, 5; length = 11), y in (-1.5, 0.0, 0.5), z in (-0.5, 0.0, 2.0)
            p = Point3(x, y, z)
            expected = min(norm(p - c1) - 1, norm(p - c2) - 1)
            @test BMO.sdf(u, p) ≈ expected
        end

        @test BMO.sdf(u, c1) ≈ -1                      # deep inside the first operand
        @test BMO.sdf(u, Point3(4.0, 0, 0)) ≈ 0 atol = 1e-12  # on the second surface
        @test BMO.sdf(u, Point3(0.0, 0, 0)) ≈ 2        # midway between both
    end

    @testset "normal3d" begin
        sphere = BMO.SphereSDF(1.0)
        probe = NormalProbeSDF(1.0)
        translate3d!(sphere, [-3.0, 0, 0])
        translate3d!(probe, [3.0, 0, 0])
        u = sphere + probe

        # closest operand is the probe -> its specialized method must be used
        @test BMO.normal3d(u, Point3(4.0, 0, 0)) == SENTINEL
        @test BMO.normal3d(u, Point3(3.0, 1.0, 0)) == SENTINEL
        # closest operand is the sphere -> regular outward normal
        @test isapprox(BMO.normal3d(u, Point3(-4.0, 0, 0)), Point3(-1.0, 0, 0); atol = 1e-6)
        @test isapprox(BMO.normal3d(u, Point3(-3.0, 1.0, 0)), Point3(0, 1.0, 0); atol = 1e-6)
    end

    @testset "thickness" begin
        a = BMO.SphereSDF(1.0)
        b = BMO.SphereSDF(0.5)
        cyl = BMO.CylinderSDF(0.5, 1.0)

        @test BMO.thickness(a + b) ≈ 3.0

        # operands without a `thickness` method are skipped
        @test !hasmethod(BMO.thickness, Tuple{typeof(cyl)})    # guard: keeps the test meaningful
        @test BMO.thickness(a + b + cyl) ≈ 3.0
        t0 = BMO.thickness(cyl + BMO.CylinderSDF(0.2, 1.0))
        @test t0 === 0.0
    end

    @testset "Kinematics (AbstractCompositeSDF)" begin
        c1 = BMO.CylinderSDF(0.5, 1.0)
        c2 = BMO.CylinderSDF(0.5, 1.0)
        translate3d!(c2, [2.0, 0, 0])
        u = c1 + c2

        # translation moves the pivot and every operand
        offset = Point3(1.0, 2.0, 3.0)
        translate3d!(u, offset)
        @test position(u) ≈ offset
        @test position(c1) ≈ offset
        @test position(c2) ≈ offset + Point3(2.0, 0, 0)

        # rotation happens around the union's own position (pivot), not the world origin
        # and not each operand's own origin
        axis = [0.0, 0.0, 1.0]
        θ = deg2rad(90)
        R = BMO.rotate3d(axis, θ)
        rotate3d!(u, axis, θ)
        @test position(u) ≈ offset
        @test position(c1) ≈ offset
        @test position(c2) ≈ offset + R * Point3(2.0, 0, 0)
        @test orientation(u) ≈ R
        @test orientation(c1) ≈ R
        @test orientation(c2) ≈ R

        # the rotated shape is where it should be: the axis of the second cylinder
        # (local y) is now at offset + R * (2, 0, 0)
        @test BMO.sdf(u, position(c2)) < 0
        @test BMO.sdf(u, offset + Point3(2.0, 0, 0)) > 0
    end
end

end # MODULE
