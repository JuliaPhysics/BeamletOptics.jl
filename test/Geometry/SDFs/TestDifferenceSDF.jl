module TestDifferenceSDF

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics

const BMO = BeamletOptics

@testset "Difference SDFs" begin
    @testset "Operator and chaining" begin
        a = BMO.SphereSDF(1.0)
        b = BMO.SphereSDF(0.5)
        c = BMO.SphereSDF(0.3)

        diff = a - b
        @test diff isa BMO.DifferenceSDF

        diff2 = diff - c
        @test diff2 isa BMO.DifferenceSDF
        @test length(diff2.tools) == 2
        @test diff2.base === a

        # `shape -= tool` accumulator idiom
        shape = a - b
        shape -= c
        @test shape isa BMO.DifferenceSDF
        @test length(shape.tools) == 2
        @test shape.base === a
    end

    @testset "Union/difference interoperability" begin
        # `+` and `-` share precedence and are left-associative, so every mixed chain is
        # evaluated strictly left to right. Layout along the x-axis:
        #   s1: r = 1   at 0      -> spans [-1, 1]
        #   s2: r = 1   at 1.5    -> spans [0.5, 2.5]
        #   s3: r = 0.5 at 0.75   -> spans [0.25, 1.25], sits in the s1 ∩ s2 lens
        #   s4: r = 0.3 at -0.8   -> spans [-1.1, -0.5], inside s1 only
        function spheres()
            s1 = BMO.SphereSDF(1.0)
            s2 = BMO.SphereSDF(1.0)
            s3 = BMO.SphereSDF(0.5)
            s4 = BMO.SphereSDF(0.3)
            translate3d!(s2, [1.5, 0, 0])
            translate3d!(s3, [0.75, 0, 0])
            translate3d!(s4, [-0.8, 0, 0])
            return s1, s2, s3, s4
        end
        s1, s2, s3, s4 = spheres()

        grid = [Point3(x, y, z) for x in -1.5:0.25:3.0, y in (-0.5, 0.0, 0.4), z in (0.0, 0.3)]
        d(s, p) = BMO.sdf(s, p)
        # composite sdf must equal the boolean formula applied to the operand sdfs
        matches(expr, ref) = all(p -> BMO.sdf(expr, p) ≈ ref(p), grid)

        @testset "s1 + s2 - s3 == (s1 ∪ s2) \\ s3" begin
            e = s1 + s2 - s3
            @test e isa BMO.DifferenceSDF
            @test e.base isa BMO.UnionSDF
            @test BMO.operands(e.base) === (s1, s2)
            @test e.tools === (s3,)
            @test matches(e, p -> max(min(d(s1, p), d(s2, p)), -d(s3, p)))
        end

        @testset "s1 - s2 + s3 == (s1 \\ s2) ∪ s3" begin
            e = s1 - s2 + s3
            @test e isa BMO.UnionSDF
            inner, added = BMO.operands(e)
            @test inner isa BMO.DifferenceSDF
            @test inner.base === s1
            @test inner.tools === (s2,)
            @test added === s3
            @test matches(e, p -> min(max(d(s1, p), -d(s2, p)), d(s3, p)))
        end

        @testset "Order changes the shape" begin
            carve_last = s1 + s2 - s3
            add_last = s1 - s2 + s3
            # centre of s3: carved away vs. added back
            @test BMO.sdf(carve_last, Point3(0.75, 0, 0)) > 0
            @test BMO.sdf(add_last, Point3(0.75, 0, 0)) < 0
            # centre of s2: added vs. subtracted (and outside s1 anyway)
            @test BMO.sdf(carve_last, Point3(1.5, 0, 0)) < 0
            @test BMO.sdf(add_last, Point3(1.5, 0, 0)) > 0
        end

        @testset "Parentheses: s1 + (s2 - s3) vs (s1 + s2) - s3" begin
            grouped = s1 + (s2 - s3)
            @test grouped isa BMO.UnionSDF
            @test BMO.operands(grouped)[1] === s1
            @test BMO.operands(grouped)[2] isa BMO.DifferenceSDF
            @test matches(grouped, p -> min(d(s1, p), max(d(s2, p), -d(s3, p))))
            # x = 0.4 lies in s1 and s3 but not in s2: only the left-to-right chain carves it
            @test BMO.sdf(grouped, Point3(0.4, 0, 0)) < 0
            @test BMO.sdf(s1 + s2 - s3, Point3(0.4, 0, 0)) > 0
        end

        @testset "Parentheses: s1 - (s2 + s3) == s1 - s2 - s3" begin
            # A \ (B ∪ C) = (A \ B) \ C: same shape, different structure
            grouped = s1 - (s2 + s3)
            chained = s1 - s2 - s3
            @test grouped.base === s1
            @test length(grouped.tools) == 1
            @test grouped.tools[1] isa BMO.UnionSDF      # a union tool is not flattened
            @test chained.tools === (s2, s3)               # a difference chain is
            @test all(p -> BMO.sdf(grouped, p) ≈ BMO.sdf(chained, p), grid)
        end

        @testset "Difference minus union" begin
            e = (s1 - s2) - (s3 + s4)
            @test e isa BMO.DifferenceSDF
            @test e.base === s1
            @test e.tools[1] === s2
            @test e.tools[2] isa BMO.UnionSDF
            @test matches(e, p -> max(d(s1, p), -d(s2, p), -min(d(s3, p), d(s4, p))))
        end

        @testset "Union of differences" begin
            e = (s1 - s4) + (s2 - s3)
            @test e isa BMO.UnionSDF
            @test all(o -> o isa BMO.DifferenceSDF, BMO.operands(e))
            @test matches(e, p -> min(max(d(s1, p), -d(s4, p)), max(d(s2, p), -d(s3, p))))
        end

        @testset "Material added after a subtraction is kept" begin
            a = BMO.SphereSDF(1.0)
            b = BMO.SphereSDF(1.0)
            translate3d!(b, [3.0, 0, 0])
            c = BMO.SphereSDF(1.0)
            translate3d!(c, [1.5, 0, 0])
            e = BMO.SphereSDF(1.0)
            translate3d!(e, [1.5, 0, 0])

            result = a + b - c + e
            @test result isa BMO.UnionSDF
            @test BMO.operands(result)[1] isa BMO.DifferenceSDF
            # probe inside both c and e (same location): must not be carved away
            @test BMO.sdf(result, Point3(1.5, 0, 0)) < 0
        end

        @testset "normal3d through nested composites" begin
            # separate layout so every surface point has a single unambiguous owner
            base = BMO.SphereSDF(1.0)
            tool = BMO.SphereSDF(0.5)
            extra = BMO.SphereSDF(0.3)
            translate3d!(tool, [1.0, 0, 0])
            translate3d!(extra, [-2.0, 0, 0])
            e = base - tool + extra

            # carved wall of `tool` inside `base`: outward normal is the flipped tool normal
            @test isapprox(BMO.normal3d(e, Point3(0.5, 0, 0)), Point3(1.0, 0, 0); atol = 1e-6)
            # remaining outer wall of `base`
            @test isapprox(BMO.normal3d(e, Point3(0, 1.0, 0)), Point3(0, 1.0, 0); atol = 1e-6)
            # surface of the union's second operand
            @test isapprox(BMO.normal3d(e, Point3(-2.3, 0, 0)), Point3(-1.0, 0, 0); atol = 1e-6)
        end

        @testset "Kinematics of nested composites" begin
            # Moving the outermost composite must move the whole nested shape rigidly:
            # sdf_after(R * p + t) == sdf_before(p) for every expression.
            axis = [0.0, 1.0, 1.0] / sqrt(2)
            θ = deg2rad(35)
            R = BMO.rotate3d(axis, θ)
            t = Point3(0.3, -1.0, 2.0)
            builders = (
                (a, b, c, _) -> a + b - c,
                (a, b, c, _) -> a - b + c,
                (a, b, c, _) -> a - (b + c),
                (a, b, c, _) -> a + (b - c),
                (a, b, c, e) -> (a - b) - (c + e),
            )
            for build in builders
                shape = build(spheres()...)
                before = [BMO.sdf(shape, p) for p in grid]
                rotate3d!(shape, axis, θ)
                translate3d!(shape, t)
                after = [BMO.sdf(shape, R * p + t) for p in grid]
                @test after ≈ before
            end
        end
    end

    @testset "sdf signs" begin
        base = BMO.SphereSDF(2.0)
        tool = BMO.SphereSDF(1.0)
        translate3d!(tool, [1.0, 0, 0])
        diff = base - tool

        # inside base, outside all tools
        p1 = Point3(-1.5, 0.0, 0.0)
        @test BMO.sdf(diff, p1) < 0

        # inside a tool
        p2 = Point3(1.0, 0.0, 0.0)
        @test BMO.sdf(diff, p2) > 0

        # far exterior matches the base
        p3 = Point3(10.0, 0.0, 0.0)
        @test BMO.sdf(diff, p3) ≈ BMO.sdf(base, p3)
    end

    @testset "Conservativeness" begin
        # A cylinder minus a longer coaxial bore is a tube. Its exact signed distance is
        # the 2D distance to the cross-section rectangle [r,R] x [-H,H] in the (ρ, y)
        # half-plane, revolved about the y-axis (CylinderSDF is y-aligned and takes a
        # HALF-height). The bore overshoots both caps so there are no coincident faces.
        R = 2.0
        r = 0.7
        H = 5.0
        outer = BMO.CylinderSDF(R, H)
        inner = BMO.CylinderSDF(r, H + 1.0)
        tube = outer - inner

        function d_tube(ρ, y)
            qx = abs(ρ - (r + R) / 2) - (R - r) / 2
            qy = abs(y) - H
            return norm(max.(Point2(qx, qy), 0.0)) + min(max(qx, qy), 0.0)
        end

        # Grid crosses both radial rims AND both end caps, so it includes the concave
        # corner regions where the max()-based bound is weakest — not just the axial
        # midplane, where the composite happens to be exact.
        worst = -Inf   # largest (sdf - true): must stay <= 0, an overestimate tunnels
        slack = Inf    # smallest (sdf - true): proves the grid reaches the degraded region
        for ρ in range(0.0, 1.5R; length = 60), y in range(-1.5H, 1.5H; length = 60)
            gap = BMO.sdf(tube, Point3(ρ, y, 0.0)) - d_tube(ρ, y)
            worst = max(worst, gap)
            slack = min(slack, gap)
        end
        # The property both ray marchers depend on: never an overestimate.
        @test worst <= 1e-12
        # Guard against the test going vacuous if the geometry is ever made exact.
        @test slack < -1e-3
    end

    @testset "normal3d: bore wall vs. outer wall" begin
        R = 2.0
        r = 0.7
        H = 5.0
        outer = BMO.CylinderSDF(R, H)
        inner = BMO.CylinderSDF(r, H)
        diff = outer - inner

        θ = deg2rad(30)
        zrotate3d!(diff, θ)

        # outer wall
        p_outer = Point3(R * cos(θ), R * sin(θ), 0.0)
        n_outer = BMO.normal3d(diff, p_outer)
        n_outer_expected = Point3(cos(θ), sin(θ), 0.0)
        @test isapprox(n_outer, n_outer_expected; atol = 1e-6)

        # bore wall: outward normal of the result is the negated tool normal
        p_inner = Point3(r * cos(θ), r * sin(θ), 0.0)
        n_inner = BMO.normal3d(diff, p_inner)
        n_inner_expected = Point3(-cos(θ), -sin(θ), 0.0)
        @test isapprox(n_inner, n_inner_expected; atol = 1e-6)
    end

    @testset "Kinematics parity" begin
        c1 = BMO.CylinderSDF(1.0, 2.0)
        c2 = BMO.CylinderSDF(0.5, 2.0)
        translate3d!(c2, [1.0, 0, 0])
        u = c1 + c2

        c1b = BMO.CylinderSDF(1.0, 2.0)
        c2b = BMO.CylinderSDF(0.5, 2.0)
        translate3d!(c2b, [1.0, 0, 0])
        diff = c1b - c2b

        offset = [2.0, -1.0, 0.5]
        axis = [0.0, 0.0, 1.0]
        θ = deg2rad(40)

        translate3d!(u, offset)
        rotate3d!(u, axis, θ)
        translate3d!(diff, offset)
        rotate3d!(diff, axis, θ)

        u_ops = BMO.operands(u)
        diff_ops = BMO.operands(diff)
        for i in 1:2
            @test position(u_ops[i]) ≈ position(diff_ops[i])
            @test orientation(u_ops[i]) ≈ orientation(diff_ops[i])
        end
    end

    @testset "align3d!, reset_rotation3d! and reset_translation3d! propagation" begin
        c1 = BMO.CylinderSDF(1.0, 2.0)
        c2 = BMO.CylinderSDF(0.3, 1.0)
        translate3d!(c2, [0.0, 0.5, 0.0])
        diff = c1 - c2

        # 1. align3d!
        target_axis = [1.0, 0.0, 0.0]
        align3d!(diff, target_axis)
        @test isapprox(orientation(diff)[:, 2], target_axis; atol = 1e-12)
        @test isapprox(orientation(diff.base)[:, 2], target_axis; atol = 1e-12)
        @test isapprox(orientation(diff.tools[1])[:, 2], target_axis; atol = 1e-12)
        # Operand position should rotate around composite pivot
        @test isapprox(position(diff.tools[1]) - position(diff), Point3(0.0, 0.0, -0.5); atol = 1e-12) ||
              isapprox(norm(position(diff.tools[1]) - position(diff)), 0.5; atol = 1e-12)

        # 2. reset_translation3d!
        translate3d!(diff, [2.0, 3.0, -1.0])
        @test position(diff) ≈ Point3(2.0, 3.0, -1.0)
        reset_translation3d!(diff)
        @test position(diff) == Point3(0.0, 0.0, 0.0)

        # 3. reset_rotation3d!
        reset_rotation3d!(diff)
        @test isapprox(orientation(diff), Matrix(1.0I, 3, 3); atol = 1e-12)
        @test isapprox(orientation(diff.base), Matrix(1.0I, 3, 3); atol = 1e-12)
        @test isapprox(orientation(diff.tools[1]), Matrix(1.0I, 3, 3); atol = 1e-12)
        # Initial relative offset was along local y: [0, 0.5, 0]
        @test isapprox(position(diff.tools[1]), Point3(0.0, 0.5, 0.0); atol = 1e-12)
    end
end

end # MODULE
