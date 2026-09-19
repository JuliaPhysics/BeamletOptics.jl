module TestSDFs

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics
using Random

const BMO = BeamletOptics

@testset "SDFs" begin
    @testset "Testing type definitions" begin
        @test isdefined(BMO, :AbstractSDF)
        @test isdefined(BMO, :SphereSDF)
        @test isdefined(BMO, :CylinderSDF)
        @test isdefined(BMO, :CutSphereSDF)
        @test isdefined(BMO, :ThinLensSDF)
    end

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

        @testset "Mixed chain nesting" begin
            a = BMO.SphereSDF(1.0)
            b = BMO.SphereSDF(1.0)
            translate3d!(b, [3.0, 0, 0])
            c = BMO.SphereSDF(1.0)
            translate3d!(c, [1.5, 0, 0])
            d = BMO.SphereSDF(1.0)
            translate3d!(d, [1.5, 0, 0])

            result = a + b - c + d
            @test result isa BMO.UnionSDF
            @test BMO.operands(result)[1] isa BMO.DifferenceSDF

            # probe point inside both c and d (same location): must not be carved away
            probe = Point3(1.5, 0, 0)
            @test BMO.sdf(result, probe) < 0
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
    end

    @testset "Conic SDF" begin
        @testset "Parabola equivalence" begin
            Random.seed!(1)
            f = 0.1
            x_off = 0.02
            D = 0.05
            t = 0.03
            c = BMO.ConicSDF(2f, -1, x_off, D, t)
            o = BMO.OffAxisParaboloidSDF(f, x_off, D, t)
            for _ in 1:100
                p = Point3(0.1 * (rand() - 0.5), 0.1 * (rand() - 0.5), 0.1 * (rand() - 0.5))
                @test BMO.sdf(c, p) == BMO.sdf(o, p)
            end
        end

        @testset "Sphere equivalence (k = 0)" begin
            R = 0.5
            D = 0.2
            t = 0.05
            cs = BMO.ConicSDF(R, 0, 0, D, t)
            sm = BMO.shape(SphericalMirror(R, t, D))
            for ρ in range(0, D / 2; length = 50)
                y = -(R - sqrt(R^2 - ρ^2))
                p = Point3(ρ, y, 0.0)
                @test abs(BMO.sdf(cs, p)) < 1e-12
                @test abs(BMO.sdf(sm, p)) < 1e-12
                n1 = BMO.normal3d(cs, p)
                n2 = BMO.normal3d(sm, p)
                @test isapprox(n1, n2; atol = 1e-6)
            end
        end

        @testset "Conservativeness" begin
            # Three segments spanning the conic family: parabolic, prolate ellipsoidal and
            # convex hyperbolic. Precompute a dense (ρ, θ) sampling of the front face in
            # local coordinates (the SDFs are left at the origin so no world transform is
            # needed) and bound the true distance from below by the minimum distance to
            # that grid. The grid overestimates the true distance to the surface by at most
            # half the grid spacing (~ (D/2)/600 ≈ 4e-5 m here), which the -1e-3 vacuity
            # guard comfortably clears.
            D = 0.05
            segments = [
                (0.2, -1.0, 0.0),     # parabolic
                (0.2, -0.5, 0.1),     # ellipsoidal, off-axis
                (-0.15, -2.5, 0.0),   # convex hyperbolic
            ]
            nρ, nθ = 600, 720
            Random.seed!(2)
            for (R, k, x_off) in segments
                t = BMO._conic_auto_thickness(R, k, x_off, D, Float64)
                s = BMO.ConicSDF(R, k, x_off, D, t)
                Z_off = BMO._conic_sag(x_off, R, k)
                r_max = D / 2

                gx = Vector{Float64}(undef, nρ * nθ)
                gy = Vector{Float64}(undef, nρ * nθ)
                gz = Vector{Float64}(undef, nρ * nθ)
                idx = 1
                for ρ in range(0, r_max; length = nρ),
                    θ in range(0, 2π; length = nθ + 1)[1:(end - 1)]

                    x = ρ * cos(θ)
                    z = ρ * sin(θ)
                    r_p = sqrt((x + x_off)^2 + z^2)
                    y = -(BMO._conic_sag(r_p, R, k) - Z_off)
                    gx[idx] = x
                    gy[idx] = y
                    gz[idx] = z
                    idx += 1
                end

                worst = -Inf
                slack = Inf
                for _ in 1:2000
                    ρ = 0.9 * r_max * rand()
                    θ = 2π * rand()
                    x = ρ * cos(θ)
                    z = ρ * sin(θ)
                    r_p = sqrt((x + x_off)^2 + z^2)
                    y_surf = -(BMO._conic_sag(r_p, R, k) - Z_off)
                    y = y_surf + (D / 4) * (2 * rand() - 1)
                    p = Point3(x, y, z)

                    dvals = @. sqrt((gx - x)^2 + (gy - y)^2 + (gz - z)^2)
                    dgrid = minimum(dvals)
                    sdfval = BMO.sdf(s, p)

                    gap = abs(sdfval) - dgrid
                    worst = max(worst, gap)
                    slack = min(slack, gap)
                end
                # The property both ray marchers depend on: never an overestimate.
                @test worst <= 1e-12
                # Guard against the test going vacuous if the geometry is ever made exact.
                @test slack < -1e-3
            end
        end

        @testset "Domain errors" begin
            @test_throws ArgumentError BMO.ConicSDF(0.1, 0.5, 0.0, 0.2, 0.05)
            @test_throws ArgumentError BMO.ConicSDF(0.1, 0.5, 0.09, 0.02, 0.05)
            @test_throws ArgumentError BMO.ConicSDF(0.0, -1.0, 0.0, 0.05, 0.01)
            # k <= -1 accepts an arbitrarily large aperture
            @test BMO.ConicSDF(0.1, -1.0, 0.0, 1000.0, 0.01) isa BMO.ConicSDF
            @test BMO.ConicSDF(0.1, -2.0, 0.0, 1000.0, 0.01) isa BMO.ConicSDF
        end
    end
end

end # MODULE