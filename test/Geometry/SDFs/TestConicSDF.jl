module TestConicSDF

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics
using Random

const BMO = BeamletOptics

@testset "Conic SDF" begin
    @testset "OffAxisParaboloidSDF alias" begin
        # `OffAxisParaboloidSDF` is a thin constructor alias for `ConicSDF(2f, -1, ...)`, so
        # comparing parameters is sufficient; the sdf itself is covered below.
        f = 0.1
        x_off = 0.02
        D = 0.05
        t = 0.03
        o = BMO.OffAxisParaboloidSDF(f, x_off, D, t)
        @test o isa BMO.ConicSDF
        @test o.f == f
        @test o.k == -1
        @test o.x_off == x_off
        @test o.diameter == D
        @test o.thickness == t
    end

    @testset "Bounding sphere (on-axis paraboloid)" begin
        f = 0.1
        r_max = 0.03
        thickness = 0.02
        sag_max = r_max^2 / (4f)

        sdf = BMO.ConicSDF(2f, -1, 0, 2r_max, thickness)
        center, r = BMO.bounding_sphere(sdf)

        @test center ≈ Point3(0, (thickness - sag_max) / 2, 0)
        @test r ≈ sqrt(r_max^2 + ((thickness + sag_max) / 2)^2) + 0.05

        # bounding_box must transform the sphere into a symmetric, aperture-covering box
        xmin, xmax, ymin, ymax, zmin, zmax = BMO.bounding_box(sdf)
        @test xmin ≈ -xmax
        @test zmin ≈ -zmax
        @test xmax ≈ zmax # rotationally symmetric about the y-axis since x_off = 0
        @test xmax >= r_max
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

end # MODULE
