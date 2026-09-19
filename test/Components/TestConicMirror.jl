module TestConicMirror

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics
const mm = 1e-3

@testset "Conic Mirrors" begin
    @testset "Parabola via ConicMirror" begin
        f = 100mm
        D = 60mm
        m = ConicMirror(2f, -1, D)

        @test BMO.shape(ConicMirror(2f, -1, D)).thickness ≈ BMO.shape(ParabolicMirror(f, D)).thickness

        F = [0, -f, 0]
        for (hx, hz) in [(0, 0), (5mm, 0), (0, 20mm), (18mm, -18mm), (29mm, 0)]
            beam = Beam(Ray([hx, -50mm, hz], [0, 1, 0]))
            solve_system!(StaticSystem([m]), beam)
            @test length(BMO.rays(beam)) == 2
            r = BMO.rays(beam)[end]
            p = position(r)
            dir = BMO.direction(r)
            @test norm(cross(F - p, dir)) < 1e-9
            @test dot(F - p, dir) > 0
        end
    end

    @testset "Ellipsoid, on-axis" begin
        s1 = 300mm
        s2 = 150mm
        D = 40mm
        m = EllipsoidalMirror(s1, s2, D)
        sh = BMO.shape(m)

        @test BMO.radius(sh) ≈ 2 * s1 * s2 / (s1 + s2)
        @test -1 < sh.k < 0

        F1 = [0, -s1, 0]
        F2 = [0, -s2, 0]
        for ang in range(0, 2π; length = 7)[1:6]
            rr = 0.4 * D / 2
            aim = [rr * cos(ang), 0, rr * sin(ang)]
            dir0 = normalize(aim - F1)
            beam = Beam(Ray(F1, dir0))
            solve_system!(StaticSystem([m]), beam)
            @test length(BMO.rays(beam)) == 2
            r = BMO.rays(beam)[end]
            p = position(r)
            dir = BMO.direction(r)
            @test norm(cross(F2 - p, dir)) < 1e-9
            @test dot(F2 - p, dir) > 0
        end
    end

    @testset "Ellipsoid, off-axis" begin
        s1 = 300mm
        s2 = 150mm
        x_off = 60mm
        D = 30mm
        m = OffAxisEllipsoidalMirror(s1, s2, x_off, D)
        sh = BMO.shape(m)
        R = BMO.radius(sh)
        k = sh.k
        Z_off = BMO._conic_sag(x_off, R, k)
        vertex = [-x_off, Z_off, 0]
        F1 = vertex + [0, -s1, 0]
        F2 = vertex + [0, -s2, 0]

        for ang in range(0, 2π; length = 7)[1:6]
            rr = 0.4 * D / 2
            x = rr * cos(ang)
            z = rr * sin(ang)
            r_p = sqrt((x + x_off)^2 + z^2)
            y_true = -(BMO._conic_sag(r_p, R, k) - Z_off)
            aim = [x, y_true, z]
            dir0 = normalize(aim - F1)
            beam = Beam(Ray(F1, dir0))
            solve_system!(StaticSystem([m]), beam)
            @test length(BMO.rays(beam)) == 2
            r = BMO.rays(beam)[end]
            p = position(r)
            dir = BMO.direction(r)
            @test norm(cross(F2 - p, dir)) < 1e-9
            @test dot(F2 - p, dir) > 0
        end
    end

    @testset "Hyperboloid, convex" begin
        # Aim each incoming ray through a known point Q on the true surface, extended
        # backward away from the virtual focus, so it approaches from in front of the
        # mirror and converges towards the virtual focus at Q.
        b = 40mm
        a = 200mm
        D = 25mm
        m = HyperbolicMirror(-b, a, D)
        sh = BMO.shape(m)

        @test BMO.radius(sh) < 0
        @test sh.k < -1

        R = BMO.radius(sh)
        k = sh.k
        Fv = [0, b, 0]
        Fr = [0, -a, 0]
        for ang in range(0, 2π; length = 7)[1:6]
            rr = 0.4 * D / 2
            x = rr * cos(ang)
            z = rr * sin(ang)
            r_p = sqrt(x^2 + z^2)
            y_true = -BMO._conic_sag(r_p, R, k)
            Q = [x, y_true, z]
            dirq = normalize(Fv - Q)
            p0 = Q .- 10 .* dirq
            beam = Beam(Ray(p0, dirq))
            solve_system!(StaticSystem([m]), beam)
            @test length(BMO.rays(beam)) == 2
            r = BMO.rays(beam)[end]
            p = position(r)
            dir = BMO.direction(r)
            @test norm(cross(Fr - p, dir)) < 1e-9
            @test dot(Fr - p, dir) > 0
        end

        # OffAxis variant, vertex-shifted focus positions. The segment's local aperture is
        # centred on the origin while the parent vertex (and virtual focus) is offset by
        # x_off2 along x, so Q is computed from the parent radius r_p, not the local ρ.
        x_off2 = 20mm
        D2 = 20mm
        m2 = OffAxisHyperbolicMirror(-b, a, x_off2, D2)
        sh2 = BMO.shape(m2)
        R2 = BMO.radius(sh2)
        k2 = sh2.k
        Z_off2 = BMO._conic_sag(x_off2, R2, k2)
        vertex2 = [-x_off2, Z_off2, 0]
        Fv2 = vertex2 + [0, b, 0]
        Fr2 = vertex2 + [0, -a, 0]

        for ang in range(0, 2π; length = 7)[1:6]
            rr = 0.4 * D2 / 2
            x = rr * cos(ang)
            z = rr * sin(ang)
            r_p = sqrt((x + x_off2)^2 + z^2)
            y_true = -(BMO._conic_sag(r_p, R2, k2) - Z_off2)
            Q = [x, y_true, z]
            dirq = normalize(Fv2 - Q)
            p0 = Q .- 10 .* dirq
            beam = Beam(Ray(p0, dirq))
            solve_system!(StaticSystem([m2]), beam)
            @test length(BMO.rays(beam)) == 2
            r = BMO.rays(beam)[end]
            p = position(r)
            dir = BMO.direction(r)
            @test norm(cross(Fr2 - p, dir)) < 1e-9
            @test dot(Fr2 - p, dir) > 0
        end
    end

    @testset "Kinematics" begin
        s1 = 300mm
        s2 = 150mm
        D = 40mm
        m = EllipsoidalMirror(s1, s2, D)
        F1 = [0, -s1, 0]
        F2 = [0, -s2, 0]

        θ = π / 7
        yrotate3d!(m, θ)
        offset = [0.1, -0.05, 0.02]
        translate3d!(m, offset)

        Rm = BMO.rotate3d([0, 1, 0], θ)
        F1t = Rm * F1 + offset
        F2t = Rm * F2 + offset

        for ang in range(0, 2π; length = 7)[1:6]
            rr = 0.4 * D / 2
            aim_local = [rr * cos(ang), 0, rr * sin(ang)]
            aim = Rm * aim_local + offset
            dir0 = normalize(aim - F1t)
            beam = Beam(Ray(F1t, dir0))
            solve_system!(StaticSystem([m]), beam)
            @test length(BMO.rays(beam)) == 2
            r = BMO.rays(beam)[end]
            p = position(r)
            dir = BMO.direction(r)
            @test norm(cross(F2t - p, dir)) < 1e-9
            @test dot(F2t - p, dir) > 0
        end
    end

    @testset "Errors and promotion" begin
        @test_throws ArgumentError ConicMirror(0.1, 0.5, 0.2)
        @test_throws ArgumentError OffAxisConicMirror(0.1, 0.5, 0.09, 0.02)
        @test_throws ArgumentError EllipsoidalMirror(0.1, -0.1, 0.02)
        @test_throws ArgumentError EllipsoidalMirror(-0.04, 0.2, 0.02)
        @test_throws ArgumentError HyperbolicMirror(0.3, 0.15, 0.02)

        @test BMO.shape(OffAxisConicMirror(36, -1, 8, 8)).f isa Float64
        @test BMO.shape(ConicMirror(1, 0, 1)).f isa Float64
        @test BMO.shape(EllipsoidalMirror(2, 3, 1)).f isa Float64
        @test BMO.shape(OffAxisHyperbolicMirror(-1, 3, 1, 1)).f isa Float64
    end
end

end # module
