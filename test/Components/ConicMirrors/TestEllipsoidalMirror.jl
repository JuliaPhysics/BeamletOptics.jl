module TestEllipsoidalMirror

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics
const mm = 1e-3

@testset "Ellipsoidal Mirrors" begin
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

    @testset "Errors and promotion" begin
        @test_throws ArgumentError EllipsoidalMirror(0.1, -0.1, 0.02)
        @test_throws ArgumentError EllipsoidalMirror(-0.04, 0.2, 0.02)

        @test BMO.shape(EllipsoidalMirror(2, 3, 1)).f isa Float64
    end

    @testset "Off-axis Ellipsoidal mirror with through-hole" begin
        s1 = 300mm
        s2 = 150mm
        x_off = 60mm
        D = 30mm
        hd = 6mm

        m = OffAxisEllipsoidalMirror(s1, s2, x_off, D; hole_diameter = hd)
        m0 = OffAxisEllipsoidalMirror(s1, s2, x_off, D)
        @test BMO.shape(m) isa BMO.DifferenceSDF

        beam = Beam(Ray([0, -50mm, 0], [0, 1.0, 0]))
        solve_system!(StaticSystem([m]), beam)
        @test length(BMO.rays(beam)) == 1

        beam0 = Beam(Ray([0, -50mm, 0], [0, 1.0, 0]))
        solve_system!(StaticSystem([m0]), beam0)
        @test length(BMO.rays(beam0)) == 2

        beam_off = Beam(Ray([hd / 2 + 2mm, -50mm, 0], [0, 1.0, 0]))
        solve_system!(StaticSystem([m]), beam_off)
        @test length(BMO.rays(beam_off)) == 2

        @test_throws ArgumentError OffAxisEllipsoidalMirror(s1, s2, x_off, D; hole_diameter = 0)
        @test_throws ArgumentError OffAxisEllipsoidalMirror(s1, s2, x_off, D; hole_diameter = D)
    end
end

end # module
