module TestHyperbolicMirror

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics
const mm = 1e-3

@testset "Hyperbolic Mirrors" begin
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

    @testset "Errors and promotion" begin
        @test_throws ArgumentError HyperbolicMirror(0.3, 0.15, 0.02)

        @test BMO.shape(OffAxisHyperbolicMirror(-1, 3, 1, 1)).f isa Float64
    end

    @testset "Off-axis Hyperbolic mirror with through-hole" begin
        b = 40mm
        a = 200mm
        x_off = 20mm
        D = 20mm
        hd = 4mm

        m = OffAxisHyperbolicMirror(-b, a, x_off, D; hole_diameter = hd)
        m0 = OffAxisHyperbolicMirror(-b, a, x_off, D)
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

        @test_throws ArgumentError OffAxisHyperbolicMirror(-b, a, x_off, D; hole_diameter = 0)
        @test_throws ArgumentError OffAxisHyperbolicMirror(-b, a, x_off, D; hole_diameter = D)
    end
end

end # module
