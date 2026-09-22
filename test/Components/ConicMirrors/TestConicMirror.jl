module TestConicMirror

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics
const mm = 1e-3

@testset "Conic Mirrors" begin
    @testset "Kinematics" begin
        # Prolate ellipsoid with conjugate foci at (0, -s1, 0) and (0, -s2, 0), built from
        # the general (R, k) parameterization so the test covers every conic mirror type.
        s1 = 300mm
        s2 = 150mm
        D = 40mm
        R = 2 * s1 * s2 / (s1 + s2)
        k = -((s2 - s1) / (s2 + s1))^2
        m = ConicMirror(R, k, D)
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

        @test BMO.shape(OffAxisConicMirror(36, -1, 8, 8)).f isa Float64
        @test BMO.shape(ConicMirror(1, 0, 1)).f isa Float64
    end

    @testset "Off-axis Conic mirror with through-hole" begin
        # Concave, valid off-axis mirror (same parameterization as the "Kinematics" testset
        # above, offset like the off-axis testset in TestEllipsoidalMirror.jl).
        s1 = 300mm
        s2 = 150mm
        D = 30mm
        x_off = 60mm
        R = 2 * s1 * s2 / (s1 + s2)
        k = -((s2 - s1) / (s2 + s1))^2
        hd = 6mm

        m = OffAxisConicMirror(R, k, x_off, D; hole_diameter = hd)
        m0 = OffAxisConicMirror(R, k, x_off, D)
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

        @test_throws ArgumentError OffAxisConicMirror(R, k, x_off, D; hole_diameter = 0)
        @test_throws ArgumentError OffAxisConicMirror(R, k, x_off, D; hole_diameter = D)
    end
end

end # module
