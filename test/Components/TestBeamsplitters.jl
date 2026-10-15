module TestBeamsplitters

using BeamletOptics
using LinearAlgebra
using Test

const BMO = BeamletOptics

const mm = 1e-3

@testset "Beamsplitters" begin
    N0 = 1.5
    @testset "Testing RectangularPlateBeamsplitter with Beam" begin
        # Init splitter
        N0 = 1.5
        pbs = RectangularPlateBeamsplitter(36mm, 25mm, 1mm, n -> N0)
        system = System([pbs])
        beam = Beam([0, -50mm, 0], [0, 1, 0], 1e-6)
        # Trace normally
        zrotate3d!(pbs, deg2rad(45))
        # Keep geometry regressions from growing an unbounded beam tree in CI.
        solve_system!(system, beam; depth_max=4)

        @testset "Test pos/dir" begin
            @test position(pbs) == zeros(3)
            @test orientation(pbs) ≈ orientation(pbs.substrate)
        end

        @testset "Test children after tracing" begin
            p = beam.rays
            t = beam.children[1].rays
            r = beam.children[2].rays
            # no of rays
            @test length(p) == 1
            @test length(t) == 2
            @test length(r) == 1
            # correct ref. index
            @test all(BMO.refractive_index.(p) .== 1)
            @test all(BMO.refractive_index.(t) .== [N0, 1])
            @test all(BMO.refractive_index.(r) .== 1)
            # correct dir
            @test BMO.direction(first(p)) ≈ BMO.direction(last(t))
            @test BMO.direction(first(r)) ≈ [1, 0, 0]
        end

        # Solve again with the backside
        zrotate3d!(pbs, π)
        solve_system!(system, beam; depth_max=4)

        @testset "Test children after solving again" begin
            p = beam.rays
            t = beam.children[1].rays
            r = beam.children[2].rays
            # no of rays
            @test length(p) == 2
            @test length(t) == 1
            @test length(r) == 2
            # correct ref. index
            @test all(BMO.refractive_index.(p) .== [1, N0])
            @test all(BMO.refractive_index.(t) .== 1)
            @test all(BMO.refractive_index.(r) .== [N0, 1])
            # correct dir
            @test BMO.direction(first(p)) ≈ BMO.direction(last(t))
            @test BMO.direction(last(r)) ≈ [1, 0, 0]
        end
    end

    @testset "Testing CubeBeamsplitter with Beam" begin
        # Init splitter
        cbs = CubeBeamsplitter(25e-3, n -> N0)
        translate3d!(cbs, [0, 50mm, 0])
        system = System([cbs])
        beam = Beam([0, 0, 0], [0, 1, 0], 1e-6)

        @testset "Initial CBS tracing" begin
            # Trace normally
            solve_system!(system, beam)
            # Test correct ray length, ref. indices, dirs
            p = BMO.rays(beam)
            t = BMO.rays(beam.children[1])
            r = BMO.rays(beam.children[2])

            @test length(p) == 2
            @test length(r) == 2
            @test length(t) == 2
            @test BMO.refractive_index.(p) == [1, N0]
            @test BMO.refractive_index.(t) == [N0, 1]
            @test BMO.refractive_index.(r) == [N0, 1]
            @test BMO.direction(last(t)) ≈ BMO.direction(first(p))
            @test BMO.direction(last(r)) ≈ [-1, 0, 0]
        end

        @testset "Solve again after 45° CBS rotation" begin
            # Solve again
            zrotate3d!(cbs, π / 2)
            solve_system!(system, beam)

            # Test correct ray dirs
            p = BMO.rays(beam)
            t = BMO.rays(beam.children[1])
            r = BMO.rays(beam.children[2])

            @test BMO.direction(last(t)) ≈ BMO.direction(first(p))
            @test BMO.direction(last(t)) ≈ [0, 1, 0]
            @test BMO.position(p[2]) ≈ [0, 50mm - 12.5mm, 0]
        end

        @testset "Solve again, CBS backside" begin
            # Solve again, backside
            zrotate3d!(cbs, π / 2)
            solve_system!(system, beam)

            # Test correct ray length, ref. indices, dirs
            p = BMO.rays(beam)
            t = BMO.rays(beam.children[1])
            r = BMO.rays(beam.children[2])

            @test length(p) == 2
            @test length(r) == 2
            @test length(t) == 2
            @test BMO.refractive_index.(p) == [1, N0]
            @test BMO.refractive_index.(t) == [N0, 1]
            @test BMO.refractive_index.(r) == [N0, 1]
            @test BMO.direction(last(t)) ≈ BMO.direction(first(p))
            @test BMO.direction(last(r)) ≈ [-1, 0, 0]
        end
    end

    @testset "CubeBeamsplitter and RightAnglePrism with a constant refractive index" begin
        prism = RightAnglePrism(25e-3, 20e-3, N0)
        @test prism isa Prism
        @test BMO.refractive_index(prism, 1e-6) == N0
        @test BMO.refractive_index(prism, 500e-9) == N0

        cbs = CubeBeamsplitter(25e-3, N0; reflectance = 0.3)
        @test BMO.refractive_index(cbs, 1e-6) == N0
        @test BMO.refractive_index(cbs, 500e-9) == N0
        @test cbs.coating.reflectance ≈ CubeBeamsplitter(25e-3, n -> N0; reflectance = 0.3).coating.reflectance
        translate3d!(cbs, [0, 50mm, 0])
        beam = Beam([0, 0, 0], [0, 1, 0], 1e-6)
        solve_system!(System([cbs]), beam)
        @test BMO.refractive_index.(BMO.rays(beam)) == [1, N0]
        @test length(beam.children) == 2
    end

    @testset "depth_max branch limiting" begin
        beamsplitter = CubeBeamsplitter(25e-3, n -> N0)
        translate3d!(beamsplitter, [0, 50mm, 0])
        system = System([beamsplitter])

        beam = Beam([0, 0, 0], [0, 1, 0], 1e-6)
        solve_system!(system, beam; depth_max=0)
        @test isempty(beam.children)

        beam = Beam([0, 0, 0], [0, 1, 0], 1e-6)
        solve_system!(system, beam; depth_max=1)
        @test length(beam.children) == 2
    end

    @testset "depth_max default" begin
        beamsplitter = CubeBeamsplitter(25e-3, n -> N0)
        translate3d!(beamsplitter, [0, 50mm, 0])
        system = System([beamsplitter])
        beam = Beam([0, 0, 0], [0, 1, 0], 1e-6)

        @test BMO.get_default_depth_max() == 100
        solve_system!(system, beam)
        @test length(beam.children) == 2
    end

    @testset "Field of a polarized ray behind the plate" begin
        # the coating transmits half of the power into the glass, the uncoated back takes its
        # Fresnel loss: the field of the transmitted ray is transverse and has this amplitude
        n = 1.5
        θ_i = π / 4
        θ_t = asin(sin(θ_i) / n)
        r_s = (cos(θ_i) - n * cos(θ_t)) / (cos(θ_i) + n * cos(θ_t))
        r_p = (n * cos(θ_i) - cos(θ_t)) / (n * cos(θ_i) + cos(θ_t))
        for (e, R_back) in (([0.0, 0, 1], r_s^2), ([1.0, 0, 0], r_p^2))
            pbs = RectangularPlateBeamsplitter(25e-3, 25e-3, 5e-3, λ -> n)
            zrotate3d!(pbs, π / 4)
            beam = Beam(PolarizedRay([0.0, -0.1, 0], [0.0, 1, 0], 1e-6, e))
            solve_system!(System([pbs]), beam)
            transmitted, reflected = BMO.children(beam)
            inside, outside = BMO.rays(transmitted)
            @test BMO.refractive_index(inside) == n
            for ray in (inside, outside, first(BMO.rays(reflected)))
                @test abs(dot(BMO.direction(ray), BMO.polarization(ray))) < 1e-12
            end
            @test sum(abs2, BMO.polarization(first(BMO.rays(reflected)))) ≈ 1 / 2
            @test sum(abs2, BMO.polarization(outside)) ≈ (1 - R_back) / 2
            # in the glass the power is n cos(θ) |E|² per area of the surface
            @test n * cos(θ_t) * sum(abs2, BMO.polarization(inside)) ≈ cos(θ_i) / 2
        end
    end
end

end
