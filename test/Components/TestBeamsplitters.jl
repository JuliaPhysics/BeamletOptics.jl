module TestBeamsplitters

using BeamletOptics
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
        solve_system!(system, beam)

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

        # Retrace backside
        zrotate3d!(pbs, π)
        solve_system!(system, beam)

        @testset "Test children after retracing" begin
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

        @testset "Retrace after 45° CBS rotation" begin
            # Retrace
            zrotate3d!(cbs, π / 2)
            solve_system!(system, beam)

            # Test correct ray dirs
            p = BMO.rays(beam)
            t = BMO.rays(beam.children[1])
            r = BMO.rays(beam.children[2])

            @test BMO.direction(last(t)) == BMO.direction(first(p))
            @test BMO.direction(last(t)) == [0, 1, 0]
        end

        @testset "Retrace CBS backside" begin
            # Retrace backside
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
end

end