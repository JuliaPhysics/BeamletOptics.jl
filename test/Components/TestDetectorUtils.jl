module TestDetectorUtils

using BeamletOptics
using Test

const BMO = BeamletOptics

const mm = 1e-3

@testset "DetectorUtils" begin
    @testset "Testing calc_center_point with Centroid" begin
        @testset "Uniform weights" begin
            xs = [1.0, 2.0, 3.0]
            zs = [4.0, 5.0, 6.0]
            projection_factor = [1.0, 1.0, 1.0]
            x0, z0 = BMO.calc_center_point(BMO.Centroid(), xs, zs, projection_factor)
            @test x0 ≈ 2.0
            @test z0 ≈ 5.0
        end

        @testset "Biased weights" begin
            # biased weights
            xs = [1.0, 2.0, 3.0]
            zs = [4.0, 5.0, 6.0]
            projection_factor = [10.0, 1.0, 1.0] 
            x0, z0 = BMO.calc_center_point(BMO.Centroid(), xs, zs, projection_factor)
            @test x0 ≈ 1.25
            @test z0 ≈ 4.25
        end

        @testset "iszero fallback" begin
            xs = [1.0, 10.0, 100.0]
            zs = [2.0, 20.0, 200.0]
            projection_factor = [0.0, 0.0, 0.0]
            x0, z0 = BMO.calc_center_point(BMO.Centroid(), xs, zs, projection_factor)
            @test x0 ≈ 50.5
            @test z0 ≈ 101.0
        end
    end
end

end # MODULE






