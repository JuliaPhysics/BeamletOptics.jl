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
    @testset "Testing calc_local_pos with GaussianBeamlet" begin
        # setup test beam
        offset = 50mm
        gb = GaussianBeamlet([0,0,0], [0,1,0])
        w, ~, ~, ~ = gauss_parameters(gb, offset)
        # setup untilted detector
        pd = Detector(5mm)
        translate_to3d!(pd, [0, offset, 0])
        system = System([pd])
        @testset "Untilted detector projection" begin
            solve_system!(system, gb)
            crop_factor = 2
            num_spots = 1000
            pts = BMO.calc_local_pos(pd; num_spots, crop_factor)
            # find local bounding box
            xmin = minimum(getindex.(pts, 1))
            xmax = maximum(getindex.(pts, 1))
            ymin = minimum(getindex.(pts, 2))
            ymax = maximum(getindex.(pts, 2))
            # test num spots and correct circle dims
            @test length(pts) == num_spots
            @test isapprox(w * crop_factor, -xmin, atol=1e-8)
            @test isapprox(w * crop_factor, xmax, atol=1e-8)
            @test isapprox(w * crop_factor, -ymin, atol=1e-8)
            @test isapprox(w * crop_factor, ymax, atol=1e-8)
        end
        @testset "Tilted detector projection" begin
            zrotate3d!(pd, deg2rad(45))
            solve_system!(system, gb)
            crop_factor = 1.5
            num_spots = 980
            pts = BMO.calc_local_pos(pd; num_spots, crop_factor)
            xmin = minimum(getindex.(pts, 1))
            xmax = maximum(getindex.(pts, 1))    
            ymin = minimum(getindex.(pts, 2))
            ymax = maximum(getindex.(pts, 2))    
            @test isapprox(w * crop_factor * sqrt(2), -xmin, atol=1e-8)
            @test isapprox(w * crop_factor * sqrt(2), xmax, atol=1e-8)
            @test isapprox(w * crop_factor, -ymin, atol=1e-8)
            @test isapprox(w * crop_factor, ymax, atol=1e-8)
        end
    end
end

end # MODULE