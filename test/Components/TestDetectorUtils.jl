module TestDetectorUtils

using BeamletOptics
using LinearAlgebra
using Test

const BMO = BeamletOptics

const mm = 1e-3

struct DummyHit <: BMO.AbstractDetectorHit end

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
            empty!(pd)
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

    @testset "Testing calc_local_pos with AstigmaticGaussianBeamlet" begin
        # setup test beam
        offset = 150mm
        agb = AstigmaticGaussianBeamlet([0,0,0], [0,1,0])
        w, ~, ~, ~ = gauss_parameters(agb, offset)
        # setup untilted detector
        pd = Detector(5mm)
        translate_to3d!(pd, [0, offset, 0])
        xrotate3d!(pd, deg2rad(45))
        system = System([pd])
        @testset "Untilted detector projection" begin
            solve_system!(system, agb)
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
            @test isapprox(w * crop_factor, -xmin, atol=5e-6)
            @test isapprox(w * crop_factor, xmax, atol=5e-6)
            @test isapprox(w * crop_factor * sqrt(2), -ymin, atol=5e-6)
            @test isapprox(w * crop_factor * sqrt(2), ymax, atol=5e-6)
        end
    end

    @testset "local_axes" begin
        pd = Detector(5mm)
        xrotate3d!(pd, deg2rad(30))
        zrotate3d!(pd, deg2rad(-50))
        x̂, ẑ = BMO.local_axes(pd)
        @test cross(x̂, ẑ) ≈ BMO.direction(pd)
        @test x̂ isa BMO.Point3 && ẑ isa BMO.Point3
    end

    @testset "calc_local_lims" begin
        pd = Detector(20mm)
        translate_to3d!(pd, [0, 50mm, 0])
        @testset "Ray hits" begin
            cs = CollimatedSource([0, 0, 0], [0, 1, 0], 5mm, 1e-6; num_rings = 2, num_rays = 41)
            solve_system!(System([pd]), cs)
            ray_hits = BMO.hits(pd)
            lims = BMO.calc_local_lims(pd)
            @test lims == BMO.calc_local_lims(pd, ray_hits, 1)
            @test lims == BMO.calc_local_lims(pd, ray_hits, nothing)
            x_min, x_max, z_min, z_max = BMO.calc_local_lims(pd, ray_hits, 2)
            @test x_max - x_min ≈ 2 * (lims[2] - lims[1])
            @test z_max - z_min ≈ 2 * (lims[4] - lims[3])
            @test (x_min + x_max) / 2 ≈ (lims[1] + lims[2]) / 2 atol = 1e-12
            @test (z_min + z_max) / 2 ≈ (lims[3] + lims[4]) / 2 atol = 1e-12
        end
        @testset "Beamlet hits" begin
            empty!(pd)
            solve_system!(System([pd]), GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0]))
            beamlet_hits = BMO.hits(pd)
            @test BMO.calc_local_lims(pd) == BMO.calc_local_lims(pd, beamlet_hits, 3)
            @test BMO.calc_local_lims(pd, beamlet_hits, nothing) ==
                  BMO.calc_local_lims(pd, beamlet_hits, 3)
        end
        @testset "Unknown hit type" begin
            err = try
                BMO.calc_local_lims(pd, [DummyHit()], nothing)
            catch e
                e
            end
            @test err isa ErrorException
            @test occursin("calc_local_lims not defined for", err.msg)
        end
    end
end

end # MODULE