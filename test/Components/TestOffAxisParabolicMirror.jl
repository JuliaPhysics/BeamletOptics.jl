module TestOffAxisParabolicMirror

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics
const mm = 1e-3
const nm = 1e-9

@testset "Off-Axis Parabolic Mirror" begin
    @testset "Constructor Geometry (90°)" begin
        rfl = 200mm
        d = 60mm
        oap = OffAxisParabolicMirror(rfl, d; angle = 90)
        sdf = BMO.shape(oap)

        @test sdf.f ≈ 100mm
        @test sdf.x_off ≈ 200mm
        @test sdf.diameter == 60mm
        @test sdf.thickness >= 37.75mm # max sag is 27.75mm + 10mm margin
    end

    @testset "Constructor Geometry (60°)" begin
        rfl = 200mm
        d = 60mm
        oap = OffAxisParabolicMirror(rfl, d; angle = 60)
        sdf = BMO.shape(oap)

        @test sdf.f ≈ 200mm * (cos(deg2rad(30))^2)
        @test sdf.x_off ≈ 200mm * sin(deg2rad(60))
    end

    @testset "Raytracing Focus Precision (90°)" begin
        rfl = 200mm
        d = 60mm
        λ = 1550nm
        angle = 90

        oap = OffAxisParabolicMirror(rfl, d; angle = angle)
        sdf = BMO.shape(oap)
        f_parent = sdf.f
        x_off = sdf.x_off

        translate_to3d!(oap, [x_off, f_parent, 0])
        beam = CollimatedSource([x_off, -50mm, 0], [0, 1, 0], 25mm, λ; num_rings = 3, num_rays = 80)

        pd = Detector(60mm, true)
        zrotate3d!(pd, deg2rad(angle))
        translate_to3d!(pd, [0, f_parent, 0])

        sys = StaticSystem([oap, pd])
        solve_system!(sys, beam)

        @test length(pd.hits) == 80
        spots = spot_diagram(pd)
        max_spot_radius = maximum(norm.(spots))
        @test max_spot_radius < 1e-9 # Sub-nanometer geometric point focus
    end

    @testset "Raytracing Focus Precision (60°)" begin
        rfl = 200mm
        d = 60mm
        λ = 1550nm
        angle = 60

        oap = OffAxisParabolicMirror(rfl, d; angle = angle)
        sdf = BMO.shape(oap)
        f_parent = sdf.f
        x_off = sdf.x_off

        translate_to3d!(oap, [x_off, f_parent, 0])
        beam = CollimatedSource([x_off, -50mm, 0], [0, 1, 0], 25mm, λ; num_rings = 3, num_rays = 80)

        y_focus = f_parent - rfl * cos(deg2rad(angle))
        pd = Detector(60mm, true)
        zrotate3d!(pd, deg2rad(angle))
        translate_to3d!(pd, [0, y_focus, 0])

        sys = StaticSystem([oap, pd])
        solve_system!(sys, beam)

        @test length(pd.hits) == 80
        spots = spot_diagram(pd)
        max_spot_radius = maximum(norm.(spots))
        @test max_spot_radius < 1e-9
    end
end

end # module
