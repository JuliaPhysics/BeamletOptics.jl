module TestParabolicMirror

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

    @testset "Integer Input Promotion" begin
        oap_int = OffAxisParabolicMirror(1, 1; angle = 90)
        sdf_int = BMO.shape(oap_int)
        @test sdf_int.f isa Float64
        @test sdf_int.f ≈ 0.5
        @test sdf_int.x_off ≈ 1.0
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

    @testset "On-Axis Parabolic Mirror" begin
        f = 100mm
        d = 60mm
        m = ParabolicMirror(f, d)
        sdf = BMO.shape(m)

        @test sdf.f ≈ f
        @test sdf.x_off == 0
        @test sdf.thickness ≈ 2.25mm + 10mm
        @test BMO.shape(ParabolicMirror(1, 1)).f isa Float64

        # Rays parallel to the optical axis must pass through the focus after reflection
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

        # Rays outside of the aperture miss the mirror
        beam = Beam(Ray([35mm, -50mm, 0], [0, 1, 0]))
        solve_system!(StaticSystem([m]), beam)
        @test length(BMO.rays(beam)) == 1
    end

    @testset "On-Axis Parabolic Mirror with bore (Cassegrain primary)" begin
        f = 100mm
        d = 60mm
        hd = 10mm
        m = ParabolicMirror(f, d; hole_diameter = hd)
        sdf = BMO.shape(m)

        @test sdf isa BMO.DifferenceSDF

        # A ray down the bore axis misses (passes straight through the hole).
        beam = Beam(Ray([0, -50mm, 0], [0, 1, 0]))
        solve_system!(StaticSystem([m]), beam)
        @test length(BMO.rays(beam)) == 1

        # A ray just outside the bore radius hits, with t matching the analytic sag.
        r_hit = hd / 2 + 1mm
        beam2 = Beam(Ray([r_hit, -50mm, 0], [0, 1, 0]))
        solve_system!(StaticSystem([m]), beam2)
        @test length(BMO.rays(beam2)) == 2
        r = BMO.rays(beam2)[1]
        sag = r_hit^2 / (4f)
        @test length(r) ≈ 50mm - sag atol=1e-9

        # Invalid hole sizes are rejected
        @test_throws ArgumentError ParabolicMirror(f, d; hole_diameter = 0)
        @test_throws ArgumentError ParabolicMirror(f, d; hole_diameter = d)
        @test_throws ArgumentError ParabolicMirror(f, d; hole_diameter = 2d)
    end

    @testset "Bounding Sphere" begin
        f = 100mm
        r_max = 30mm
        thickness = 20mm
        sag_max = r_max^2 / (4f)

        sdf = BMO.OffAxisParaboloidSDF(f, 0.0, 2r_max, thickness)
        center, r = BMO.bounding_sphere(sdf)

        @test center ≈ BMO.Point3(0, (thickness - sag_max) / 2, 0)
        @test r ≈ sqrt(r_max^2 + ((thickness + sag_max) / 2)^2) + 0.05

        # bounding_box must transform the sphere into a symmetric, aperture-covering box
        xmin, xmax, ymin, ymax, zmin, zmax = BMO.bounding_box(sdf)
        @test xmin ≈ -xmax
        @test zmin ≈ -zmax
        @test xmax ≈ zmax # rotationally symmetric about the y-axis since x_off = 0
        @test xmax >= r_max
    end
end

end # module
