module TestParabolicMirror

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics
const mm = 1e-3
const nm = 1e-9

@testset "Parabolic Mirrors" begin
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

        # Transformation of pierced mirror (testing AbstractCompositeSDF align3d! and reset_*)
        m_trans = ParabolicMirror(f, d; hole_diameter = hd)
        align3d!(m_trans, [1.0, 0.0, 0.0])
        # Ray along the aligned local y-axis ([1, 0, 0]) must miss
        beam_align = Beam(Ray([-50mm, 0, 0], [1.0, 0, 0]))
        solve_system!(StaticSystem([m_trans]), beam_align)
        @test length(BMO.rays(beam_align)) == 1

        # Ray just outside bore radius in aligned orientation must hit
        beam_align_hit = Beam(Ray([-50mm, 0, r_hit], [1.0, 0, 0]))
        solve_system!(StaticSystem([m_trans]), beam_align_hit)
        @test length(BMO.rays(beam_align_hit)) == 2

        # Reset rotation restores original alignment
        reset_rotation3d!(m_trans)
        beam_reset = Beam(Ray([0, -50mm, 0], [0, 1.0, 0]))
        solve_system!(StaticSystem([m_trans]), beam_reset)
        @test length(BMO.rays(beam_reset)) == 1

        # Invalid hole sizes are rejected
        @test_throws ArgumentError ParabolicMirror(f, d; hole_diameter = 0)
        @test_throws ArgumentError ParabolicMirror(f, d; hole_diameter = d)
        @test_throws ArgumentError ParabolicMirror(f, d; hole_diameter = 2d)
    end

    @testset "Off-Axis Parabolic Mirror with through-hole (Thorlabs OAP)" begin
        rfl = 100mm
        d = 50mm
        hd = 10mm

        # 1. Collimated through-hole
        oap_coll = OffAxisParabolicMirror(rfl, d; hole_diameter = hd, hole_axis = :collimated)
        @test BMO.shape(oap_coll) isa BMO.DifferenceSDF

        # Ray along local +y through aperture center passes straight through hole
        beam_coll = Beam(Ray([0, -50mm, 0], [0, 1.0, 0]))
        solve_system!(StaticSystem([oap_coll]), beam_coll)
        @test length(BMO.rays(beam_coll)) == 1

        # Ray outside bore hits the mirror
        beam_coll_hit = Beam(Ray([hd / 2 + 2mm, -50mm, 0], [0, 1.0, 0]))
        solve_system!(StaticSystem([oap_coll]), beam_coll_hit)
        @test length(BMO.rays(beam_coll_hit)) == 2

        # 2. Focused through-hole, including a non-right-angle segment.
        @testset "Focused bore ($(angle)°)" for angle in (90, 60)
            oap_foc = OffAxisParabolicMirror(rfl, d; angle,
                hole_diameter = hd, hole_axis = :focused)
            @test BMO.shape(oap_foc) isa BMO.DifferenceSDF

            solid = OffAxisParabolicMirror(rfl, d; angle)
            parent = BMO.shape(solid)
            # The segment origin is on the surface, not at the parent vertex.
            # Its parent vertex is (-x_off, sag_off, 0), so the focus is
            # (-x_off, sag_off - f, 0); at 90° this is (-rfl, 0, 0).
            sag_off = parent.x_off^2 / (4 * parent.f)
            focus = [-parent.x_off, sag_off - parent.f, 0.0]
            @test norm(focus) ≈ rfl
            f_dir = normalize(focus)

            # Rays along the bore pass in both directions. The solid mirror
            # must intercept them, so this cannot pass merely by missing it.
            for side in (-1, 1)
                ray = Ray(side * rfl * f_dir, -side * f_dir)
                @test !isnothing(BMO.intersect3d(solid, ray))
                beam_foc = Beam(ray)
                solve_system!(StaticSystem([oap_foc]), beam_foc)
                @test length(BMO.rays(beam_foc)) == 1
            end
        end

        # Error handling
        @test_throws ArgumentError OffAxisParabolicMirror(rfl, d; hole_diameter = 0)
        @test_throws ArgumentError OffAxisParabolicMirror(rfl, d; hole_diameter = d)
        @test_throws ArgumentError OffAxisParabolicMirror(rfl, d; hole_diameter = hd, hole_axis = :invalid)
    end

    @testset "Cassegrain and other mirrors with through-holes" begin
        d = 60mm
        hd = 15mm

        # ConicMirror with hole
        m_conic = ConicMirror(200mm, -1.0, d; hole_diameter = hd)
        @test BMO.shape(m_conic) isa BMO.DifferenceSDF
        b_conic = Beam(Ray([0, -50mm, 0], [0, 1.0, 0]))
        solve_system!(StaticSystem([m_conic]), b_conic)
        @test length(BMO.rays(b_conic)) == 1

        # HyperbolicMirror (Ritchey-Chrétien primary with hole)
        m_hyp = HyperbolicMirror(100mm, -200mm, d; hole_diameter = hd)
        @test BMO.shape(m_hyp) isa BMO.DifferenceSDF
        b_hyp = Beam(Ray([0, -50mm, 0], [0, 1.0, 0]))
        solve_system!(StaticSystem([m_hyp]), b_hyp)
        @test length(BMO.rays(b_hyp)) == 1

        # EllipsoidalMirror (Dall-Kirkham primary with hole)
        m_ell = EllipsoidalMirror(100mm, 200mm, d; hole_diameter = hd)
        @test BMO.shape(m_ell) isa BMO.DifferenceSDF
        b_ell = Beam(Ray([0, -50mm, 0], [0, 1.0, 0]))
        solve_system!(StaticSystem([m_ell]), b_ell)
        @test length(BMO.rays(b_ell)) == 1

        # SphericalMirror with hole
        m_sph = SphericalMirror(200mm, 10mm, d; hole_diameter = hd)
        @test BMO.shape(m_sph) isa BMO.DifferenceSDF
        b_sph = Beam(Ray([0, -50mm, 0], [0, 1.0, 0]))
        solve_system!(StaticSystem([m_sph]), b_sph)
        @test length(BMO.rays(b_sph)) == 1

        # RoundPlanoMirror with hole
        m_plano = RoundPlanoMirror(d, 10mm; hole_diameter = hd)
        @test BMO.shape(m_plano) isa BMO.DifferenceSDF
        b_plano = Beam(Ray([0, -50mm, 0], [0, 1.0, 0]))
        solve_system!(StaticSystem([m_plano]), b_plano)
        @test length(BMO.rays(b_plano)) == 1
    end

    @testset "Mirror constructors return Mirror" begin
        d = 25mm
        t = 5mm
        hd = 5mm

        m_plano = RoundPlanoMirror(d, t)
        @test m_plano isa Mirror
        @test BMO.shape(m_plano) isa BMO.PlanoSurfaceSDF

        m_plano_hole = RoundPlanoMirror(d, t; hole_diameter = hd)
        @test m_plano_hole isa Mirror

        m_sph = SphericalMirror(0.2, t, d)
        @test m_sph isa Mirror
        @test BMO.shape(m_sph) isa BMO.UnionSDF

        m_sph_hole = SphericalMirror(0.2, t, d; hole_diameter = hd)
        @test m_sph_hole isa Mirror

        m_rap = RightAnglePrismMirror(d, d)
        @test m_rap isa Mirror
        @test BMO.shape(m_rap) isa BMO.RightAnglePrismSDF
    end

    @testset "Parabola via ConicMirror" begin
        f = 100mm
        D = 60mm
        # ConicMirror(2f, -1, D) must build the same shape as ParabolicMirror(f, D), whose
        # focusing is traced in "On-Axis Parabolic Mirror" above
        sc = BMO.shape(ConicMirror(2f, -1, D))
        sp = BMO.shape(ParabolicMirror(f, D))
        for field in (:f, :k, :x_off, :diameter, :thickness)
            @test getfield(sc, field) ≈ getfield(sp, field)
        end
    end
end

end # module
