using Test
using BeamletOptics
using StaticArrays
using LinearAlgebra
using GeometryBasics

@testset "Coatings" begin
    @testset "Uncoated (Fresnel)" begin
        c = Uncoated()
        n1 = 1.0
        n2 = 1.5
        λ = 1000e-9
        rs, rp, ts, tp = BeamletOptics.coating_coefficients(c, n1, n2, λ, 0.0)
        expected_r = (n1 - n2) / (n1 + n2)
        @test isapprox(real(rs), expected_r, atol = 1e-6)
    end

    @testset "SimpleCoating" begin
        c = SimpleCoating(0.5, 0.5)
        rs, rp, ts, tp = BeamletOptics.coating_coefficients(c, 1.0, 1.0, 500e-9, 0.0)
        @test isapprox(abs2(rs), 0.5, atol = 1e-6)
    end

    @testset "Spatial Coatings" begin
        # 1. Test SpatiallyVaryingCoating logic directly
        # Define a coating that is Mirror for x>0, AR for x<=0
        c_spatial = BeamletOptics.SpatiallyVaryingCoating((p, n) -> begin
            if p[1] > 0
                return SimpleCoating(1.0) # Mirror
            else
                return SimpleCoating(1.0, 0.0) # AR (T=1, R=0)
            end
        end)

        # Point at x=1
        # Mock optic with identity transform (Global = Local)
        # We need a struct with a shape field that is an SDF.
        # BoxSDF at origin works.
        mock_shape = BeamletOptics.BoxSDF(2.0, 2.0, 2.0) # centered at 0
        mock_optic = (; shape = mock_shape)

        c1 = BeamletOptics.get_coating_at(
            c_spatial, mock_optic, Point3(1, 0, 0), Point3(0, 0, 1))
        rs, _, _, _ = BeamletOptics.coating_coefficients(c1, 1.0, 1.5, 500e-9, 0.0)
        @test isapprox(abs2(rs), 1.0, atol = 1e-6)

        # Point at x=-1
        c2 = BeamletOptics.get_coating_at(
            c_spatial, mock_optic, Point3(-1, 0, 0), Point3(0, 0, 1))
        rs, _, _, _ = BeamletOptics.coating_coefficients(c2, 1.0, 1.5, 500e-9, 0.0)
        @test isapprox(abs2(rs), 0.0, atol = 1e-6)
    end

    @testset "Lens Constructor (Spatial)" begin
        # Just verify the constructor runs and produces a lens with SpatiallyVaryingCoating
        front = SphericalSurface(0.05, 0.02)
        back = SphericalSurface(-0.05, 0.02)

        l = Lens(front, back, 0.01, x -> 1.5, SimpleCoating(0.0, 1.0), SimpleCoating(1.0))

        @test BeamletOptics.coating(l) isa BeamletOptics.SpatiallyVaryingCoating
    end

    @testset "Meniscus Lens Spatial Coating" begin
        # Create a meniscus lens (Positive Meniscus)
        # Front convex (R1 > 0), Back convex (R2 > 0)
        r1 = 0.05
        r2 = 0.10
        thick = 0.005

        front = SphericalSurface(r1, 0.02)
        back = SphericalSurface(r2, 0.02)

        c_front = SimpleCoating(0.1)
        c_back = SimpleCoating(0.9)

        # Use our new constructor
        l = Lens(front, back, thick, x -> 1.5, c_front, c_back)

        # Check coating at Front Apex area (approx Y=0 or shifted depending on exact construction SDF logic)
        # get_coating_at will check distances to the surfaces.
        # (0,0,0) should be much closer to front than back (min dist).
        p_front = Point3(0.0, 0.0, 0.0)
        c_f = BeamletOptics.get_coating_at(
            BeamletOptics.coating(l), l, p_front, Point3(0, 0, 1))
        # Note: We rely on the generic selector finding the min distance.
        @test c_f.reflection(500e-9, 0.0) == 0.1

        # Check near back (Y = thickness)
        p_back = Point3(0.0, thick, 0.0)
        c_b = BeamletOptics.get_coating_at(
            BeamletOptics.coating(l), l, p_back, Point3(0, 0, 1))
        @test c_b.reflection(500e-9, 0.0) == 0.9
    end

    @testset "Side Coating" begin
        front = SphericalSurface(0.05, 0.02) # D=0.02 -> R_cyl = 0.01
        back = SphericalSurface(-0.05, 0.02)
        thick = 0.01

        c_front = SimpleCoating(0.1)
        c_back = SimpleCoating(0.1)
        c_side = SimpleCoating(0.5)

        l = Lens(front, back, thick, x -> 1.5, c_front, c_back, c_side)

        # Point on the cylinder.
        # Radius is 0.01.
        p_side = Point3(0.01, thick / 2, 0.0) # On the surface

        c_s = BeamletOptics.get_coating_at(
            BeamletOptics.coating(l), l, p_side, Point3(1, 0, 0))
        @test c_s.reflection(500e-9, 0.0) == 0.5
    end

    @testset "Bragg Mirror (TMM)" begin
        # Design a Bragg mirror for 1000nm
        λ0 = 1000e-9
        nH = 2.3 # TiO2 approx
        nL = 1.46 # SiO2 approx
        dH = λ0 / (4 * nH)
        dL = λ0 / (4 * nL)

        # 10 pairs
        layers = ThinFilmLayer[]
        for i in 1:10
            push!(layers, ThinFilmLayer(nH, dH))
            push!(layers, ThinFilmLayer(nL, dL))
        end

        bragg = MultilayerCoating(layers, λ0)

        # At design wavelength, R should be very high
        rs, _, _, _ = BeamletOptics.coating_coefficients(bragg, 1.0, 1.5, λ0, 0.0)
        R = abs2(rs)
        @test R > 0.99

        # At far wavelength, R should drop
        rs_far, _, _, _ = BeamletOptics.coating_coefficients(bragg, 1.0, 1.5, 1.5 * λ0, 0.0)
        @test abs2(rs_far) < 0.5 # Should be significantly lower
    end

    @testset "Prism Coating" begin
        # Verify Prism accepts coating and stores it
        p = RightAnglePrism(0.01, 0.01, x -> 1.5, SimpleCoating(0.1))
        @test BeamletOptics.coating(p) isa SimpleCoating
        # SimpleCoating stores functions for reflection/transmission.
        # We need to call it to check the value.
        @test BeamletOptics.coating(p).reflection(500e-9, 0.0) == 0.1
    end
end
