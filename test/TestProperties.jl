module TestProperties

using Test
using BeamletOptics

const BMO = BeamletOptics

"""Returns the value of the property `name` in the list `props`, or `nothing`."""
prop(props, name) = (i = findfirst(p -> p.first == name, props); isnothing(i) ? nothing : props[i].second)

@testset "properties" begin
    @testset "Conventions for all component types" begin
        objs = Any[
            SphericalLens(50e-3, -50e-3, 5e-3, 25.4e-3),
            SphericalLens(50e-3, 50e-3, 0, 25.4e-3),
            Lens(EvenAsphericalSurface(20e-3, 25e-3, -1.0, [0, 1e-5]), 8e-3, λ -> 1.5),
            Lens(CylindricalSurface(5.2e-3, 10e-3, 20e-3), 5.9e-3, λ -> 1.517),
            SphericalDoubletLens(50e-3, -40e-3, -100e-3, 4e-3, 2e-3, 25.4e-3, λ -> 1.5, λ -> 1.6),
            SquarePlanoMirror2D(25e-3), RectangularPlanoMirror(25e-3, 20e-3, 5e-3),
            RoundPlanoMirror(25e-3, 5e-3), SphericalMirror(100e-3, 25e-3, 5e-3),
            ParabolicMirror(25e-3, 50e-3), OffAxisParabolicMirror(50e-3, 25e-3; angle = 90),
            RightAnglePrismMirror(10e-3, 10e-3), Retroreflector(1.0),
            Detector(10e-3), ThinBeamsplitter(10e-3, 10e-3), RoundThinBeamsplitter(10e-3),
            RoundPlateBeamsplitter(25e-3, 3e-3, λ -> 1.5),
            RectangularPlateBeamsplitter(25e-3, 25e-3, 3e-3, λ -> 1.5),
            CubeBeamsplitter(10e-3, λ -> 1.5), RectangularCompensatorPlate(25e-3, 25e-3, 3e-3, λ -> 1.5),
            PolarizationFilter(10e-3), RoundPolarizationFilter(10e-3),
            RoundLinearPolarizer(25e-3, 1e-3, 1e-3, λ -> 1.5), RightAnglePrism(10e-3, 10e-3, λ -> 1.5),
            ObjectGroup([RoundPlanoMirror(25e-3, 5e-3), Detector(10e-3)]),
            Beam([0, 0, 0.0], [0, 1.0, 0], 1e-6), GaussianBeamlet([0, 0, 0.0], [0, 1.0, 0], 1e-6, 1e-3),
            CollimatedSource([0, 0, 0.0], [0, 1.0, 0], 5e-3; num_rings = 2),
            PointSource([0, 0, 0.0], [0, 1.0, 0], 0.1; num_rings = 2),
            CollimatedGaussianBeamletSource([0, 0, 0.0], [0, 1.0, 0], 5e-3, 1e-6, 1e-3; n_grid = 4)]
        for obj in objs
            props = BMO.properties(obj)
            @test props isa Vector{Pair{String, Any}}
            @test first(props) == ("Type" => string(nameof(typeof(obj))))
            # unique names, a unit in brackets at the end, if any
            names = first.(props)
            @test allunique(names)
            @test all(n -> !occursin('[', n) || occursin(r" \[[^\]]+\]$", n), names)
            @test prop(props, "Position [m]") ≈ collect(position(obj))
        end
    end

    @testset "Pose by the kinematic trait" begin
        lens = SphericalLens(50e-3, -50e-3, 5e-3, 25.4e-3)
        translate3d!(lens, [1e-3, 2e-3, 3e-3])
        zrotate3d!(lens, π / 2)
        props = BMO.properties(lens)
        @test prop(props, "Position [m]") ≈ [1e-3, 2e-3, 3e-3]
        @test prop(props, "Optical axis") ≈ [-1, 0, 0] atol = 1e-12
        beam = Beam([0, 0, 0.0], [1.0, 0, 0], 1e-6)
        @test prop(BMO.properties(beam), "Direction") ≈ [1, 0, 0]
        @test isnothing(prop(BMO.properties(beam), "Optical axis"))
        # a shape: its type and dimensions, its pose by default_properties
        cyl = BMO.CylinderSDF(5e-3, 1e-3)
        @test BMO.properties(cyl) == ["Shape" => "CylinderSDF"]
        @test first.(BMO.default_properties(cyl)) == ["Type", "Position [m]", "Optical axis"]
        # no root type of BMO
        @test_throws MethodError BMO.properties(42)
        @test_throws MethodError BMO.default_properties(42)
    end

    @testset "Component values" begin
        lens = SphericalLens(50e-3, -50e-3, 5e-3, 25.4e-3, λ -> 1.5)
        props = BMO.properties(lens)
        @test prop(props, "Shape") == "UnionSDF"
        @test prop(props, "Thickness [m]") ≈ 5e-3
        @test prop(props, "n(λ₀)") ≈ 1.5
        @test prop(props, "λ₀ [m]") == get_default_wavelength()
        @test prop(props, "Index model") == "function of λ"
        glass = DiscreteRefractiveIndex([1e-6, 2e-6], [1.5, 1.4])
        @test prop(BMO.properties(SphericalLens(50e-3, -50e-3, 5e-3, 25.4e-3, glass)), "Index model") ==
              "DiscreteRefractiveIndex"
        # Only defined at its wavelengths, which need not include λ₀
        glass = DiscreteRefractiveIndex([633e-9, 532e-9], [1.515, 1.519])
        props = BMO.properties(SphericalLens(50e-3, -50e-3, 5e-3, 25.4e-3, glass))
        @test prop(props, "Wavelengths [m]") == [532e-9, 633e-9]
        @test prop(props, "Indices") == [1.519, 1.515]
        @test !any(p -> first(p) == "n(λ₀)", props)

        mirror = RoundPlanoMirror(25e-3, 5e-3)
        @test prop(BMO.properties(mirror), "Diameter [m]") ≈ 25e-3
        @test prop(BMO.properties(mirror), "Thickness [m]") ≈ 5e-3
        @test prop(BMO.properties(SphericalMirror(100e-3, 25e-3, 5e-3)), "Shape") == "UnionSDF"
        @test prop(BMO.properties(ParabolicMirror(25e-3, 50e-3)), "Conic constant") == -1
        pierced = ConicMirror(200e-3, -0.5, 25e-3; thickness = 4e-3, hole_diameter = 5e-3)
        @test prop(BMO.properties(pierced), "Thickness [m]") == 4e-3

        bs = ThinBeamsplitter(10e-3, 20e-3; reflectance = 0.3)
        @test prop(BMO.properties(bs), "Reflectance") ≈ 0.3
        @test prop(BMO.properties(bs), "Transmittance") ≈ 0.7
        @test prop(BMO.properties(bs), "Size [m]") ≈ [10e-3, 20e-3]
        pbs = RoundPlateBeamsplitter(25e-3, 3e-3, λ -> 1.5; reflectance = 0.4)
        @test prop(BMO.properties(pbs), "Reflectance") ≈ 0.4
        @test prop(BMO.properties(pbs), "Thickness [m]") ≈ 3e-3
        @test prop(BMO.properties(pbs), "Parts") == 2
        @test prop(BMO.properties(CubeBeamsplitter(10e-3, λ -> 1.5; reflectance = 0.2)), "Reflectance") ≈ 0.2

        pf = PolarizationFilter(10e-3)
        @test prop(BMO.properties(pf), "Transmission axis") ≈ collect(transmission_axis(pf))
        lp = RoundLinearPolarizer(25e-3, 1e-3, 1e-3, λ -> 1.5)
        @test prop(BMO.properties(lp), "Transmission axis") ≈ collect(transmission_axis(lp))

        group = ObjectGroup([RoundPlanoMirror(25e-3, 5e-3), Detector(10e-3), Detector(5e-3)])
        @test prop(BMO.properties(group), "Parts") == 3
    end

    @testset "Detector" begin
        pd = Detector(10e-3)
        zrotate3d!(pd, π / 3)
        props = BMO.properties(pd)
        @test prop(props, "Size [m]") ≈ [10e-3, 10e-3]
        @test prop(props, "Hits") == 0
        @test prop(props, "Stops beams") == true
        translate3d!(pd, [0, 0.2, 0])
        system = System([pd])
        beam = Beam([0, 0, 0.0], [0, 1.0, 0], 1e-6)
        solve_system!(system, beam)
        @test prop(BMO.properties(pd), "Hits") == 1
    end

    @testset "Sources" begin
        gb = GaussianBeamlet([0, 0, 0.0], [0, 1.0, 0], 1e-6, 1e-3)
        props = BMO.properties(gb)
        @test prop(props, "Wavelength [m]") == 1e-6
        @test prop(props, "Waist radius [m]") == 1e-3
        @test prop(props, "Rayleigh range [m]") ≈ rayleigh_range(gb)
        cs = CollimatedSource([0, 0, 0.0], [0, 1.0, 0], 5e-3, 633e-9; num_rings = 2)
        @test prop(BMO.properties(cs), "Diameter [m]") == 5e-3
        @test prop(BMO.properties(cs), "Beams") == length(BMO.beams(cs))
        @test prop(BMO.properties(cs), "Wavelength [m]") == 633e-9
        @test prop(BMO.properties(cs), "Optical axis") ≈ [0, 1, 0]
        ps = PointSource([0, 0, 0.0], [0, 1.0, 0], 0.1; num_rings = 2)
        @test prop(BMO.properties(ps), "Numerical aperture") ≈ sin(0.1)
        ag = CollimatedGaussianBeamletSource([0, 0, 0.0], [0, 1.0, 0], 5e-3, 1e-6, 1e-3; n_grid = 4)
        @test prop(BMO.properties(ag), "Wavelength [m]") ≈ 1e-6
    end

    @testset "Extension by a custom component" begin
        # defined in a module like a user package
        m = Module()
        Core.eval(m, quote
            import BeamletOptics
            struct MyFilter{T, S <: BeamletOptics.AbstractShape{T}} <: BeamletOptics.AbstractObject{T}
                shape::S
                od::T
            end
            BeamletOptics.properties(x::MyFilter) =
                [BeamletOptics.default_properties(x); "Optical density" => x.od]
        end)
        f = Base.invokelatest(m.MyFilter, BMO.CylinderSDF(5e-3, 1e-3), 2.0)
        props = Base.invokelatest(BMO.properties, f)
        @test first.(props)[1:4] == ["Type", "Position [m]", "Optical axis", "Shape"]
        @test prop(props, "Optical density") == 2.0
    end
end

end
