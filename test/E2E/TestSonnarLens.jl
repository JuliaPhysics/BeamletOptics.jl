module TestSonnarLens

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics

const mm = 1e-3

@testset "Sonnar lens" begin
    # Sonnar 50 mm F1.5, source: pencilofrays.com, Sonnar_50mmF1p5_FR837616.zmx,
    # scaled to f = 100 mm
    l1 = SphericalLens(69.21mm, 433.84mm, 9.33mm, 70mm, λ -> 1.671)
    # front triplet: last surface only 40 mm clear aperture -> assembled from individual lenses
    l2 = SphericalLens(35.86mm, 85.87mm, 11.81mm, 60mm, λ -> 1.671)
    l3 = SphericalLens(85.87mm, -646.31mm, 7.05mm, 60mm, λ -> 1.4892)
    l4 = Lens(SphericalSurface(-646.31mm, 60mm), SphericalSurface(23.51mm, 40mm), 1.9mm, λ -> 1.7394)
    translate3d!(l3, [0, thickness(l2), 0])
    translate3d!(l4, [0, thickness(l2) + thickness(l3), 0])
    l234 = TripletLens(l2, l3, l4)
    l567 = SphericalTripletLens(Inf, 51.09mm, -22.12mm, -103.13mm, 2.48mm, 19.81mm, 4.57mm, 42mm,
        λ -> 1.5232, λ -> 1.6578, λ -> 1.5894)
    l_234 = thickness(l1) + 0.38mm
    l_567 = l_234 + thickness(l234) + 13.0mm + 2.24mm
    translate3d!(l234, [0, l_234, 0])
    translate3d!(l567, [0, l_567, 0])
    y_end = l_567 + thickness(l567)

    λ = 587.6e-9

    """Axis crossing y-coordinate of the last ray of `beam`."""
    function axis_crossing(beam)
        r = last(rays(beam))
        p = position(r)
        d = direction(r)
        return p[2] - p[3] / d[3] * d[2]
    end

    system = System([l1, l234, l567])

    @testset "Back focal length at h = 0.01 mm" begin
        beam = Beam(Ray([0, -0.05, 0.01mm], [0, 1, 0], λ))
        solve_system!(system, beam)
        y_c = axis_crossing(beam)
        @test abs((y_c - y_end) - 44.902mm) ≤ 5e-6
    end

    @testset "Ray count for h in -33mm:0.25mm:33mm" begin
        for h in -33mm:0.25mm:33mm
            beam = Beam(Ray([0, -0.05, h], [0, 1, 0], λ))
            solve_system!(system, beam)
            @test length(rays(beam)) == 11
        end
    end

    @testset "Back focal length vs. reference (spherical aberration)" begin
        references = Dict(0.01mm => 44.902mm, 5mm => 44.86mm, 15mm => 44.591mm)
        for (h, bfl_ref) in references
            beam = Beam(Ray([0, -0.05, h], [0, 1, 0], λ))
            solve_system!(system, beam)
            y_c = axis_crossing(beam)
            @test abs((y_c - y_end) - bfl_ref) ≤ 0.01mm
        end
    end
end

end # MODULE
