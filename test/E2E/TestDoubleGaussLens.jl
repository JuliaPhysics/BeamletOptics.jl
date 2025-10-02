module TestDoubleGaussLens

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics

const mm = 1e-3

@testset "Double Gauss lens" begin
    """Tests hit distance to detector origin"""
    function test_coma(pd::Detector; atol=7e-5)
        hits = spot_diagram(pd)
        for hit in hits
            if norm(hit) > atol
                return false
            end
        end
        return true
    end

    # Based on https://www.pencilofrays.com/double-gauss-sonnar-comparison/
    l1 = SphericalLens(48.88mm, 182.96mm, 8.89mm, 52.3mm, λ -> 1.62286)
    l23 = SphericalDoubletLens(36.92mm, Inf, 23.06mm, 15.11mm, 2.31mm,
    45.11mm, λ -> 1.58565, λ -> 1.67764)
    l45 = SphericalDoubletLens(-23.91mm, Inf, -36.92mm, 1.92mm, 7.77mm,
    40.01mm, λ -> 1.57046, λ -> 1.64128)
    l6 = SphericalLens(1063.24mm, -48.88mm, 6.73mm, 45.11mm, λ -> 1.62286)
    # Calculate translation distances
    l_23 = thickness(l1) + 0.38mm
    l_45 = l_23 + thickness(l23) + 9.14mm + 13.36mm
    l_6 = l_45 + thickness(l45) + 0.38mm
    # Corresponds to back focal length of f=59.21 mm on y-axis from link above + "error" δf
    δf = 7e-4
    f_z = l_6 + thickness(l6.shape) + 58.21mm + δf
    translate3d!(l23, [0, l_23, 0])
    translate3d!(l45, [0, l_45, 0])
    translate3d!(l6, [0, l_6, 0])
    # Create and move group - this tests a bunch of kinematic correctness
    double_gauss = ObjectGroup([l1, l23, l45, l6])

    # Spawn detector, move to "focus"
    detector = Detector(5mm)
    translate3d!(detector, [0,f_z,0])
    # Setup test, rotate and translate to test correct kinematics
    test_setup = ObjectGroup([double_gauss, detector])
    translate3d!(test_setup, [0.05, 0.05, 0.05])
    xrotate3d!(test_setup, deg2rad(60))
    zrotate3d!(test_setup, deg2rad(45))
    system = System([test_setup])

    @testset "Test with collimated source" begin
        ## Test against back focal length as per source above
        dir = orientation(double_gauss)[:, 2] # rotated collimated ray direction
        pos = position(l1) - 0.05 * dir # rotated collimated ray position
        source = CollimatedSource(pos, dir, 0.04, 486.0e-9, num_rays=1000, num_rings=10)
        solve_system!(system, source)
        @test test_coma(detector, atol=2e-5)
    end

    # Reset test scenario back to origin, y-axis alignment
    reset_rotation3d!(test_setup)
    reset_translation3d!(test_setup)
    translate_to3d!(detector, [0, 0.147, 0])

    @testset "Test with point source (wide)" begin
        source = PointSource([0, -0.5, 0], [0, 1 ,0], deg2rad(2), 486.0e-9, num_rays=1000, num_rings=10)
        empty!(detector)
        solve_system!(system, source)
        @test test_coma(detector, atol=6e-5)
    end

    @testset "Test with point source (narrow)" begin
        # Tests regression for https://github.com/JuliaPhysics/BeamletOptics.jl/issues/11
        source = PointSource([0, -0.5, 0], [0, 1 ,0], 5e-5, 486.0e-9, num_rays=1000, num_rings=10)
        empty!(detector)
        solve_system!(system, source)
        @test test_coma(detector, atol=2e-7)
    end
end

end # MODULE