module TestAstigmaticGaussianSources

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics
const mm = 1e-3
const nm = 1e-9

@testset "Astigmatic Gaussian Beamlet Sources" begin
    λ = 532nm
    dir = [0.0, 1.0, 0.0]
    pos = [0.0, 0.0, 0.0]

    @testset "CollimatedGaussianBeamletSource" begin
        D = 2mm
        w0s = D / 10
        bg = CollimatedGaussianBeamletSource(pos, dir, D, λ, w0s; n_grid = 10)
        @test bg isa AstigmaticBeamGroup
        beams = BMO.beams(bg)

        @test length(beams) == 100
        # Check all point in the right direction
        @test all(normalize(BMO.direction(b)) ≈ dir for b in bg)
        # Check all have the correct sub-waist
        for b in beams
            _, w1, _ = BMO.waist_parameters(b, 0.0)
            @test norm(w1) ≈ w0s
        end
    end

    @testset "SphericalGaussianBeamletSource" begin
        θ = 0.1 # radians
        bg = SphericalGaussianBeamletSource(pos, dir, θ, λ; num_rings = 5, num_rays = 500)
        @test bg isa AstigmaticBeamGroup

        @test length(bg) == 500
        # Check they all originate from pos
        @test all(BMO.position(b) ≈ pos for b in bg)

        @testset "basis kwarg" begin
            dir_n = normalize(dir)
            b0 = BMO.normal3d(dir_n)
            rotated = SphericalGaussianBeamletSource(pos, dir, θ, λ; num_rings = 5,
                num_rays = 500, basis = BMO.rotate3d(dir_n, deg2rad(30)) * b0)

            dirs = [BMO.direction(b) for b in bg]
            dirs_rot = [BMO.direction(b) for b in rotated]
            # spinning the basis about dir leaves the polar angle distribution intact
            @test sort(round.(BMO.angle3d.(Ref(dir_n), dirs), digits = 11)) ≈
                  sort(round.(BMO.angle3d.(Ref(dir_n), dirs_rot), digits = 11))
            @test !all(dirs .≈ dirs_rot)
            # passing the default basis reproduces the default
            @test all(dirs .≈ [BMO.direction(b) for b in SphericalGaussianBeamletSource(
                pos, dir, θ, λ; num_rings = 5, num_rays = 500, basis = b0)])
            # a degenerate basis would silently yield NaN directions, so it throws
            @test_throws ArgumentError SphericalGaussianBeamletSource(pos, dir, θ, λ;
                num_rings = 5, num_rays = 500, basis = dir)
        end

        @testset "invalid num_rays and θ" begin
            @test_throws ArgumentError SphericalGaussianBeamletSource(pos, dir, θ, λ;
                num_rings = 5, num_rays = 99)
            @test_throws ArgumentError SphericalGaussianBeamletSource(pos, dir, π, λ;
                num_rings = 5, num_rays = 500)
        end
    end

    @testset "EllipticalGaussianBeamletSource" begin
        θ_x = 0.2
        θ_y = 0.1
        bg = EllipticalGaussianBeamletSource(pos, dir, θ_x, θ_y, λ; num_rings = 5, num_rays = 500)
        @test bg isa AstigmaticBeamGroup

        @test length(bg) == 500
        # Check they all originate from pos
        @test all(BMO.position(b) ≈ pos for b in bg)

        # Check directions are within the elliptical cone
        dir_n = normalize(dir)
        b1 = BMO.normal3d(dir_n)
        b2 = cross(dir_n, b1)
        
        for b in bg
            d = BMO.direction(b)
            # project on tangent plane using local basis
            u = dot(d, b1) / dot(d, dir_n)
            v = dot(d, b2) / dot(d, dir_n)
            @test (u / tan(θ_x))^2 + (v / tan(θ_y))^2 <= 1.05
        end

        @testset "basis kwarg" begin
            b0 = BMO.normal3d(dir_n)
            # the elliptical cone is not rotationally symmetric, so a rotated basis
            # genuinely reorients the pattern rather than just relabelling it
            rotated = EllipticalGaussianBeamletSource(pos, dir, θ_x, θ_y, λ; num_rings = 5,
                num_rays = 500, basis = BMO.rotate3d(dir_n, deg2rad(30)) * b0)
            dirs = [BMO.direction(b) for b in bg]
            dirs_rot = [BMO.direction(b) for b in rotated]
            @test !all(dirs .≈ dirs_rot)
            # passing the default basis reproduces the default
            @test all(dirs .≈ [BMO.direction(b) for b in EllipticalGaussianBeamletSource(
                pos, dir, θ_x, θ_y, λ; num_rings = 5, num_rays = 500, basis = b0)])
            # a degenerate basis would silently yield NaN directions, so it throws
            @test_throws ArgumentError EllipticalGaussianBeamletSource(pos, dir, θ_x, θ_y, λ;
                num_rings = 5, num_rays = 500, basis = dir)
        end

        @testset "invalid num_rays and θ" begin
            @test_throws ArgumentError EllipticalGaussianBeamletSource(pos, dir, θ_x, θ_y, λ;
                num_rings = 5, num_rays = 99)
            @test_throws ArgumentError EllipticalGaussianBeamletSource(pos, dir, π / 2, θ_y, λ;
                num_rings = 5, num_rays = 500)
            @test_throws ArgumentError EllipticalGaussianBeamletSource(pos, dir, θ_x, π / 2, λ;
                num_rings = 5, num_rays = 500)
        end
    end

    @testset "AstigmaticBeamGroup from dir and from orientation" begin
        bs = BMO.beams(SphericalGaussianBeamletSource(pos, dir, 0.1, λ; num_rings = 2, num_rays = 40))
        p = [1.0, 2.0, 3.0]
        d = normalize([0.3, 1.0, -0.2])
        e1 = BMO.normal3d(d)
        g_dir = AstigmaticBeamGroup(bs, p, d)
        g_mat = AstigmaticBeamGroup(bs, p, hcat(e1, d, cross(e1, d)))
        @test BMO.position(g_dir) == BMO.position(g_mat) == p
        @test BMO.orientation(g_dir) ≈ BMO.orientation(g_mat) atol = 1e-15
        @test BMO.beams(g_dir) === bs
        @test BMO.beams(g_mat) === bs
    end
end

@testset "Power and profile of beamlet sources" begin
    # the power of a source on a detector that captures all of it
    function power(source, L, half_width)
        pd = Detector(2half_width)
        translate3d!(pd, [0, L, 0])
        solve_system!(System([pd]), source; progress = false)
        return BeamletOptics.optical_power(pd, RegularGrid(num_points = 200); progress = false)
    end
    P0 = 1e-3
    @testset "SphericalGaussianBeamletSource emits P0" begin
        for (θ, num_rings, num_rays) in ((0.02, 5, 200), (0.02, 10, 1000))
            source = SphericalGaussianBeamletSource([0.0, 0, 0], [0.0, 1, 0], θ, 1e-6; num_rings, num_rays, P0)
            @test power(source, 0.2, 3θ * 0.2) ≈ P0 rtol = 1e-2
        end
    end
    @testset "EllipticalGaussianBeamletSource emits P0" begin
        source = EllipticalGaussianBeamletSource([0.0, 0, 0], [0.0, 1, 0], 0.02, 0.01, 1e-6; num_rings = 10, num_rays = 1000, P0)
        @test power(source, 0.2, 3 * 0.02 * 0.2) ≈ P0 rtol = 1e-2
    end
    @testset "CollimatedGaussianBeamletSource emits P0" begin
        for n_grid in (5, 10)
            D = 10e-3
            source = CollimatedGaussianBeamletSource([0.0, 0, 0], [0.0, 1, 0], D, 1e-6, 1.2D / n_grid; n_grid, P0)
            @test power(source, 0.05, D) ≈ P0 rtol = 1e-2
        end
    end
    @testset "E0 sets the polarization, P0 the power" begin
        # E0 of any magnitude, linear or elliptical
        for E0 in ([3.0, 0, 0], [7.0, 0, 3.5im])
            spherical = SphericalGaussianBeamletSource([0.0, 0, 0], [0.0, 1, 0], 0.02, 1e-6; num_rings = 10, num_rays = 1000, P0, E0)
            @test power(spherical, 0.2, 0.012) ≈ P0 rtol = 1e-2
            elliptical = EllipticalGaussianBeamletSource([0.0, 0, 0], [0.0, 1, 0], 0.02, 0.01, 1e-6; num_rings = 10, num_rays = 1000, P0, E0)
            @test power(elliptical, 0.2, 0.012) ≈ P0 rtol = 1e-2
            decomposition = GaussianBeamletDecomposition([0.0, 0, 0], [0.0, 1, 0], 1e-6, 2e-3; P0, E0)
            @test power(decomposition, 0.1, 0.012) ≈ P0 rtol = 1e-2
            # the polarization of a beamlet is that of E0
            p = BeamletOptics.electric_field(first(BeamletOptics.beams(spherical)))
            @test p / norm(p) ≈ E0 / norm(E0)
        end
    end
    @testset "num_rings counts the center as a ring" begin
        # the center and num_rings - 1 rings around it, as for PointSource
        dir = [0.0, 1, 0]
        spherical = SphericalGaussianBeamletSource([0.0, 0, 0], dir, 0.02, 1e-6; num_rings = 5, num_rays = 500)
        angles = unique(round.([acos(clamp(dot(BeamletOptics.direction(b), dir), -1, 1)) for b in BeamletOptics.beams(spherical)], digits = 6))
        @test length(angles) == 5
        θx, θy = 0.02, 0.01
        elliptical = EllipticalGaussianBeamletSource([0.0, 0, 0], dir, θx, θy, 1e-6; num_rings = 5, num_rays = 500)
        b1 = BeamletOptics.orientation(elliptical)[:, 1]
        b2 = normalize(cross(dir, b1))
        ρ(d) = (u = dot(d, b1) / dot(d, dir); v = dot(d, b2) / dot(d, dir); sqrt((u / tan(θx))^2 + (v / tan(θy))^2))
        radii = unique(round.([ρ(BeamletOptics.direction(b)) for b in BeamletOptics.beams(elliptical)], digits = 6))
        @test sort(radii) ≈ collect(LinRange(0, 1, 5))
    end
end

end # MODULE
