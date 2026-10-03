module TestClearAperture

using BeamletOptics
using Test

const BMO = BeamletOptics

const mm = 1e-3

@testset "Clear aperture" begin
    pos = [0.0, -0.1, 0.0]
    dir = [0.0, 1.0, 0.0]
    # weak lens followed by a smaller one: the second lens limits the bundle
    l1 = SphericalLens(0.5, -0.5, 4mm, 40mm, 1.5)
    l2 = SphericalLens(0.5, -0.5, 4mm, 10mm, 1.5)
    translate3d!(l2, [0, 0.05, 0])
    system = System([l1, l2])

    @testset "single lens: lens diameter" begin
        D = clear_aperture(System([l1]), pos, dir)
        @test D ≈ 40mm rtol = 2e-3
        # the search limit is returned if the bundle is not vignetted
        @test clear_aperture(System([l1]), pos, dir; d_max = 10mm) == 10mm
        # an axial ray that does not hit anything has no clear aperture
        @test_throws ArgumentError clear_aperture(System([l1]), pos .+ [0.1, 0, 0], dir)
    end

    @testset "second smaller lens: analytic clear aperture" begin
        # the radius of a marginal ray at the second lens follows from tracing it
        b = Beam([5mm, pos[2], 0], dir)
        solve_system!(system, b)
        ρ = abs(position(rays(b)[4])[1]) / 5mm
        @test BMO.object(BMO.intersection(rays(b)[3])) === l2
        D = clear_aperture(system, pos, dir)
        @test D ≈ 10mm / ρ rtol = 5e-3
        # the sampled directions do not matter: other azimuth sampling and a rotated axis
        @test clear_aperture(system, pos, dir; rings = 3, azimuths = 7) ≈ D rtol = 5e-3
        @test clear_aperture(system, pos, dir; λ = 500e-9) ≈ D rtol = 5e-3

        # vignetted beams of a source that is slightly larger or smaller
        axis = Beam(pos, dir)
        solve_system!(system, axis)
        for (scale, nvig) in ((0.95, false), (1.05, true))
            src = UniformDiscSource(pos, dir, scale * D; num_rays = 200)
            solve_system!(system, src; progress = false)
            @test isempty(vignetted(src, axis)) == !nvig
        end
        src = UniformDiscSource(pos, dir, 40mm; num_rays = 200)
        solve_system!(system, src; progress = false)
        idx = vignetted(src, axis)
        @test !isempty(idx)
        @test all(i -> hypot(position(src[i])[1], position(src[i])[3]) > 0.9 * D / 2, idx)
    end

    @testset "geometry-agnostic: tilted axis" begin
        # the same system rotated about z: no assumed world orientation
        s = System([l1, l2])
        D0 = clear_aperture(s, pos, dir)
        θ = deg2rad(30)
        rotate3d!(l1, [0, 0, 1], θ, [0, 0, 0])
        rotate3d!(l2, [0, 0, 1], θ, [0, 0, 0])
        R = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
        @test clear_aperture(s, R * pos, R * dir) ≈ D0 rtol = 5e-3
    end

    @testset "mechanical rim is not clear aperture" begin
        # the annular rim ring is entered and exited by the same lens, like the optical faces
        lens = Lens(BMO.SphericalSurface{Float64}(0.5, 20mm, 40mm),
            BMO.SphericalSurface{Float64}(-0.5, 20mm, 40mm), 4mm, λ -> 1.5)
        s = System([lens])
        axis = Beam(pos, dir)
        solve_system!(s, axis)
        rim = Beam([15mm, pos[2], 0], dir)
        solve_system!(s, rim)
        # same object sequence, different boundary parts
        @test all(r -> BMO.object(BMO.intersection(r)) === lens, rays(rim)[1:2])
        @test BMO.hit_sequence(rim) != BMO.hit_sequence(axis)
        @test clear_aperture(s, pos, dir) ≈ 20mm rtol = 5e-3
    end

    @testset "detectors are not modified" begin
        det = Detector(0.2)
        translate3d!(det, [0, 0.2, 0])
        s = System([l1, det])
        D0 = clear_aperture(System([l1]), pos, dir)
        # a fresh detector stays empty
        @test clear_aperture(s, pos, dir) ≈ D0 rtol = 5e-3
        @test isnothing(BMO.hits(det))
        # existing beamlet hits are neither overwritten nor does the scalar probe throw
        solve_system!(s, GaussianBeamlet(pos, dir, 1e-6, 1mm))
        old = BMO.hits(det)
        @test old isa Vector{<:BMO.GaussianBeamletHit}
        @test clear_aperture(s, pos, dir) ≈ D0 rtol = 5e-3
        @test BMO.hits(det) === old
        @test length(old) == 1
        # also for detectors inside object groups
        det2 = Detector(0.2)
        translate3d!(det2, [0, 0.2, 0])
        sg = System([l1, ObjectGroup([det2])])
        solve_system!(sg, GaussianBeamlet(pos, dir, 1e-6, 1mm))
        old2 = BMO.hits(det2)
        @test !isnothing(old2)
        @test clear_aperture(sg, pos, dir) ≈ D0 rtol = 5e-3
        @test BMO.hits(det2) === old2
        # the state is restored if the function throws
        @test_throws ArgumentError clear_aperture(s, pos .+ [0.3, 0, 0], dir)
        @test BMO.hits(det) === old
    end

    @testset "beamlet groups" begin
        # the chief rays of Gaussian beamlets are compared
        axis = AstigmaticGaussianBeamlet(pos, dir, 1e-6, 1mm)
        solve_system!(System([l1]), axis)
        src = CollimatedGaussianBeamletSource(pos, dir, 60mm, 1e-6, 1mm; n_grid = 10)
        solve_system!(System([l1]), src; progress = false)
        idx = vignetted(src, axis)
        @test !isempty(idx)
        @test length(idx) < length(src)
        @test all(i -> max(abs(position(src[i])[1]), abs(position(src[i])[3])) > 15mm, idx)
        g = GaussianBeamlet(pos, dir, 1e-6, 1mm)
        solve_system!(System([l1]), g)
        @test length(BMO.hit_sequence(g)) == length(BMO.hit_sequence(axis))
    end

    @testset "argument checks" begin
        @test_throws ArgumentError clear_aperture(system, pos, dir; rings = 0)
        @test_throws ArgumentError clear_aperture(system, pos, dir; azimuths = 2)
        @test_throws ArgumentError clear_aperture(system, pos, dir; d_start = 0)
        @test_throws ArgumentError clear_aperture(system, pos, dir; d_start = 1, d_max = 0.5)
    end
end

end # module
