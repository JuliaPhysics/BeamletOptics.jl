module TestDetectorSampling

using BeamletOptics
using Test

const BMO = BeamletOptics

const mm = 1e-3

# A sampling of one's own: points of any arrangement and their weights, as in #168
struct PointCloud{P, W} <: BMO.AbstractDetectorSampling
    points::P
    weights::W
end

BMO.sample_points(::Detector, s::PointCloud) = BMO.DetectorSamples((s.points,), s.points, s.weights)

# A sampling without `sample_points`
struct NoPoints <: BMO.AbstractDetectorSampling end

# A Gaussian beam of 1 mW on a detector tilted by 30°
function tilted_detector()
    pd = Detector(50mm)
    zrotate3d!(pd, deg2rad(30))
    translate3d!(pd, [0, 100mm, 0])
    solve_system!(System([pd]), GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 1e-6, 1mm; P0 = 1e-3))
    return pd
end

@testset "Detector sampling" begin
    pd = tilted_detector()

    @testset "RegularGrid gives the values of the former keyword readouts" begin
        # recorded with `intensity(pd; n = 50)`, `optical_power(pd; ...)` before the samplings
        x, z, I = intensity(pd, RegularGrid(num_points = 50); progress = false)
        @test sum(I) ≈ 49970.20330881326 rtol = 1e-12
        @test I[20, 27] ≈ 203.20885922342407 rtol = 1e-12
        @test first(x) ≈ -0.0034658561053307387 rtol = 1e-12
        @test last(x) ≈ 0.0034658561053307387 rtol = 1e-12
        @test optical_power(pd, RegularGrid(num_points = 50); progress = false) ≈ 0.0010000000061118099 rtol = 1e-12
        @test optical_power(pd; progress = false) ≈ 0.0010000000063770916 rtol = 1e-12
        s = RegularGrid(num_points = 31, crop_factor = 2, x_min = -1e-3, x_max = 2e-3)
        x, z, I = intensity(pd, s; progress = false)
        @test sum(I) ≈ 62674.78458469183 rtol = 1e-12
        @test (first(x), last(x)) == (-1e-3, 2e-3)
        @test last(z) ≈ 0.002310570736887159 rtol = 1e-12
    end

    @testset "A single bound replaces only its side" begin
        x_auto, z_auto, _ = electric_field(pd, RegularGrid(num_points = 11); progress = false)
        x, z, _ = electric_field(pd, RegularGrid(num_points = 11, x_min = -1mm); progress = false)
        @test first(x) == -1mm
        @test last(x) == last(x_auto)
        @test z == z_auto
        x, z, _ = electric_field(pd, RegularGrid(num_points = 11, z_max = 2mm, x0_shift = 1mm); progress = false)
        @test x ≈ x_auto .+ 1mm
        @test (first(z), last(z)) == (first(z_auto), 2mm)
    end

    @testset "RectilinearGrid" begin
        # non-uniform: dense at the edges, coarse in the center
        xs = 4mm .* sinpi.(range(-0.5, 0.5, length = 601))
        zs = 3mm .* sinpi.(range(-0.5, 0.5, length = 401))
        g = RectilinearGrid(xs, zs)
        x, z, E = electric_field(pd, g; progress = false)
        @test x === xs && z === zs
        @test size(E) == (601, 401)
        @test optical_power(pd, g; progress = false) ≈ 1e-3 rtol = 1e-4
        # a single node has no area
        @test optical_power(pd, RectilinearGrid([0.0], [0.0]); progress = false) == 0
        # the weights stay lazy and add up to the area of the window
        @test BMO.sample_points(pd, g) isa BMO.DetectorSamples
        w = BMO.sample_points(pd, g).weights
        @test !(w isa AbstractArray)
        @test sum(w) ≈ 8mm * 6mm rtol = 1e-12
        # a wide and a tall grid give the same values as their transposes
        xs1, zs1 = range(-3mm, 3mm, length = 3), range(-2mm, 2mm, length = 301)
        I_tall = intensity(pd, RectilinearGrid(zs1, xs1); progress = false)[3]
        I_wide = intensity(pd, RectilinearGrid(xs1, zs1); progress = false)[3]
        @test size(I_tall) == (301, 3) && size(I_wide) == (3, 301)
    end

    @testset "A sampling of one's own" begin
        g = RectilinearGrid(range(-3mm, 3mm, length = 41), range(-2mm, 2mm, length = 31))
        grid = BMO.sample_points(pd, g)
        cloud = PointCloud(vec(grid.points), vec(copy(grid.weights)))
        pts, I = intensity(pd, cloud; progress = false)
        @test pts === cloud.points
        @test I isa Vector
        @test I ≈ vec(intensity(pd, g; progress = false)[3])
        @test optical_power(pd, cloud; progress = false) ≈ sum(cloud.weights .* I)
        @test optical_power(pd, cloud; progress = false) ≈ optical_power(pd, g; progress = false)
        _, E = electric_field(pd, cloud; progress = false)
        @test E ≈ vec(electric_field(pd, g; progress = false)[3])
        # the values at local points alone, without a sampling
        @test intensity(pd, cloud.points; progress = false) == I
        @test electric_field(pd, cloud.points; progress = false) == E
    end

    @testset "Errors" begin
        @test_throws ArgumentError RegularGrid(num_points = 0)
        @test_throws "sample_points not defined for $(NoPoints)" intensity(pd, NoPoints())
        # a bound beyond its opposite side would give a descending grid and a negative power
        x_auto, z_auto, _ = intensity(pd, RegularGrid(num_points = 5); progress = false)
        @test_throws ArgumentError intensity(pd, RegularGrid(x_min = last(x_auto) + 1mm))
        @test_throws ArgumentError optical_power(pd, RegularGrid(z_max = first(z_auto) - 1mm))
        @test_throws ArgumentError RectilinearGrid([1.0, 0.0], [0.0, 1.0])
        @test_throws ArgumentError RectilinearGrid(Float64[], [0.0])
        # a single point is not an array of points
        @test_throws MethodError intensity(pd, BMO.Point2(0.0, 0.0))
    end
end

@testset "A detector without hits reads as nothing" begin
    # A beam that misses a tilted detector, as in a scan (#189)
    empty = Detector(10mm)
    zrotate3d!(empty, deg2rad(30))
    translate3d!(empty, [20mm, 100mm, 0])
    solve_system!(System([empty]), GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 1e-6, 1mm))
    @test length(empty) == 0
    @test spot_diagram(empty) == BMO.Point2{Float64}[]
    @test optical_power(empty; progress = false) == 0
    # The default window is the active surface
    x, z, I = intensity(empty, RegularGrid(num_points = 11); progress = false)
    @test [first(x), last(x), first(z), last(z)] ≈ [-5mm, 5mm, -5mm, 5mm] rtol = 1e-12
    @test size(I) == (11, 11) && all(iszero, I)
    x, z, E = electric_field(empty, RegularGrid(num_points = 11, x_min = -1mm); progress = false)
    @test first(x) == -1mm
    @test eltype(E) == ComplexF64 && all(iszero, E)
    g = RectilinearGrid(range(-1mm, 1mm, length = 5), range(0, 2mm, length = 3))
    @test optical_power(empty, g; progress = false) == 0
    @test all(iszero, intensity(empty, g; progress = false)[3])
    # A detector that was hit and emptied, e.g. by moving it off the beam
    pd = tilted_detector()
    @test optical_power(pd; progress = false) > 0
    translate3d!(pd, [20mm, 0, 0])
    @test optical_power(pd; progress = false) == 0
    @test isempty(spot_diagram(pd))
end

end # module
