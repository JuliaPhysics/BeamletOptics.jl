module TestBeams

using BeamletOptics
using Test

const BMO = BeamletOptics

@testset "Beams" begin
    is = BMO.Intersection(1.0, zeros(3))
    r1 = Ray([0.0, 0, 0], [1, 0, 0])
    r2 = Ray([1.0, 0, 0], [0, 1, 0])
    r3 = Ray([1.0, 1, 0], [0, 0, 1])
    r4 = Ray([1.0, 1, 1], [1, 0, 0])
    BMO.intersection!(r1, is)
    BMO.intersection!(r2, is)
    BMO.intersection!(r3, is)
    BMO.intersection!(r4, nothing)
    # Test beam
    beam = Beam(r1)
    push!(beam, r2)
    push!(beam, r3)
    push!(beam, r4)
    @test length(beam) == 3
    @test BMO.point_on_beam(beam, 0) == ([0, 0, 0], 1)
    @test BMO.point_on_beam(beam, 1) == ([1, 0, 0], 2)
    @test BMO.point_on_beam(beam, 2) == ([1, 1, 0], 3)
    @test BMO.point_on_beam(beam, 3) == ([1, 1, 1], 4)
    @test BMO.point_on_beam(beam, 10) == ([8, 1, 1], 4)
    @test BMO.isparentbeam(beam, r2) == true
end

@testset "Rays and first ray of every beam type" begin
    b = Beam([0.0, 0, 0], [0, 1, 0])
    g = GaussianBeamlet([0.0, 0, 0], [0, 1, 0], 1e-6, 1e-3)
    a = AstigmaticGaussianBeamlet([0.0, 0, 0], [0, 1, 0], 1e-6, 1e-3)
    # a beamlet returns the rays of its chief beam
    @test BMO.rays(g) === BMO.rays(g.chief)
    @test BMO.rays(a) === BMO.rays(a.c)
    @test BMO.first_ray(b) === first(BMO.rays(b))
    @test BMO.first_ray(g) === first(BMO.rays(g.chief))
    @test BMO.first_ray(a) === first(BMO.rays(a.c))
end

@testset "Length before a ray" begin
    # beam along +y, folded to +x by m1 at y = 0.1, back to +y by m2 at x = 0.2
    m1 = SquarePlanoMirror2D(25e-3)
    m2 = SquarePlanoMirror2D(25e-3)
    zrotate3d!(m1, deg2rad(45))
    zrotate3d!(m2, deg2rad(45))
    translate3d!(m1, [0, 0.1, 0])
    translate3d!(m2, [0.2, 0.1, 0])
    beam = Beam([0.0, 0, 0], [0, 1, 0])
    solve_system!(System([m1, m2]), beam)
    @test length(BMO.rays(beam)) == 3
    @test BMO.length_before(beam, 1) == 0
    @test BMO.length_before(beam, 2) ≈ 0.1 atol = 1e-12
    @test BMO.length_before(beam, 3) ≈ 0.3 atol = 1e-12
    # the transmitted child behind a beamsplitter at y = 0.05 includes the parent length
    bs = ThinBeamsplitter(25e-3)
    zrotate3d!(bs, deg2rad(45))
    translate3d!(bs, [0, 0.05, 0])
    beam = Beam([0.0, 0, 0], [0, 1, 0])
    solve_system!(System([m1, m2, bs]), beam)
    child = only(filter(c -> BMO.direction(c) ≈ [0, 1, 0], BMO.children(beam)))
    @test length(BMO.rays(child)) == 3
    @test BMO.length_before(child, 1) ≈ 0.05 atol = 1e-12
    @test BMO.length_before(child, 2) ≈ 0.1 atol = 1e-12
    @test BMO.length_before(child, 3) ≈ 0.3 atol = 1e-12
end

@testset "Beam of PolarizedRays with a complex field" begin
    # a complex `E0`, e.g. circular polarization, gives a beam of real number type
    for E0 in ([1.0, 0, 0], [1, 0, im] / sqrt(2), ComplexF32[1, 0, 0], [1, 0, 0])
        b = Beam([0.0, 0, 0], [0.0, 1, 0], 1e-6, E0)
        @test b isa Beam{Float64, PolarizedRay{Float64}}
        @test BMO.polarization(first(BMO.rays(b))) ≈ E0
    end
end

end # MODULE