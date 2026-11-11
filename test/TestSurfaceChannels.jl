module TestSurfaceChannels

using BeamletOptics
using Test
using LinearAlgebra
using AbstractTrees

const BMO = BeamletOptics
const mm = 1e-3
const λ0 = 1e-6

# The transition of the first hit of the `ray` in the `system`
function first_transition(system, ray)
    return BMO.with_bounding_spheres(system) do
        BMO.intersection!(ray, BMO.first_to_interact(BMO.trace_all(system, ray), ray))
        BMO.transition(system, ray)
    end
end

# A 50:50 thin beamsplitter at 45° in air, at y = 50 mm
function splitter_system()
    bs = ThinBeamsplitter(25mm)
    zrotate3d!(bs, deg2rad(45))
    translate3d!(bs, [0, 50mm, 0])
    return System([bs])
end

@testset "Channels of a bare boundary" begin
    n = 1.5
    plate = RectangularCompensatorPlate(40mm, 40mm, 5mm, ConstantRefractiveIndex(n))
    translate3d!(plate, [0, 50mm, 0])
    system = System([plate])

    ray = Ray([0, 0, 0], [0, 1, 0], λ0)
    transmitted, reflected = BMO.channels(first_transition(system, ray), ray)
    @test BMO.side(transmitted) isa BMO.Transmitted
    @test BMO.side(reflected) isa BMO.Reflected
    @test BMO.direction(transmitted) ≈ [0, 1, 0]
    @test BMO.direction(reflected) ≈ [0, -1, 0]
    # A ray without a field continues without Fresnel loss
    @test BMO.amplitude(transmitted) == 1
    @test BMO.amplitude(reflected) == 0

    # Fresnel amplitudes at normal incidence, Hecht eq. 4.47-4.48
    E = [1.0 + 0im, 0, 0]
    pray = PolarizedRay([0, 0, 0], [0, 1, 0], λ0, E)
    pt, pr = BMO.channels(first_transition(system, pray), pray)
    @test norm(pt.P * E) ≈ 2 / (1 + n)
    @test norm(pr.P * E) ≈ (n - 1) / (n + 1)
    # Uncoated traces the dominant channel: the transmitted one
    @test BMO.select_channels(BMO.Dominant(), (pt, pr)) === pt
end

@testset "Splitting keeps the power" begin
    system = splitter_system()

    beam = Beam(PolarizedRay([0, 0, 0], [0, 1, 0], λ0, [1.0 + 0im, 0, 0]))
    solve_system!(system, beam)
    @test length(children(beam)) == 2
    @test sum(BMO.power_fraction(c, beam) for c in children(beam)) ≈ 1

    for beamlet in (GaussianBeamlet([0, 0, 0], [0, 1, 0], λ0, 1mm),
            AstigmaticGaussianBeamlet([0, 0, 0], [0, 1, 0], λ0, 1mm, 0.5mm))
        solve_system!(system, beamlet)
        @test length(children(beamlet)) == 2
        @test sum(optical_power, children(beamlet)) ≈ optical_power(beamlet)
    end
end

@testset "The power cutoff bounds the beam tree of a cavity" begin
    # Two partial mirrors at normal incidence, R = 0.9: each pass loses a factor of 0.9 or 0.1
    mirrors = [ThinBeamsplitter(25mm; reflectance = 0.9) for _ in 1:2]
    translate3d!(mirrors[1], [0, 50mm, 0])
    translate3d!(mirrors[2], [0, 60mm, 0])
    system = System(mirrors)
    num_beams(beam) = count(_ -> true, PreOrderDFS(beam))
    new_beam() = Beam(PolarizedRay([0, 0, 0], [0, 1, 0], λ0, [1.0 + 0im, 0, 0]))
    for beam in (new_beam(), GaussianBeamlet([0, 0, 0], [0, 1, 0], λ0, 1mm))
        solve_system!(system, beam; depth_max = 30)
        all_beams = num_beams(beam)
        solve_system!(system, beam; depth_max = 30, power_cutoff = 1e-3)
        @test 1 < num_beams(beam) < all_beams
        @test all(b -> BMO.power_fraction(b, beam) ≥ 1e-3, PreOrderDFS(beam))
    end
    # The tree of rays without a field is bounded by depth_max only
    rays_only = Beam(Ray([0, 0, 0], [0, 1, 0], λ0))
    solve_system!(system, rays_only; depth_max = 6)
    n = num_beams(rays_only)
    solve_system!(system, rays_only; depth_max = 6, power_cutoff = 0.5)
    @test num_beams(rays_only) == n
    for cutoff in (-0.1, 1.5)
        @test_throws ArgumentError solve_system!(system, new_beam(); power_cutoff = cutoff)
    end
end

# A dummy surface model that implements only its channels: a lossy half mirror that transmits and
# reflects a quarter of the power each
struct HalfMirror <: BMO.AbstractSurfaceModel end
BMO.splitting_of(::HalfMirror) = BMO.AllChannels()
BMO.channels(::HalfMirror, ::BMO.Transition, ray::BMO.AbstractRay, (t, r)) =
    (BMO.transfer(t, ray, BMO.SPBasis(0.5, 0, 0, 0.5), 0.5),
        BMO.transfer(r, ray, BMO.SPBasis(-0.5, 0, 0, 0.5), 0.5))

# A dummy thin layer that carries the half mirror
struct HalfMirrorFilm{S <: BMO.AbstractShape{Float64}} <: BMO.AbstractObject{Float64}
    shape::S
end
BMO.surface_model(::HalfMirrorFilm) = HalfMirror()

@testset "A surface model needs only its channels" begin
    film = HalfMirrorFilm(BMO.RectangularFlatMesh(25mm, 25mm))
    zrotate3d!(film, deg2rad(45))
    translate3d!(film, [0, 50mm, 0])
    system = System([film])

    beam = Beam([0, 0, 0], [0, 1, 0], λ0)
    solve_system!(system, beam)
    @test length(children(beam)) == 2
    @test BMO.direction(first(rays(children(beam)[1]))) ≈ [0, 1, 0]
    @test abs(BMO.direction(first(rays(children(beam)[2])))[1]) ≈ 1

    gauss = GaussianBeamlet([0, 0, 0], [0, 1, 0], λ0, 1mm)
    solve_system!(system, gauss)
    @test length(children(gauss)) == 2
    @test all(optical_power(c) ≈ optical_power(gauss) / 4 for c in children(gauss))

    agb = AstigmaticGaussianBeamlet([0, 0, 0], [0, 1, 0], λ0, 1mm, 1mm)
    solve_system!(system, agb)
    @test sum(optical_power, children(agb)) ≈ optical_power(agb) / 2
end

@testset "The hooks of bodies, thin layers and surface models are public" begin
    for name in (:AbstractMedium, :IsotropicMedium, :medium, :AbstractSurfaceModel, :splitting_of,
            :Transition, :transition, :channels, :BoundaryChannel, :transfer, :amplitude)
        @test Base.ispublic(BeamletOptics, name)
    end
    for name in (:Uncoated, :Dominant, :AllChannels)
        @test Base.isexported(BeamletOptics, name)
    end
end

end # module
