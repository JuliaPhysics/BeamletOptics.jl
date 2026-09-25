module TestStorageSources

using BeamletOptics
using Test
using Random: MersenneTwister

const BMO = BeamletOptics

tmpfile() = joinpath(mktempdir(), "sources.bmo")

function roundtrip(src)
    path = tmpfile()
    BMO.save_setup(path; sources = [src])
    loaded = BMO.load_setup(path)
    @test length(loaded.sources) == 1
    return only(loaded.sources), path
end

# All component beams of a source, first ray only, in a fixed order
component_beams(b::Beam) = [b]
component_beams(g::GaussianBeamlet) = [g.chief, g.waist, g.divergence]
component_beams(a::AstigmaticGaussianBeamlet) = [getfield(a, f) for f in (:c, BMO._ASTIGMATIC_AUX...)]
component_beams(bg::BMO.AbstractBeamGroup) = reduce(vcat, [component_beams(b) for b in BMO.beams(bg)])

ray_state(r::Ray) = (r.pos, r.dir, r.λ, r.n)
ray_state(r::PolarizedRay) = (r.pos, r.dir, r.λ, r.n, r.E0)

source_state(::Beam) = ()
source_state(g::GaussianBeamlet) = (g.λ, g.w0, g.E0)
source_state(::AstigmaticGaussianBeamlet) = ()
source_state(ps::PointSource) = (ps.NA,)
source_state(cs::CollimatedSource) = (cs.diameter,)
source_state(::AstigmaticBeamGroup) = ()

function test_identical(a, b)
    @test typeof(a) == typeof(b)
    @test source_state(a) == source_state(b)
    ba, bb = component_beams(a), component_beams(b)
    @test length(ba) == length(bb)
    @test all(ray_state(first(rays(x))) == ray_state(first(rays(y))) for (x, y) in zip(ba, bb))
    # Fresh, untraced state
    @test all(length(rays(y)) == 1 && isnothing(BMO.intersection(first(rays(y)))) for y in bb)
    @test all(isnothing(y.parent) && isempty(y.children) for y in bb)
end

function test_system()
    mirror = RoundPlanoMirror(0.05, 0.005)
    zrotate3d!(mirror, deg2rad(40))
    xrotate3d!(mirror, deg2rad(5))
    translate3d!(mirror, [0.001, 0.2, -0.002])
    return System([mirror])
end

traced_state(src) = [[(r.pos, r.dir) for r in rays(b)] for b in component_beams(src)]

@testset "Storage sources" begin
    pos = [0.001, -0.02, 0.003]
    dir = [0.01, 1.0, -0.02]
    λ = 633e-9

    sources = Pair{String, Any}[
        "Beam(Ray)" => Beam(pos, dir, λ),
        "Beam(PolarizedRay)" => Beam(PolarizedRay(pos, [0, 1.0, 0], λ, [1.0 + 0.5im, 0, 0.3im])),
        "Beam(Ray) with n" => Beam(BMO.Ray{Float64}(BMO.Point3(pos...), BMO.Point3(0.0, 1.0, 0.0), nothing, λ, 1.33)),
        "GaussianBeamlet" => GaussianBeamlet(pos, dir, λ, 1.3e-3; M2 = 1.4, P0 = 2e-3, z0 = 0.05),
        "GaussianBeamlet support" => GaussianBeamlet(pos, [0, 1.0, 0], λ, 0.5e-3; support = [1.0, 0, 0]),
        "AstigmaticGaussianBeamlet" => AstigmaticGaussianBeamlet(pos, [0, 1.0, 0], λ, 1e-3, 2e-3;
            M2_x = 1.2, E0 = [1, 0, 0.5im], z0_x = 0.01, z0_y = -0.02),
        "PointSource" => PointSource(pos, dir, deg2rad(5), λ; num_rings = 5, num_rays = 200),
        "CollimatedSource" => CollimatedSource(pos, dir, 0.01, λ; num_rings = 5, num_rays = 200, basis = [1.0, 0, 0]),
        "UniformDiscSource" => UniformDiscSource(pos, dir, 0.01, λ; num_rays = 300),
        "CollimatedSource(PolarizedRay)" => CollimatedSource(
            [Beam(PolarizedRay([x, 0, 0], [0, 1.0, 0], λ, [0, 0, 1.0 + 0im])) for x in range(-1e-3, 1e-3, 5)], 2e-3),
        "CollimatedGaussianBeamletSource" => CollimatedGaussianBeamletSource(pos, dir, 5e-3, λ, 0.5e-3;
            n_grid = 6, randomize_axes = true, rng = MersenneTwister(1)),
        "GaussianBeamletDecomposition" => GaussianBeamletDecomposition(pos, [0, 1.0, 0], λ, 1e-3;
            n_grid = 7, randomize_axes = true, rng = MersenneTwister(2), E0 = [1, 0, 1im]),
        "SphericalGaussianBeamletSource" => SphericalGaussianBeamletSource(pos, dir, deg2rad(3), λ;
            num_rings = 3, num_rays = 60, randomize_axes = true, rng = MersenneTwister(3)),
        "EllipticalGaussianBeamletSource" => EllipticalGaussianBeamletSource(pos, dir, deg2rad(3), deg2rad(1), λ;
            num_rings = 3, num_rays = 60, randomize_axes = true, rng = MersenneTwister(4)),
        "WavefrontBeamletDecomposition" => let x = range(-1e-3, 1e-3, 9), y = range(-1e-3, 1e-3, 9)
            amp = [exp(-(xi^2 + yi^2) / (0.6e-3)^2) for xi in x, yi in y]
            phase = [1e3 * xi + 5e5 * yi^2 for xi in x, yi in y]
            WavefrontBeamletDecomposition(collect(x), collect(y), amp, phase, dir, λ;
                randomize_axes = true, rng = MersenneTwister(5))
        end,
    ]

    @testset "UniformDiscSource is a CollimatedSource" begin
        @test UniformDiscSource(pos, dir, 0.01, λ; num_rays = 10) isa CollimatedSource
    end

    @testset "Round trip: $name" for (name, src) in sources
        loaded, _ = roundtrip(src)
        test_identical(src, loaded)
    end

    @testset "Traced state is not written, loaded source can be traced: $name" for (name, src) in sources
        system = test_system()
        solve_system!(system, src)
        traced = traced_state(src)
        @test any(length(t) > 1 for t in traced)
        loaded, _ = roundtrip(src)
        test_identical(src, loaded)
        solve_system!(test_system(), loaded)
        @test traced_state(loaded) == traced
    end

    @testset "Large groups are stored as assets" begin
        src = UniformDiscSource(pos, dir, 0.01, λ; num_rays = 10_000)
        loaded, path = roundtrip(src)
        test_identical(src, loaded)
        setup, _ = BMO._read_archive(path)
        packed = only(setup["sources"])["beams"]
        @test packed["count"] == 10_000
        for key in ("pos", "dir", "λ", "n")
            @test haskey(packed[key], "asset")
            @test !haskey(packed[key], "data")
        end

        bg = GaussianBeamletDecomposition(pos, dir, λ, 1e-3; n_grid = 40, randomize_axes = true, rng = MersenneTwister(6))
        loaded, path = roundtrip(bg)
        test_identical(bg, loaded)
        setup, _ = BMO._read_archive(path)
        entry = only(setup["sources"])
        @test haskey(entry["aux"]["pos"], "asset")
        @test haskey(entry["chief"]["E0"], "asset")
    end

    @testset "Shared sources stay shared" begin
        src = CollimatedSource(pos, dir, 0.01, λ; num_rings = 2, num_rays = 40)
        path = tmpfile()
        # Empty systems: object storage is tested elsewhere
        BMO.save_setup(path, System(BMO.AbstractObject[]) => src, System(BMO.AbstractObject[]) => src)
        loaded = BMO.load_setup(path)
        @test length(loaded.sources) == 1
        @test last(loaded.pairs[1]) === last(loaded.pairs[2]) === only(loaded.sources)
    end

    @testset "Only Float64 sources" begin
        @test_throws ArgumentError BMO.save_setup(tmpfile(); sources = [Beam(Float32[0, 0, 0], Float32[0, 1, 0], 1f-6)])
        @test_throws ArgumentError BMO.save_setup(tmpfile();
            sources = [CollimatedSource([Beam(Float32[0, 0, 0], Float32[0, 1, 0], 1f-6)], 1f-2)])
    end
end

end # module
