module TestStorageE2E

using BeamletOptics
using Test

const BMO = BeamletOptics

const mm = 1e-3

tmpfile() = joinpath(mktempdir(), "setup.bmo")

N_DG1(λ) = 1.62286 + 1e3 * (λ - 486e-9)

# Hit order may depend on threading, compare sorted spots
sorted_spots(pd) = sort([(p[1], p[2]) for p in eachrow(reduce(hcat, spot_diagram(pd))')])

function michelson()
    l_0 = 0.1
    m1 = SquarePlanoMirror2D(BMO.inch)
    m2 = SquarePlanoMirror2D(BMO.inch)
    bs = ThinBeamsplitter(BMO.inch, reflectance = 0.5)
    pd = Detector(BMO.inch / 5)
    translate3d!(m1, [l_0, 0, 0])
    translate3d!(m2, [0, l_0 + 1e-7, 0])
    translate3d!(pd, [-l_0, 0, 0])
    zrotate3d!(bs, deg2rad(45))
    zrotate3d!(m1, deg2rad(90))
    zrotate3d!(pd, deg2rad(90))
    beam = GaussianBeamlet([0, -l_0, 0], [0, 1.0, 0], 635e-9, 1e-4, P0 = 5mm)
    return System([m1, m2, bs, pd]), beam, pd
end

function double_gauss()
    # A user-defined dispersion function must be registered by name to be storable
    register_material!("DG-glass-1", N_DG1)
    l1 = SphericalLens(48.88mm, 182.96mm, 8.89mm, 52.3mm, N_DG1)
    l23 = SphericalDoubletLens(36.92mm, Inf, 23.06mm, 15.11mm, 2.31mm,
        45.11mm, ConstantRefractiveIndex(1.58565), ConstantRefractiveIndex(1.67764))
    l45 = SphericalDoubletLens(-23.91mm, Inf, -36.92mm, 1.92mm, 7.77mm,
        40.01mm, ConstantRefractiveIndex(1.57046), ConstantRefractiveIndex(1.64128))
    l6 = SphericalLens(1063.24mm, -48.88mm, 6.73mm, 45.11mm, 1.62286)
    l_23 = thickness(l1) + 0.38mm
    l_45 = l_23 + thickness(l23) + 9.14mm + 13.36mm
    l_6 = l_45 + thickness(l45) + 0.38mm
    f_z = l_6 + thickness(l6.shape) + 58.21mm
    translate3d!(l23, [0, l_23, 0])
    translate3d!(l45, [0, l_45, 0])
    translate3d!(l6, [0, l_6, 0])
    lenses = ObjectGroup([l1, l23, l45, l6])
    detector = Detector(5mm)
    translate3d!(detector, [0, f_z, 0])
    setup = ObjectGroup([lenses, detector])
    translate3d!(setup, [0.05, 0.05, 0.05])
    xrotate3d!(setup, deg2rad(60))
    zrotate3d!(setup, deg2rad(45))
    dir = orientation(lenses)[:, 2]
    pos = position(l1) - 0.05 * dir
    source = CollimatedSource(pos, dir, 0.04, 486.0e-9, num_rays = 500, num_rings = 10)
    return System([setup]), source, detector
end

@testset "Storage end-to-end" begin
    @testset "Michelson interferometer" begin
        system, beam, pd = michelson()
        path = tmpfile()
        save_setup(path, system => beam; names = IdDict{Any, String}(pd => "pd"))
        setup = load_setup(path)
        lsystem, lbeam = only(setup.pairs)
        lpd = setup.names["pd"]
        @test lpd isa Detector
        @test lpd in lsystem.objects

        solve_system!(system, beam)
        solve_system!(lsystem, lbeam)
        kw = (; n = 100, x_min = -BMO.inch / 10, x_max = BMO.inch / 10, z_min = -BMO.inch / 10, z_max = BMO.inch / 10)
        _, _, I = intensity(pd; kw...)
        _, _, lI = intensity(lpd; kw...)
        @test maximum(abs.(I .- lI)) <= 1e-10 * maximum(abs.(I))
        @test isapprox(optical_power(pd; kw...), optical_power(lpd; kw...); rtol = 1e-10)
    end

    @testset "Double Gauss lens in nested groups" begin
        system, source, detector = double_gauss()
        path = tmpfile()
        save_setup(path, system => source; names = IdDict{Any, String}(detector => "detector"))
        setup = load_setup(path)
        lsystem, lsource = only(setup.pairs)
        ldetector = setup.names["detector"]
        # The group tree survives
        @test only(lsystem.objects) isa ObjectGroup
        @test ldetector === only(lsystem.objects).objects[2]

        solve_system!(system, source)
        solve_system!(lsystem, lsource)
        spots, lspots = sorted_spots(detector), sorted_spots(ldetector)
        @test length(spots) == length(lspots) > 0
        @test maximum(maximum(abs.(a .- b)) for (a, b) in zip(spots, lspots)) <= 1e-12
    end

    @testset "Several systems and sources" begin
        mi, mi_beam, _ = michelson()
        dg, dg_source, _ = double_gauss()
        path = tmpfile()
        save_setup(path, mi => mi_beam, dg => dg_source, dg => mi_beam;
            metadata = Dict("gui" => Dict("auto_trace" => true)))
        setup = load_setup(path)
        @test length(setup.pairs) == 3
        @test setup.pairs[1][2] === setup.pairs[3][2]
        @test setup.pairs[2][1] === setup.pairs[3][1]
        @test setup.metadata["gui"]["auto_trace"] == true
    end

    @testset "Fixture of the current series" begin
        # Written by fixtures/make_fixture.jl. When the breaking series changes, the file must be regenerated.
        fixture = joinpath(@__DIR__, "fixtures", "michelson.bmo")
        setup = load_setup(fixture)
        v = pkgversion(BeamletOptics)
        @test BMO._series(setup.version) == BMO._series(v)
        lsystem, lbeam = only(setup.pairs)
        solve_system!(lsystem, lbeam)
        @test optical_power(setup.names["pd"]) > 0
        other = v.major == 0 ? VersionNumber(0, v.minor + 1, 0) : VersionNumber(v.major + 1, 0, 0)
        @test_throws r"load it with a BeamletOptics .* release" load_setup(fixture; _loader = other)
    end
end

end # module
