module TestDetector

using BeamletOptics
using Test

const BMO = BeamletOptics

const mm = 1e-3

@testset "Detector" begin

@testset "Testing spot diagram" begin
    # Init system and one ring beam
    pd = Detector(10e-3)
    translate3d!(pd, [1mm, 50mm, -2mm])
    zrotate3d!(pd, deg2rad(45))
    xrotate3d!(pd, deg2rad(-45))    
    system = System([pd])    
    aperture = 2mm
    num_rings = 2
    num_rays = 101    
    cs = CollimatedSource([0,0,0], [0,1,0], aperture, 1e-6; num_rings, num_rays)
    #
    solve_system!(system, cs)
    sd = spot_diagram(pd)
    # test correct num of hits
    @test length(BMO.hits(pd)) == num_rays
    # test correct ellipse form
    xs = getindex.(sd, 1)
    zs = getindex.(sd, 2)    
    @test minimum(xs)*1000 ≈ -2.8279    atol = 1e-3
    @test maximum(xs)*1000 ≈ 0          atol = 1e-3
    @test minimum(zs)*1000 ≈ 0.09707    atol = 1e-3
    @test maximum(zs)*1000 ≈ 3.55978    atol = 1e-3
    # test empty! fct.
    empty!(pd)
    @test isnothing(BMO.hits(pd))
end

@testset "Testing continued tracing (stop = false)" begin
    # Regression test: interact3d used to resume the next segment at the
    # incoming ray's origin instead of the actual hit point, which made
    # `stop = false` re-hit the same detector forever (bounded only by r_max).
    d1 = Detector(20mm, false)
    d2 = Detector(20mm, true)
    translate_to3d!(d1, [0.0, 0.10, 0.0])
    translate_to3d!(d2, [0.0, 0.20, 0.0])
    system = System([d1, d2])

    beam = Beam(Ray([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]))
    solve_system!(system, beam)

    @test length(BMO.hits(d1)) == 1
    @test length(BMO.hits(d2)) == 1
    @test length(rays(beam)) == 2

    hit = BMO.hits(d1)[1]
    @test BMO.hit_point(hit) ≈ [0.0, 0.10, 0.0]
    # the continued segment must start at the hit point, not the incoming ray's origin
    @test position(rays(beam)[2]) ≈ [0.0, 0.10, 0.0]
end

@testset "Testing point spread function" begin
    # parameters for an almost thin-lens
    l = 1mm
    R1 = 100mm
    R2 = Inf
    d = 25.4mm
    n = 1.5
    λ = 1e-6    
    D = 15mm
    num_rays = 1000
    
    # plane wave source
    cs = UniformDiscSource([0, -10e-3, 0], [0, 1, 0], D, λ; num_rays)
    
    # test lens
    lens = SphericalLens(R1, R2, l, d, x -> n)
    
    # PSF detector
    x_shift = y_shift = -2mm
    psfd = Detector(10e-3)
    translate3d!(psfd, [x_shift, 200e-3 + 0.13e-3, y_shift])

    @testset "Airy-disc test" begin
        # build system and solve it
        sys = System([lens, psfd])
        solve_system!(sys, cs)

        # get PSF
        x, y, I_num = intensity(psfd; n=500, crop_factor=5, center=MinMax())

        # walk from the peak to the first local minimum through the centre column
        ix_ctr, jx_ctr = Tuple(argmax(I_num))
        col = I_num[:, jx_ctr]
        i_min = ix_ctr
        while i_min < length(col) && col[i_min + 1] < col[i_min]
            i_min += 1
        end

        # compare relative to the peak, since the absolute offset (x_shift) dwarfs the Airy radius
        airy_radius = 1.22*λ*200e-3/D
        @test x[i_min] - x[ix_ctr] ≈ airy_radius rtol=2e-2
    end

    @testset "Detector reset" begin
        # Test if data in psfd
        @test length(BMO.hits(psfd)) == num_rays
        # Empty psfd and test
        empty!(psfd)
        @test isnothing(BMO.hits(psfd))
    end
end

@testset "Testing Gaussian beamlet interference" begin
    @testset "Pre-Beamsplitter tests with separate beams" begin
        # Gauss beam parameters (selected for ring fringes)
        w0 = 0.01e-3
        λ = 1000e-9
        M2 = 1
        P0 = 1e-3
        I0 = 2 * P0 / (π * w0^2)
        E0 = BMO.electric_field(I0)
        zR = BMO.rayleigh_range(λ, w0, M2)
        # Detector parameters
        z = 0.1     # distance to detector
        l = 1e-2    # detector size
        n = 1000    # detector grid resolution
        # Lens parameters
        R1 = R2 = d = 0.01
        nl = 1.5
        f = BMO.lensmakers_eq(R1, -R2, nl)
        # Raytracing system (for all tests)
        pd_l = Detector(l) # n
        pd_s = Detector(l / 10) # n ÷ 10
        ln = ThinLens(R1, R2, d, nl)
        translate3d!(pd_l, [0, z, 0])
        translate3d!(pd_s, [0, z, 0])
        translate3d!(ln, [0, z - f - thickness(ln) / 2, 0])

        @testset "Testing fringe pattern" begin
            system = System(pd_l)
            Δz = 5e-3   # arm length difference
            # Analytic solution
            xs = ys = LinRange(-l / 2, l / 2, n)
            screen = zeros(ComplexF64, length(xs), length(ys))
            for (j, y) in enumerate(ys)
                for (i, x) in enumerate(xs)
                    r = sqrt(x^2 + y^2)
                    screen[i, j] += BMO.electric_field(r, z, E0, w0, λ, M2)
                    screen[i, j] += BMO.electric_field(r, z + Δz, E0, w0, λ, M2)
                end
            end
            # Numerical solution
            empty!(pd_l)
            g_1 = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0],
                λ,
                w0,
                M2 = M2,
                P0 = P0)
            g_2 = GaussianBeamlet([0.0, -Δz, 0], [0.0, 1, 0],
                λ,
                w0,
                M2 = M2,
                P0 = P0)
            solve_system!(system, g_1)
            solve_system!(system, g_2)

            # Compare solutions
            I_analytical = intensity.(screen)
            ~, ~, I_numerical = intensity(pd_l; n, x_min=-l/2, x_max=l/2, z_min=-l/2, z_max=l/2)
            Pt = optical_power(pd_l; n, x_min=-l/2, x_max=l/2, z_min=-l/2, z_max=l/2)
            @test all(isapprox.(I_analytical, I_numerical, atol = 2e-1))
            @test isapprox(Pt, 2 * P0, atol = 3e-5)
        end

        @testset "Testing λ phase shift" begin
            system = System([pd_s, ln])
            # Numerical solution
            Δz = LinRange(0, λ, 50)
            Pt_numerical = zeros(length(Δz))
            for (i, z_i) in enumerate(Δz)
                empty!(pd_s)
                g_1 = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0],
                    λ,
                    w0,
                    M2 = M2,
                    P0 = P0)
                g_2 = GaussianBeamlet([0.0, z_i, 0], [0.0, 1, 0],
                    λ,
                    w0,
                    M2 = M2,
                    P0 = P0)
                solve_system!(system, g_1)
                solve_system!(system, g_2)
                Pt_numerical[i] = BMO.optical_power(pd_s)
                # Test length/opl function
                @test length(g_1) ≈ z
                @test length(g_1) ≈ BMO.optical_path_length(g_1) - thickness(ln) * (nl - 1)
                @test length(g_2) ≈ length(g_1) - z_i
            end
            # Analytical solution (cosine over Δz), ref. power is 4*P0 since beamsplitter is missing
            Pt_analytical = 4 * P0 * [(cos(2π * z / (maximum(Δz))) + 1) / 2 for z in Δz]

            # Compare detectors (this also tests correct behavior when focussing the beam)
            @test all(isapprox.(Pt_numerical, Pt_analytical, atol = 1e-4))
        end
    end
end

@testset "initialize!" begin
    # two detectors, one nested two levels deep in groups, plus a mirror that stores nothing
    pd1 = Detector(10mm, false)
    pd2 = Detector(10mm)
    translate3d!(pd1, [0, 50mm, 0])
    translate3d!(pd2, [0, 60mm, 0])
    mirror = RoundPlanoMirror(10mm, 5mm)
    translate3d!(mirror, [0, 200mm, 0])
    inner = ObjectGroup([pd2])
    outer = ObjectGroup([mirror, inner])
    system = System([pd1, outer])
    make_source() = CollimatedSource([0, 0, 0], [0, 1, 0], 2mm, 1e-6; num_rings = 2, num_rays = 41)
    solve_system!(system, make_source(); progress = false)
    @test BMO.hit_count(pd1) == 41
    @test BMO.hit_count(pd2) == 41
    # accumulation is intended: a second source solved afterwards superposes
    solve_system!(system, make_source(); progress = false)
    @test BMO.hit_count(pd1) == 82
    # an object that stores nothing is left alone
    @test isnothing(initialize!(mirror))
    @test isnothing(initialize!(pd2))
    @test isnothing(BMO.hits(pd2))
    solve_system!(system, make_source(); progress = false)
    # nested groups are reached from groups and from systems
    @test isnothing(initialize!(inner))
    @test isnothing(BMO.hits(pd2))
    @test BMO.hit_count(pd1) == 123
    solve_system!(system, make_source(); progress = false)
    @test isnothing(initialize!(system))
    @test isnothing(BMO.hits(pd1))
    @test isnothing(BMO.hits(pd2))
    static = StaticSystem([pd1, outer])
    solve_system!(static, make_source(); progress = false)
    @test BMO.hit_count(pd2) == 41
    @test isnothing(initialize!(static))
    @test isnothing(BMO.hits(pd1))
    @test isnothing(BMO.hits(pd2))
    @test isnothing(initialize!(System()))
end

@testset "solve_system! with initialize" begin
    @testset "Beam" begin
        pd = Detector(10mm)
        translate3d!(pd, [0, 50mm, 0])
        system = System([pd])
        beam = Beam([0.0, 0, 0], [0.0, 1, 0], 1e-6)
        solve_system!(system, beam)
        solve_system!(system, beam)
        # by default the hits of both solves add up
        @test BMO.hit_count(pd) == 2
        solve_system!(system, beam; initialize = true)
        @test BMO.hit_count(pd) == 1
        solve_system!(system, beam; initialize = false)
        @test BMO.hit_count(pd) == 2
    end

    @testset "Beam group" begin
        pd = Detector(20mm)
        translate3d!(pd, [0, 50mm, 0])
        system = System([pd])
        cs = CollimatedSource([0, 0, 0], [0, 1, 0], 10mm, 1e-6; num_rings = 20, num_rays = 2000)
        solve_system!(system, cs; progress = false)
        # every beam of the group contributes exactly one hit, none is lost to a race
        @test BMO.hit_count(pd) == length(cs)
        # solving the group again without initializing accumulates
        solve_system!(system, cs; progress = false)
        @test BMO.hit_count(pd) == 2 * length(cs)
        # the system is initialized once, not for every beam of the group
        solve_system!(system, cs; progress = false, initialize = true)
        @test BMO.hit_count(pd) == length(cs)
        solve_system!(system, cs; progress = false, initialize = true)
        @test BMO.hit_count(pd) == length(cs)
    end

    @testset "Scan loop equals a newly built setup" begin
        function setup()
            m = RoundPlanoMirror(25.4mm, 5mm)
            pd = Detector(10mm)
            translate3d!(m, [0, 20mm, 0])
            zrotate3d!(m, deg2rad(45))
            translate3d!(pd, [20mm, 20mm, 0])
            zrotate3d!(pd, deg2rad(90))
            return m, pd, System([m, pd])
        end
        beam() = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 632.8e-9, 1e-3)
        m, pd, system = setup()
        b = beam()
        solve_system!(system, b)
        n_first = BMO.hit_count(pd)
        @test n_first > 0
        m2, pd2, system2 = setup()
        zrotate3d!(m2, 1e-3)
        solve_system!(system2, beam())
        # change the setup and solve again, via the keyword ...
        zrotate3d!(m, 1e-3)
        solve_system!(system, b; initialize = true)
        @test BMO.hit_count(pd) == BMO.hit_count(pd2)
        @test optical_power(pd) ≈ optical_power(pd2)
        # ... and via the function
        initialize!(system)
        solve_system!(system, b)
        @test BMO.hit_count(pd) == BMO.hit_count(pd2)
        @test optical_power(pd) ≈ optical_power(pd2)
    end
end

@testset "Hit buffers of a beam group solve" begin
    # two detectors side by side, the beams of a group hit both
    pd1, pd2 = Detector(20mm), Detector(20mm)
    translate3d!(pd1, [-10mm, 50mm, 0])
    translate3d!(pd2, [10mm, 50mm, 0])
    system = System([pd1, pd2])
    cs = CollimatedSource([0, 0, 0], [0, 1, 0], 30mm, 1e-6; num_rings = 10, num_rays = 500)
    solve_system!(system, cs; progress = false)
    n1, n2 = BMO.hit_count(pd1), BMO.hit_count(pd2)
    @test n1 > 0 && n2 > 0
    # the group stores the same hits as its beams solved one by one
    pts(pd) = sort(map(h -> Tuple(BMO.hit_point(h)), BMO.hits(pd)))
    group = (pts(pd1), pts(pd2))
    initialize!(system)
    foreach(b -> solve_system!(system, b), BMO.beams(cs))
    @test (pts(pd1), pts(pd2)) == group
    # within `with_hit_buffers` the hits are collected, and stored also if the code throws
    initialize!(system)
    b = first(BMO.beams(cs))
    @test_throws ErrorException BMO.with_hit_buffers() do
        solve_system!(system, b)
        @test BMO.hit_count(pd1) + BMO.hit_count(pd2) == 0
        error("cancelled")
    end
    @test BMO.hit_count(pd1) + BMO.hit_count(pd2) == 1
end

@testset "Default window holds the beam power" begin
    # https://github.com/JuliaPhysics/BeamletOptics.jl/issues/127
    P0 = 1e-3
    for beam in (GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 1e-6, 1mm; P0),
        AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 1e-6, 1mm; P0))
        pd = Detector(20mm)
        translate3d!(pd, [0, 100mm, 0])
        solve_system!(System([pd]), beam)
        @test optical_power(pd) ≈ P0 rtol = 1e-4
        x_min, x_max, z_min, z_max = BMO.calc_local_lims(pd)
        w = first(BMO.hits(pd)).w_max
        @test x_max - x_min ≈ 6w
        @test z_max - z_min ≈ 6w
    end
end

end # TESTSET

end # MODULE