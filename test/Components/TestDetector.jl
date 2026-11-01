module TestDetector

using BeamletOptics
using Test
using LinearAlgebra: norm

const BMO = BeamletOptics

const mm = 1e-3

# a hit type that implements none of the readout functions
struct DummyHit <: BMO.AbstractDetectorHit end

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

@testset "Ray hits: points and directions in detector coordinates" begin
    # https://github.com/JuliaPhysics/BeamletOptics.jl/issues/191: rays from the origin onto a
    # detector at distance L, then the whole setup moved, which changes no local value
    L = 50mm
    dirs = [[sin(θ), cos(θ), 0.0] for θ in deg2rad.([-20, -5, 0, 10, 30])]
    push!(dirs, [0.0, cos(0.2), sin(0.2)])
    # crossings at distance L along the normal
    pts_exp = [[d[1], d[3]] * L / d[2] for d in dirs]
    axis = [1.0, 2, 3] / norm([1.0, 2, 3])
    for moved in (false, true)
        pd = Detector(100mm)
        translate3d!(pd, [0, L, 0])
        beams = [Beam([0.0, 0, 0], d) for d in dirs]
        if moved
            for x in (pd, beams...)
                rotate3d!(x, axis, deg2rad(40), [0.0, 0, 0])
                translate3d!(x, [10mm, -5mm, 3mm])
            end
        end
        system = System([pd])
        foreach(b -> solve_system!(system, b), beams)
        pts, ds = ray_hits(pd)
        @test pts isa Vector{BMO.Point2{Float64}} && ds isa Vector{BMO.Point3{Float64}}
        @test maximum(norm.(ds .- dirs)) < 1e-12
        @test maximum(norm.(pts .- pts_exp)) < 1e-12
        # angle to the normal, the same quantity in both placements
        θ(d) = atan(hypot(d[1], d[3]), abs(d[2]))
        @test θ.(ds) ≈ [deg2rad.([20, 5, 0, 10, 30]); 0.2] atol = 1e-12
        @test pts == spot_diagram(pd)
    end
    # A hit on the back side keeps its side in d[2], the angle is that to the normal
    pd = Detector(100mm)
    translate3d!(pd, [0, L, 0])
    d_back = [sin(0.1), -cos(0.1), 0.0]
    solve_system!(System([pd]), Beam([0.0, 2L, 0], d_back))
    pts, ds = ray_hits(pd)
    @test norm(only(ds) - d_back) < 1e-12
    @test norm(only(pts) - [tan(0.1) * L, 0]) < 1e-12
    @test atan(hypot(only(ds)[1], only(ds)[3]), abs(only(ds)[2])) ≈ 0.1 atol = 1e-12
    # Beamlets: the chief ray, not the waist contour of spot_diagram
    d = [sin(0.1), cos(0.1), 0.0]
    for beam in (GaussianBeamlet([0.0, 0, 0], d, 1e-6, 1mm),
        AstigmaticGaussianBeamlet([0.0, 0, 0], d, 1e-6, 1mm))
        pd = Detector(100mm)
        translate3d!(pd, [0, L, 0])
        solve_system!(System([pd]), beam)
        pts, ds = ray_hits(pd)
        @test length(pts) == length(ds) == 1
        @test norm(only(ds) - d) < 1e-12
        @test norm(only(pts) - [d[1], d[3]] * L / d[2]) < 1e-12
    end
    # A detector without hits
    @test ray_hits(Detector(10mm)) == (BMO.Point2{Float64}[], BMO.Point3{Float64}[])
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
        x, y, I_num = intensity(psfd, RegularGrid(num_points = 500, crop_factor=5, center=MinMax()))

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
            ~, ~, I_numerical = intensity(pd_l, RegularGrid(; num_points = n, x_min=-l/2, x_max=l/2, z_min=-l/2, z_max=l/2))
            Pt = optical_power(pd_l, RegularGrid(; num_points = n, x_min=-l/2, x_max=l/2, z_min=-l/2, z_max=l/2))
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
                # The 3 mm wide beam behind the ball-like lens leaves the model (Lagrange
                # invariant at 3 % of its start value), the phase and power are still valid
                solve_system!(system, g_1; check_invariant = false)
                solve_system!(system, g_2; check_invariant = false)
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
    # no hits before the first solve
    @test length(pd1) == 0
    solve_system!(system, make_source(); progress = false)
    @test length(pd1) == 41
    @test length(pd2) == 41
    # accumulation is intended: a second source solved afterwards superposes
    solve_system!(system, make_source(); progress = false)
    @test length(pd1) == 82
    # an object that stores nothing is left alone
    @test isnothing(initialize!(mirror))
    @test isnothing(initialize!(pd2))
    @test isnothing(BMO.hits(pd2))
    solve_system!(system, make_source(); progress = false)
    # nested groups are reached from groups and from systems
    @test isnothing(initialize!(inner))
    @test isnothing(BMO.hits(pd2))
    @test length(pd1) == 123
    solve_system!(system, make_source(); progress = false)
    @test isnothing(initialize!(system))
    @test isnothing(BMO.hits(pd1))
    @test isnothing(BMO.hits(pd2))
    static = StaticSystem([pd1, outer])
    solve_system!(static, make_source(); progress = false)
    @test length(pd2) == 41
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
        @test length(pd) == 2
        solve_system!(system, beam; initialize = true)
        @test length(pd) == 1
        solve_system!(system, beam; initialize = false)
        @test length(pd) == 2
    end

    @testset "Beam group" begin
        pd = Detector(20mm)
        translate3d!(pd, [0, 50mm, 0])
        system = System([pd])
        cs = CollimatedSource([0, 0, 0], [0, 1, 0], 10mm, 1e-6; num_rings = 20, num_rays = 2000)
        solve_system!(system, cs; progress = false)
        # every beam of the group contributes exactly one hit, none is lost to a race
        @test length(pd) == length(cs)
        # solving the group again without initializing accumulates
        solve_system!(system, cs; progress = false)
        @test length(pd) == 2 * length(cs)
        # the system is initialized once, not for every beam of the group
        solve_system!(system, cs; progress = false, initialize = true)
        @test length(pd) == length(cs)
        solve_system!(system, cs; progress = false, initialize = true)
        @test length(pd) == length(cs)
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
        n_first = length(pd)
        @test n_first > 0
        m2, pd2, system2 = setup()
        zrotate3d!(m2, 1e-3)
        solve_system!(system2, beam())
        # change the setup and solve again, via the keyword ...
        zrotate3d!(m, 1e-3)
        solve_system!(system, b; initialize = true)
        @test length(pd) == length(pd2)
        @test optical_power(pd) ≈ optical_power(pd2)
        # ... and via the function
        initialize!(system)
        solve_system!(system, b)
        @test length(pd) == length(pd2)
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
    n1, n2 = length(pd1), length(pd2)
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
        @test length(pd1) + length(pd2) == 0
        error("cancelled")
    end
    @test length(pd1) + length(pd2) == 1
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

@testset "Hits are values" begin
    # https://github.com/JuliaPhysics/BeamletOptics.jl/issues/42: a hit holds no reference to the
    # beam that made it, hence a vector of hits is one block of memory
    for H in (BMO.RayHit, BMO.PolarizedRayHit, BMO.GaussianBeamletHit, BMO.AstigmaticGaussianBeamletHit)
        @test isbitstype(H{Float64})
    end
end

@testset "A new solve does not change the stored hits" begin
    # https://github.com/JuliaPhysics/BeamletOptics.jl/issues/42
    @testset "Beam rotated between solves" begin
        pd = Detector(50mm)
        translate3d!(pd, [0, 100mm, 0])
        system = System([pd])
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        for _ in 1:5
            solve_system!(system, beam)
            zrotate3d!(beam, deg2rad(1))
        end
        # the hit of the i-th solve lies 100 mm · tan((i − 1)°) to the side
        xs = sort(abs.(first.(spot_diagram(pd))))
        @test xs ≈ 100mm .* tan.(deg2rad.(0:4)) atol = 1e-12
    end
    E0 = [1.0, 0, 0]
    beams = (
        "Beam" => () -> Beam([0.0, 0, 0], [0.0, 1, 0], 1e-6),
        "Beam of PolarizedRays" => () -> Beam([0.0, 0, 0], [0.0, 1, 0], 1e-6, E0),
        "GaussianBeamlet" => () -> GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 1e-6, 1mm),
        "AstigmaticGaussianBeamlet" => () -> AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 1e-6, 1mm))
    # the hit is on the first ray, on the second ray behind a mirror, or on a child beam
    for (label, beam) in beams, path in (:direct, :mirror, :splitter)
        @testset "$label, $path" begin
            pd = Detector(20mm)
            system = System([pd])
            if path === :mirror
                m = RoundPlanoMirror(30mm, 5mm)
                zrotate3d!(m, deg2rad(45))
                translate3d!(m, [0, 50mm, 0])
                zrotate3d!(pd, deg2rad(90))
                translate3d!(pd, [100mm, 50mm, 0])
                push!(system, m)
            else
                translate3d!(pd, [0, 150mm, 0])
            end
            if path === :splitter
                bs = ThinBeamsplitter(30mm)
                zrotate3d!(bs, deg2rad(45))
                translate3d!(bs, [0, 50mm, 0])
                push!(system, bs)
            end
            b = beam()
            solve_system!(system, b)
            before = only(BMO.hits(pd))
            # move the source and the mirror or splitter, the beam still hits the detector
            translate3d!(b, [0, 0, 0.5mm])
            path === :mirror && translate3d!(m, [5mm, 0, 0])
            path === :splitter && translate3d!(bs, [0, 10mm, 0])
            solve_system!(system, b)
            @test length(pd) == 2
            @test first(BMO.hits(pd)) === before
            @test BMO.hit_point(last(BMO.hits(pd)))[3] ≈ 0.5mm
        end
    end
    @testset "Field of a stored GaussianBeamletHit" begin
        # free space from the waist: E0 w0 / w(L) · exp(i (k L − atan(L / z_R))) on the axis
        λ, w0, L = 1e-6, 0.5mm, 300mm
        pd = Detector(20mm)
        translate3d!(pd, [0, L, 0])
        system = System([pd])
        gb = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0)
        solve_system!(system, gb)
        hit = only(BMO.hits(pd))
        translate3d!(gb, [0, 0, 1mm])
        solve_system!(system, gb)
        z_R = π * w0^2 / λ
        w = w0 * sqrt(1 + (L / z_R)^2)
        E = electric_field(gb) * w0 / w * cis(2π / λ * L - atan(L / z_R))
        # the waist and divergence rays give the Gouy phase to about 1e-7 rad
        @test electric_field(first(BMO.hits(pd)), BMO.Point3(0.0, L, 0.0)) ≈ E rtol = 1e-6
        @test first(BMO.hits(pd)) === hit
    end

    @testset "Field of a beamlet hit as from the former scalar field" begin
        # the former scalar field of a beamlet hit, from the stored fields of the hit
        function scalar_field(hit::BMO.GaussianBeamletHit, p)
            v = p - hit.p0
            l1 = BMO._pseudo_dot(v, hit.dir)
            r = norm(v - l1 * hit.dir)
            w, curv, ψ, _ = BMO.gauss_parameters(hit.segment, l1)
            ψ += hit.ψ0
            k0 = BMO.wavenumber(hit.segment.λ)
            n = hit.segment.n
            opl = hit.l0 + l1 + hit.Δl0 + (n - 1) * l1
            return hit.E_amp / (sqrt(n) * w) * exp(-r^2 / w^2) *
                   cis(k0 * opl + ψ + n * k0 * r^2 * curv / 2)
        end
        # an astigmatic hit: the phasor along the polarization of the chief ray
        scalar_field(hit::BMO.AstigmaticGaussianBeamletHit, p) = BMO.astigmatic_phasor(hit, p) .* hit.E0
        relative_error(E, E_ref) = norm(E - E_ref) / norm(E_ref)

        λ, w0, L = 1e-6, 0.5mm, 300mm
        p = BMO.Point3(0.2mm, L, -0.1mm)
        for beam in (GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0),
                     AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0, 2w0;
                         support = [1.0, 0, 0]))
            pd = Detector(20mm)
            translate3d!(pd, [0, L, 0])
            solve_system!(System([pd]), beam)
            hit = only(BMO.hits(pd))
            @test relative_error(electric_field(hit, p), scalar_field(hit, p)) < 1e-12
        end
    end
end

@testset "Moving a detector empties it" begin
    pivot = [1mm, 2mm, 3mm]
    R = BMO.rotate3d([0.0, 0, 1], 0.1)
    moves = (
        "translate3d!" => pd -> translate3d!(pd, [0, 1mm, 0]),
        "translate_to3d!" => pd -> translate_to3d!(pd, [0, 120mm, 0]),
        "reset_translation3d!" => pd -> reset_translation3d!(pd),
        "rotate3d!(R)" => pd -> rotate3d!(pd, R),
        "rotate3d!(R, pivot)" => pd -> rotate3d!(pd, R, pivot),
        "rotate3d!(axis, θ)" => pd -> rotate3d!(pd, [0.0, 0, 1], 0.1),
        "rotate3d!(axis, θ, pivot)" => pd -> rotate3d!(pd, [0.0, 0, 1], 0.1, pivot),
        "xrotate3d!" => pd -> xrotate3d!(pd, 0.1),
        "yrotate3d!" => pd -> yrotate3d!(pd, 0.1),
        "zrotate3d!" => pd -> zrotate3d!(pd, 0.1),
        "align3d!" => pd -> align3d!(pd, [0.1, 1, 0]),
        "reset_rotation3d!" => pd -> reset_rotation3d!(pd),
        "translate3d! of a group" => pd -> translate3d!(ObjectGroup([pd]), [0, 1mm, 0]),
        "zrotate3d! of a group" => pd -> zrotate3d!(ObjectGroup([pd]), 0.1))
    for (label, move!) in moves
        @testset "$label" begin
            pd = Detector(50mm)
            translate3d!(pd, [0, 100mm, 0])
            solve_system!(System([pd]), Beam([0.0, 0, 0], [0.0, 1, 0]))
            @test length(pd) == 1
            move!(pd)
            @test isnothing(BMO.hits(pd))
        end
    end
    @testset "Moving other objects keeps the hits" begin
        pd = Detector(50mm)
        translate3d!(pd, [0, 100mm, 0])
        lens = ThinLens(80mm, 80mm, 4mm, 1.5)
        translate3d!(lens, [0, 50mm, 0])
        beam = Beam([0.0, 0, 0], [0.0, 1, 0])
        solve_system!(System([lens, pd]), beam)
        translate3d!(lens, [1mm, 0, 0])
        translate3d!(beam, [0, 0, 1mm])
        @test length(pd) == 1
    end
end

@testset "Flux through a tilted detector" begin
    P0, λ, w0 = 1e-3, 1e-6, 1e-3
    function tilted(θ)
        pd = Detector(50mm)
        zrotate3d!(pd, θ)
        translate3d!(pd, [0, 100mm, 0])
        return pd
    end
    gb(pos = [0.0, 0, 0], dir = [0.0, 1, 0]) = GaussianBeamlet(pos, dir, λ, w0; P0)
    agb() = AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; P0)

    @testset "$name holds the beam power at $(θd)°" for (name, beam) in (("GaussianBeamlet", gb),
            ("AstigmaticGaussianBeamlet", agb)), θd in (0, 30, 60)
        pd = tilted(deg2rad(θd))
        solve_system!(System([pd]), beam())
        @test optical_power(pd, RegularGrid(num_points = 300); progress = false) ≈ P0 rtol = 1e-4
    end

    @testset "The field is not projected onto the detector" begin
        pd = tilted(deg2rad(60))
        beam = gb()
        solve_system!(System([pd]), beam)
        E = electric_field(pd, RectilinearGrid([0.0], [0.0]); progress = false)[3][1]
        @test abs(E) ≈ abs(BMO.electric_field(beam, 0.0, 100mm)) rtol = 1e-9
    end

    @testset "Irradiance of a ray falls with cos θ" begin
        θ = deg2rad(60)
        d, p = [sin(θ), cos(θ), 0], [-100mm * tan(θ), 0, 0]
        function irradiance(ray)
            pd = tilted(0)
            solve_system!(System([pd]), Beam(ray))
            return intensity(pd, RectilinearGrid([0.0], [0.0]); progress = false)[3][1]
        end
        # a ray has the amplitude 1, a polarized ray that of its field vector
        @test irradiance(Ray(p, d, λ)) ≈ cos(θ) / (2 * BMO.Z_vacuum)
        @test irradiance(PolarizedRay(p, d, λ, [0, 0, 2.0])) ≈ 4cos(θ) / (2 * BMO.Z_vacuum)
    end

    @testset "The phase of a ray grows with the optical distance" begin
        # in a medium, the phase behind the hit point grows with n k₀, as along the traced path
        n, opl, δ = 1.5, 0.1, 0.3λ
        p0, s = BMO.Point3(0.0, 0.1, 0.0), BMO.Point3(sin(0.5), cos(0.5), 0.0)
        expected = cis(2π / λ * (opl + n * δ))
        hit = BMO.RayHit{Float64}(p0, s, opl, λ, n, 1.0)
        @test BMO.electric_field(hit, p0 + δ * s) ≈ expected
        hit = BMO.PolarizedRayHit{Float64}(p0, s, opl, λ, n, 1.0, BMO.Point3{ComplexF64}(0, 0, 1))
        @test BMO.electric_field(hit, p0 + δ * s) ≈ [0, 0, expected]
    end

    @testset "Crossing beams add their powers" begin
        # beams at ±2 mrad give fringes of 0.25 mm, whose cross term integrates to zero
        pd = tilted(deg2rad(45))
        α = 2e-3
        system = System([pd])
        for s in (-1, 1)
            dir = [s * sin(α), cos(α), 0]
            solve_system!(system, gb(-100mm / cos(α) * dir + [0, 100mm, 0], dir))
        end
        @test optical_power(pd, RegularGrid(num_points = 800); progress = false) ≈ 2P0 rtol = 1e-3
    end

    @testset "The grid methods sample the automatic grid" begin
        pd = tilted(deg2rad(30))
        solve_system!(System([pd]), agb())
        xs, zs, E = electric_field(pd, RegularGrid(num_points = 20); progress = false)
        @test electric_field(pd, RectilinearGrid(xs, zs); progress = false)[3] == E
        xs, zs, I = intensity(pd, RegularGrid(num_points = 20); progress = false)
        @test intensity(pd, RectilinearGrid(xs, zs); progress = false)[3] == I
    end
end

@testset "A new hit type without readout methods" begin
    p, n̂ = BMO.Point3(0.0, 0.1, 0.0), BMO.Point3(0.0, 1.0, 0.0)
    @test_throws "electric_field not defined for" electric_field(DummyHit(), p)
    @test_throws "intensity not defined for" intensity([DummyHit()], p, n̂)
    @test_throws ErrorException electric_field(DummyHit(), p)
    @test_throws ErrorException intensity([DummyHit()], p, n̂)
    @test_throws "hit_point not defined for" ray_hits(Detector(10mm), [DummyHit()])
end

@testset "Path lengths of a child beamlet equal those of its chief beam" begin
    λ, w0 = 1e-6, 0.5mm
    for (beam, chief) in (
            (GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0), b -> b.chief),
            (AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0), b -> b.c))
        # lens in front of and behind the splitter
        lens = ThinLens(200mm, 200mm, 5mm, 1.5)
        translate3d!(lens, [0, 30mm, 0])
        bs = ThinBeamsplitter(25mm)
        zrotate3d!(bs, π / 4)
        translate3d!(bs, [0, 80mm, 0])
        lens2 = ThinLens(200mm, 200mm, 5mm, 1.5)
        translate3d!(lens2, [0, 130mm, 0])
        solve_system!(System([lens, bs, lens2]), beam)
        @test length(BMO.children(beam)) == 2
        for b in BMO.children(beam), id in eachindex(BMO.rays(b))
            @test BMO.path_excess_before(b, id) == BMO.path_excess_before(chief(b), id)
            @test BMO.length_before(b, id) == BMO.length_before(chief(b), id)
        end
    end
end

end # TESTSET

end # MODULE