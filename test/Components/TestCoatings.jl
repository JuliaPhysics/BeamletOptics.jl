using Test
using BeamletOptics
using LinearAlgebra
const BMO = BeamletOptics

@testset "Coatings and Coated Components" begin
    # Test Uncoated coefficients
    rs_unc, rp_unc, ts_unc, tp_unc = BMO.fresnel_coefficients(
        Uncoated(), 0.0, 1000e-9, 1.0, 1.5)
    @test rs_unc ≈ -0.2
    @test rp_unc ≈ -0.2
    @test ts_unc ≈ 0.8
    @test tp_unc ≈ 0.8

    # Test get_jones_matrix for Uncoated
    J_refl = get_jones_matrix(Uncoated(), 0.0, 1000e-9, 1.0, 1.5, true)
    J_trans = get_jones_matrix(Uncoated(), 0.0, 1000e-9, 1.0, 1.5, false)
    @test J_refl[1, 1] ≈ 0.2
    @test J_refl[2, 2] ≈ -0.2
    @test J_trans[1, 1] ≈ 0.8
    @test J_trans[2, 2] ≈ 0.8

    # Test quarter-wave AR coating at normal incidence
    # n1 = 1.0, n2 = 1.5. Optimal nc = sqrt(1.5).
    # d = lambda / (4 * nc)
    λ = 1000e-9
    nc = sqrt(1.5)
    d = λ / (4 * nc)
    coat_thin = ThinFilmCoating(nc, d)

    rs, rp, ts, tp = BMO.fresnel_coefficients(coat_thin, 0.0, λ, 1.0, 1.5)
    @test abs(rs) < 1e-12
    @test abs(rp) < 1e-12
    @test abs(ts) ≈ sqrt(1.0 / 1.5)  # Power transmission is 100%, amplitude ratio is sqrt(n1/n2)
    @test abs(tp) ≈ sqrt(1.0 / 1.5)

    # Test Lens Interaction with Coating (Bidirectional & Custom Filters)
    n_lens = λ -> 1.5
    base_lens = SphericalLens(100e-3, -100e-3, 5e-3, 10e-3, n_lens)
    coated_lens = base_lens |> with_coatings(front = coat_thin)

    system = System([coated_lens])

    # Trace left-to-right (enters coated front, exits uncoated back)
    ray_lr = PolarizedRay([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0], λ, [1.0, 0.0, 0.0])
    beam_lr = Beam(ray_lr)
    solve_system!(system, beam_lr; retrace = false)

    # Ray 1: before lens, Ray 2: inside (AR coated), Ray 3: outside (uncoated exit)
    ray_inside_lr = rays(beam_lr)[2]
    @test ray_inside_lr.n ≈ 1.5
    @test abs(BMO.polarization(ray_inside_lr)[1]) ≈ sqrt(1.0 / 1.5)

    ray_exit_lr = rays(beam_lr)[3]
    @test ray_exit_lr.n ≈ 1.0
    # Inside polarization is sqrt(1/1.5). Uncoated exit: multiplied by 2 * 1.5 / (1.5 + 1.0) = 1.2.
    @test abs(BMO.polarization(ray_exit_lr)[1]) ≈ sqrt(1.0 / 1.5) * 1.2

    # Trace right-to-left (enters uncoated back, exits coated front)
    ray_rl = PolarizedRay([0.0, 10e-3, 0.0], [0.0, -1.0, 0.0], λ, [1.0, 0.0, 0.0])
    beam_rl = Beam(ray_rl)
    solve_system!(system, beam_rl; retrace = false)

    # Ray 1: before lens, Ray 2: inside (uncoated entry), Ray 3: outside (AR coated exit)
    ray_inside_rl = rays(beam_rl)[2]
    @test ray_inside_rl.n ≈ 1.5
    # Uncoated entry: multiplied by 2 * 1.0 / (1.0 + 1.5) = 0.8
    @test abs(BMO.polarization(ray_inside_rl)[1]) ≈ 0.8

    ray_exit_rl = rays(beam_rl)[3]
    @test ray_exit_rl.n ≈ 1.0
    # AR coated exit: transmittance is 1.0, amplitude ratio factor is sqrt(1.5 / 1.0).
    @test abs(BMO.polarization(ray_exit_rl)[1]) ≈ 0.8 * sqrt(1.5)

    # Test custom normal filter construction
    coating_custom = Coating(
        BMO.shape(base_lens), coat_thin, normal_filter = [0.0, -1.0, 0.0])
    @test coating_custom.normal_filter ≈ [0.0, -1.0, 0.0]

    # Test Mirror Interaction with Coating
    # Coat a flat mirror with a SimpleHRCoating (reflectance = 99%)
    mirror_base = SquarePlanoMirror2D(10e-3)
    coated_mirror = mirror_base |> with_coatings(front = SimpleHRCoating(0.99))

    system_mirror = System([coated_mirror])
    ray_mirror = PolarizedRay([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0], λ, [1.0, 0.0, 0.0])
    beam_mirror = Beam(ray_mirror)
    solve_system!(system_mirror, beam_mirror; retrace = false)

    ray_reflected = rays(beam_mirror)[2]
    @test abs(BMO.polarization(ray_reflected)[1]) ≈ sqrt(0.99)

    # Test Waveplate Coating
    # Construct a JonesCoating representing a quarter-wave plate (retardation pi/2 along fast axis)
    # Global Jones Basis for XZBasis: Ex is fast axis (retardation 0), Ez is slow axis (retardation pi/2)
    wp_jones = BMO.XZBasis(1.0, 0.0, 0.0, exp(im * π / 2))
    coat_wp = JonesCoating(wp_jones, BMO.XZBasis(0, 0, 0, 0))

    # Plano-planar lens of index 1.0 (air) coated on its front surface
    n_air = λ -> 1.0
    base_planar = SphericalLens(Inf, Inf, 1e-3, 10e-3, n_air)
    coated_planar = base_planar |> with_coatings(front = coat_wp)

    system_wp = System([coated_planar])
    # Trace a linearly polarized ray with equal components in Ex and Ez (45 degrees)
    ray_wp = PolarizedRay(
        [0.0, -5e-3, 0.0], [0.0, 1.0, 0.0], λ, [1 / sqrt(2), 0.0, 1 / sqrt(2)])
    beam_wp = Beam(ray_wp)
    solve_system!(system_wp, beam_wp; retrace = false)

    ray_out = rays(beam_wp)[2] # after first surface
    pol_out = BMO.polarization(ray_out)
    @test pol_out[1] ≈ 1 / sqrt(2)
    @test pol_out[3] ≈ im / sqrt(2)

    # Test Coated Cube Beamsplitter
    # Custom coating on the splitter interface (transmittance 80%, reflectance 20%)
    coat_bs = SimpleBeamsplitterCoating(sqrt(0.2), sqrt(0.2), sqrt(0.8), sqrt(0.8))

    flat_mesh = BMO.QuadraticFlatMesh(10e-3)
    coating_obj = Coating(flat_mesh, coat_bs)

    system_bs = System([coating_obj])
    ray_bs = PolarizedRay([0.0, -5e-3, 0.0], [0.0, 1.0, 0.0], λ, [1.0, 0.0, 0.0])
    beam_bs = Beam(ray_bs)
    solve_system!(system_bs, beam_bs; retrace = false)

    # We should have child 1 (transmitted) and child 2 (reflected) from the splitter
    @test length(beam_bs.children) == 2
    ray_t = first(rays(beam_bs.children[1]))
    @test abs(BMO.polarization(ray_t)[1]) ≈ sqrt(0.8)

    ray_r = first(rays(beam_bs.children[2]))
    @test abs(BMO.polarization(ray_r)[1]) ≈ sqrt(0.2)

    # Test Coated Doublet Lens (Cemented Coating)
    n1_doublet = λ -> 1.5
    n2_doublet = λ -> 1.6

    dl_front = SphericalLens(Inf, -100e-3, 5e-3, 10e-3, n1_doublet)
    dl_back = SphericalLens(-100e-3, Inf, 5e-3, 10e-3, n2_doublet)

    # Coat the cemented boundary between dl_front and dl_back
    coat_cement = SimpleARCoating(0.01) # 1% reflectance

    dl_front_coated = with_coatings(dl_front, back = coat_cement)

    # Move back lens flush with front lens
    translate3d!(dl_back, [0, thickness(dl_front), 0])

    system_dl = System([dl_front_coated, dl_back])
    ray_dl = PolarizedRay([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0], λ, [1.0, 0.0, 0.0])
    beam_dl = Beam(ray_dl)
    solve_system!(system_dl, beam_dl; retrace = false)

    # Verify index and polarization amplitude scaling at cemented interface
    @test length(rays(beam_dl)) >= 3
    ray_inside_back = rays(beam_dl)[3]
    @test ray_inside_back.n ≈ 1.6

    pol_inside_back = BMO.polarization(ray_inside_back)
    @test abs(pol_inside_back[1]) ≈ 0.8 * sqrt(0.99)

    # Test Segmented Coating (Spatial Predicates)
    # Lens with left-half AR coated, right-half HR coated
    left_half = (p, n) -> p[1] < 0.0
    right_half = (p, n) -> p[1] >= 0.0

    ar_model = SimpleARCoating(0.0)
    hr_model = SimpleHRCoating(1.0)

    base_lens_flat = SphericalLens(Inf, Inf, 1e-3, 10e-3, n_air)
    coated_segmented = with_coatings(
        base_lens_flat,
        left_half => ar_model,
        right_half => hr_model
    )

    system_seg = System([coated_segmented])

    # Ray 1: hits left half (transmissive segment)
    ray_seg_l = Ray([-2e-3, -5e-3, 0.0], [0.0, 1.0, 0.0], λ)
    beam_seg_l = Beam(ray_seg_l)
    solve_system!(system_seg, beam_seg_l; retrace = false)
    # Should transmit through: pos after lens is around y = 0.005
    @test length(rays(beam_seg_l)) == 3
    @test rays(beam_seg_l)[3].pos[2] > 0.0

    # Ray 2: hits right half (reflective segment)
    ray_seg_r = Ray([2e-3, -5e-3, 0.0], [0.0, 1.0, 0.0], λ)
    beam_seg_r = Beam(ray_seg_r)
    solve_system!(system_seg, beam_seg_r; retrace = false)
    # Should reflect back: direction is [0.0, -1.0, 0.0]
    @test length(rays(beam_seg_r)) == 2
    @test rays(beam_seg_r)[2].dir[2] ≈ -1.0

    # Test Face-Selective Coating on a Right Angle Prism
    n_prism = λ -> 1.5
    prism_base = RightAnglePrism(10e-3, 10e-3, n_prism)
    # Coat the hypotenuse with HR, keep other legs uncoated
    coated_prism = with_coatings(prism_base, :hypotenuse => hr_model)

    system_prism = System([coated_prism])
    # Send a ray entering from Leg 2 (y=-5e-3), hitting hypotenuse (x+y=0), reflecting out of Leg 1 (x=-5e-3)
    ray_prism = Ray([-2e-3, -10e-3, 0.0], [0.0, 1.0, 0.0], λ)
    beam_prism = Beam(ray_prism)
    solve_system!(system_prism, beam_prism; retrace = false)

    # Ray 1: before prism (starts at y=-10e-3)
    # Ray 2: inside prism (enters y=-5e-3, goes to y=2e-3)
    # Ray 3: inside prism (reflects to [-1, 0, 0], goes to x=-5e-3)
    # Ray 4: after prism (exits Leg 1, goes to x < -5e-3)
    @test length(rays(beam_prism)) == 4
    @test rays(beam_prism)[3].dir[1] ≈ -1.0
    @test rays(beam_prism)[4].dir[1] ≈ -1.0

    @testset "Phase Convention Validation" begin
        # Let's define a plano-planar lens of glass (n = 1.5)
        n_glass = λ -> 1.5
        planar_lens = SphericalLens(Inf, Inf, 1e-3, 10e-3, n_glass)

        # Splitting coating (bare index jump splitting)
        splitting_coating = ThinFilmCoating(
            Float64[], Float64[]; behavior = BMO.Splitting())

        # External Normal Reflection (Air to Glass)
        coated_lens_ext = with_coatings(planar_lens, front = splitting_coating)
        system_ext = System([coated_lens_ext])

        gauss_ext = GaussianBeamlet([0.0, -5e-3, 0.0], [0.0, 1.0, 0.0], 1000e-9, 1e-3)
        solve_system!(system_ext, gauss_ext; retrace = false)

        @test length(gauss_ext.children) == 2
        t_ext = gauss_ext.children[1]
        r_ext = gauss_ext.children[2]
        @test electric_field(r_ext) / electric_field(gauss_ext)≈-0.2 atol=1e-8

        # Internal Normal Reflection (Glass to Air)
        coated_lens_int = with_coatings(planar_lens, back = splitting_coating)
        system_int = System([coated_lens_int])

        gauss_int = GaussianBeamlet([0.0, -5e-3, 0.0], [0.0, 1.0, 0.0], 1000e-9, 1e-3)
        solve_system!(system_int, gauss_int; retrace = false)

        @test length(gauss_int.children) == 2
        t_int = gauss_int.children[1]
        r_int = gauss_int.children[2]
        @test electric_field(r_int) / electric_field(gauss_int)≈0.2 atol=1e-8

        # Consistency of Polarization and Amplitude for AstigmaticGaussianBeamlet
        p_ray_ext = PolarizedRay(
            [0.0, -5e-3, 0.0], [0.0, 1.0, 0.0], 1000e-9, [1.0, 0.0, 0.0])
        p_beam = Beam(p_ray_ext)
        solve_system!(system_ext, p_beam; retrace = false)

        @test length(p_beam.children) == 2
        p_r = first(rays(p_beam.children[2]))

        agb_ext = AstigmaticGaussianBeamlet(
            [0.0, -5e-3, 0.0], [0.0, 1.0, 0.0], 1000e-9, 1e-3;
            E0 = [1.0, 0.0, 0.0], support = [0.0, 0.0, 1.0])
        solve_system!(system_ext, agb_ext; retrace = false)

        @test length(agb_ext.children) == 2
        agb_r = agb_ext.children[2]
        agb_r_chief = first(rays(agb_r.c))
        @test agb_r_chief.E0≈p_r.E0 atol=1e-8
    end

    @testset "Complex Refractive Index ThinFilmCoating" begin
        # Test a thin metallic coating (e.g. chromium: n = 3.0 + 4.0im)
        λ = 1000e-9
        n_metal = 3.0 + 4.0im
        d_metal = 10e-9 # 10 nm physical thickness
        coat_metal = ThinFilmCoating(n_metal, d_metal; behavior = BMO.Reflective())

        rs, rp, ts, tp = BMO.fresnel_coefficients(coat_metal, 0.0, λ, 1.0, 1.5)
        # Verify coefficients are complex numbers
        @test rs isa Complex
        @test rp isa Complex
        @test ts isa Complex
        @test tp isa Complex
        # Non-zero reflection/transmission due to finite thickness
        @test abs(rs) > 0.0
        @test abs(ts) > 0.0
    end

    @testset "TIR on Splitting Boundaries" begin
        # Plano-planar lens (n = 1.5) with back surface coated with splitting coating
        n_glass = λ -> 1.5
        planar_lens = SphericalLens(Inf, Inf, 1e-3, 10e-3, n_glass)
        splitting_coating = ThinFilmCoating(
            Float64[], Float64[]; behavior = BMO.Splitting())
        coated_lens = with_coatings(planar_lens, back = splitting_coating)
        system = System([coated_lens])

        # Single Ray starting inside the glass lens heading towards the back surface at 45 degrees
        # Normal to back surface is [0.0, 1.0, 0.0]
        # Ray direction: normalize([1.0, 1.0, 0.0]), which has 45 deg angle of incidence
        dir_45 = normalize([1.0, 1.0, 0.0])
        ray = Ray{Float64}([0.0, -0.4e-3, 0.0], dir_45, nothing, 1000e-9, 1.5)
        beam = Beam(ray)
        solve_system!(system, beam; retrace = false)

        # Under TIR, the ray should reflect (direction [1/√2, -1/√2, 0.0]) and stay in glass (n=1.5)
        # Then it hits the front surface (normal [0.0, -1.0, 0.0]) and exits into air.
        @test length(rays(beam)) >= 3
        refl_ray = rays(beam)[3]
        @test refl_ray.dir[1] ≈ 1 / sqrt(2)
        @test refl_ray.dir[2] ≈ -1 / sqrt(2)
        @test refl_ray.n ≈ 1.5
        @test isempty(beam.children) # No children created!

        # Polarized Ray starting inside the glass lens
        pray = PolarizedRay{Float64}(
            [0.0, -0.4e-3, 0.0], dir_45, nothing, 1000e-9, 1.5, [0.0, 0.0, 1.0])
        pbeam = Beam(pray)
        solve_system!(system, pbeam; retrace = false)
        @test length(rays(pbeam)) >= 3
        refl_pray = rays(pbeam)[3]
        @test refl_pray.dir[1] ≈ 1 / sqrt(2)
        @test refl_pray.dir[2] ≈ -1 / sqrt(2)
        @test refl_pray.n ≈ 1.5
        @test isempty(pbeam.children)

        # Gaussian Beamlet starting inside the glass lens
        chief_ray = Ray{Float64}([0.0, -0.4e-3, 0.0], dir_45, nothing, 1000e-9, 1.5)
        waist_ray = Ray{Float64}([0.0, -0.4e-3, 0.1e-3], dir_45, nothing, 1000e-9, 1.5)
        div_ray = Ray{Float64}(
            [0.0, -0.4e-3, 0.0], normalize([1.0, 1.0, 0.01]), nothing, 1000e-9, 1.5)
        chief_beam = Beam(chief_ray)
        waist_beam = Beam(waist_ray)
        div_beam = Beam(div_ray)
        gauss = GaussianBeamlet(
            chief_beam, waist_beam, div_beam, 1000e-9, 1e-3, 1.0 + 0.0im)
        solve_system!(system, gauss; retrace = false)
        @test isempty(gauss.children) # Should not split under TIR!
    end

    @testset "Frustrated Total Internal Reflection (FTIR) on Splitting Boundaries" begin
        # Glass index n=1.5
        n_glass = λ -> 1.5

        # 100 nm air gap splitting coating (index 1.0, thickness 100nm)
        ftir_coating = ThinFilmCoating(1.0, 100e-9; behavior = BMO.Splitting())

        # Two plano-planar lenses made of glass
        lens_front = SphericalLens(Inf, Inf, 1e-3, 10e-3, n_glass)
        lens_back = SphericalLens(Inf, Inf, 1e-3, 10e-3, n_glass)

        # Coat the back of the front lens with our air gap coating
        dl_front_coated = with_coatings(lens_front; back = ftir_coating)

        # Place the back lens right after the front lens (thickness is 1mm)
        translate3d!(lens_back, [0.0, 1e-3, 0.0])

        system = System([dl_front_coated, lens_back])

        # Launch a ray inside the front lens hitting the coated boundary at 45 degrees
        # (Incident angle 45 deg > critical angle 41.8 deg, but because the transmitted
        # medium is also glass (n=1.5), it does not trigger the bulk TIR check).
        dir_45 = normalize([1.0, 1.0, 0.0])
        ray = Ray{Float64}([0.0, -0.4e-3, 0.0], dir_45, nothing, 1000e-9, 1.5)
        beam = Beam(ray)

        solve_system!(system, beam; retrace = false, depth_max = 1)

        # Verify that FTIR allows both reflection and transmission (splitting happens)
        @test length(beam.children) == 2

        t_child = beam.children[1]
        r_child = beam.children[2]

        @test length(BMO.rays(t_child)) >= 1
        @test length(BMO.rays(r_child)) >= 1

        t_ray = first(BMO.rays(t_child))
        r_ray = first(BMO.rays(r_child))

        # Verify weight/power distribution. T and R should be non-zero and sum to ~1.
        # Since T is non-zero, evanescent coupling/frustration is working.
        T_weight = BMO.weight(t_ray)
        R_weight = BMO.weight(r_ray)
        @test T_weight > 0.0
        @test R_weight > 0.0
        @test T_weight + R_weight≈1.0 atol=1e-5
    end

    @testset "Backwards Ray Coincident Boundary Bug" begin
        # Issue: Rays tracing backwards through a coincident coated boundary
        # had flipped refractive indices because of `isentering` relying on intersection primary object.
        n1_doublet = λ -> 1.5
        n2_doublet = λ -> 1.6

        dl_front = SphericalLens(Inf, -100e-3, 5e-3, 10e-3, n1_doublet)
        dl_back = SphericalLens(-100e-3, Inf, 5e-3, 10e-3, n2_doublet)

        coat_cement = SimpleARCoating(0.01)
        dl_front_coated = with_coatings(dl_front; back = coat_cement)

        translate3d!(dl_back, [0, thickness(dl_front), 0])

        # Reverse order to ensure dl_back is the primary coincident object
        system_rev = System([dl_back, dl_front_coated])

        # Ray starts behind the back lens and travels backwards
        ray_bwd = PolarizedRay(
            [0.0, 12e-3, 0.0], [0.0, -1.0, 0.0], 1000e-9, [1.0, 0.0, 0.0])
        beam_bwd = Beam(ray_bwd)
        solve_system!(system_rev, beam_bwd; retrace = false)

        @test length(rays(beam_bwd)) >= 3
        # Ray 1: Air to dl_back (y=12e-3 to 10e-3)
        # Ray 2: Inside dl_back (y=10e-3 to 5e-3)
        # Ray 3: Inside dl_front (y=5e-3 to 0.0)

        r2 = rays(beam_bwd)[2]
        @test r2.n ≈ 1.6

        r3 = rays(beam_bwd)[3]
        @test r3.n ≈ 1.5

        # Verify amplitude scaling: transmission from 1.0 to 1.6 (Uncoated), then 1.6 to 1.5 (SimpleARCoating)
        expected_amp = (2 * 1.0 / (1.0 + 1.6)) * sqrt(0.99)
        pol3 = BMO.polarization(r3)
        @test abs(pol3[1]) ≈ expected_amp
    end

    @testset "Flat shape coatings checks" begin
        # Test degenerate collinear vertex check in is_flat_shape
        mesh_3d = BMO.CuboidMesh(1.0, 1.0, 1.0)
        # Force the first three vertices to be collinear
        mesh_3d.vertices[1, :] = [0.0, 0.0, 0.0]
        mesh_3d.vertices[2, :] = [1.0, 0.0, 0.0]
        mesh_3d.vertices[3, :] = [2.0, 0.0, 0.0]
        @test !BMO.is_flat_shape(mesh_3d)

        # Test ArgumentError on flat shape with both front and back coatings
        flat_mesh = BMO.QuadraticFlatMesh(10e-3)
        coat1 = SimpleARCoating(0.0)
        coat2 = SimpleHRCoating(1.0)

        mirror = Mirror(flat_mesh)
        @test_throws ArgumentError with_coatings(mirror, front = coat1, back = coat2)
    end

    @testset "Show methods and Display" begin
        repr_plain(x) = repr(MIME"text/plain"(), x)
        # Test models show
        @test repr_plain(Uncoated()) == "Uncoated"
        @test repr_plain(SimpleARCoating(0.05)) == "SimpleARCoating(R = 0.05)"
        @test repr_plain(SimpleHRCoating(0.98)) == "SimpleHRCoating(R = 0.98)"
        @test repr_plain(SimpleBeamsplitterCoating(0.5, 0.5, 0.5, 0.5)) ==
              "SimpleBeamsplitterCoating(rs = 0.5 + 0.0im, rp = 0.5 + 0.0im, ts = 0.5 + 0.0im, tp = 0.5 + 0.0im)"
        @test repr_plain(JonesCoating(BMO.SPBasis(1, 0, 0, 1))) ==
              "JonesCoating(behavior = Transmissive())"
        @test repr_plain(ThinFilmCoating(1.38, 200e-9)) ==
              "ThinFilmCoating(1 layers, behavior = Transmissive())"

        # Test coated components show
        lens = SphericalLens(100e-3, -100e-3, 5e-3, 10e-3, λ -> 1.5)
        cl = with_coatings(lens, front = SimpleARCoating(0.02))
        @test cl isa Lens
        @test BMO.coatings(cl) != ()

        mirror = SquarePlanoMirror2D(10e-3)
        cm = with_coatings(mirror, front = SimpleHRCoating(0.99))
        @test cm isa Mirror
        @test BMO.coatings(cm) != ()
    end

    @testset "Absorptive behavior" begin
        # Standalone absorptive coating
        struct DummyAbsorptiveCoating
            behavior::Absorptive
        end
        BMO.coating_behavior(::DummyAbsorptiveCoating) = Absorptive()

        flat_mesh = BMO.QuadraticFlatMesh(10e-3)
        abs_coat = Coating(flat_mesh, DummyAbsorptiveCoating(Absorptive()))

        system = System([abs_coat])
        ray = Ray([0.0, -5e-3, 0.0], [0.0, 1.0, 0.0], 1000e-9)
        beam = Beam(ray)
        solve_system!(system, beam; retrace = false)
        # Ray should be terminated (absorbed) at the interface, so only 1 ray exists in the beam
        @test length(rays(beam)) == 1

        # CoatedRefractive with absorptive coating
        cl_abs = SphericalLens(100e-3, -100e-3, 5e-3, 10e-3, λ -> 1.5) |>
                 with_coatings(front = DummyAbsorptiveCoating(Absorptive()))
        system_cl = System([cl_abs])
        beam_cl = Beam(Ray([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0], 1000e-9))
        solve_system!(system_cl, beam_cl; retrace = false)
        @test length(rays(beam_cl)) == 1

        # CoatedMirror with absorptive coating
        cm_abs = SquarePlanoMirror2D(10e-3) |>
                 with_coatings(front = DummyAbsorptiveCoating(Absorptive()))
        system_cm = System([cm_abs])
        beam_cm = Beam(Ray([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0], 1000e-9))
        solve_system!(system_cm, beam_cm; retrace = false)
        @test length(rays(beam_cm)) == 1
    end

    @testset "Wavelength-dependent coatings and dynamic behavior" begin
        # Dynamic behavior based on ray properties
        struct DynamicARCoating end
        # behaves as Transmissive for λ < 800nm, Reflective for λ >= 800nm
        BMO.coating_behavior(::DynamicARCoating, ray) = BMO.wavelength(ray) < 800e-9 ?
                                                        Transmissive() : Reflective()
        BMO.get_jones_matrix(::DynamicARCoating, θi, λ, n1, n2, is_reflected; from_front = true, kwargs...) = BMO.SPBasis(
            1.0, 0, 0, 1.0)

        flat_mesh = BMO.QuadraticFlatMesh(10e-3)
        dyn_coat = Coating(flat_mesh, DynamicARCoating())
        system = System([dyn_coat])

        # λ = 600nm -> should transmit
        beam_t = Beam(Ray([0.0, -5e-3, 0.0], [0.0, 1.0, 0.0], 600e-9))
        solve_system!(system, beam_t; retrace = false)
        @test length(rays(beam_t)) == 2
        @test rays(beam_t)[2].dir[2] ≈ 1.0

        # λ = 900nm -> should reflect
        beam_r = Beam(Ray([0.0, -5e-3, 0.0], [0.0, 1.0, 0.0], 900e-9))
        solve_system!(system, beam_r; retrace = false)
        @test length(rays(beam_r)) == 2
        @test rays(beam_r)[2].dir[2] ≈ -1.0

        # ThinFilmCoating with dispersion function n(λ)
        n_disp = λ -> λ < 800e-9 ? 1.38 : 1.6
        coat_disp = ThinFilmCoating(n_disp, 200e-9)
        rs_short, _, _, _ = fresnel_coefficients(coat_disp, 0.0, 600e-9, 1.0, 1.5)
        rs_long, _, _, _ = fresnel_coefficients(coat_disp, 0.0, 900e-9, 1.0, 1.5)
        @test rs_short != rs_long
    end

    @testset "Coated Doublet Lens" begin
        n1 = 1.5
        n2 = 1.6
        ar = SimpleARCoating(0.0)

        doublet = SphericalDoubletLens(50e-3, -50e-3, -100e-3, 5e-3, 3e-3, 25e-3, n1,
            n2; front_coating = ar, back_coating = ar)
        @test doublet.front isa Lens
        @test doublet.back isa Lens
        @test BMO.coatings(doublet.front) != ()
        @test BMO.coatings(doublet.back) != ()

        sys = System([doublet])
        beam = Beam(Ray([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0], 589e-9))
        solve_system!(sys, beam; retrace = false)
        @test length(rays(beam)) == 4 # Initial ray + 3 surface transitions
    end

    @testset "CSG UnionSDF Coating Resolution" begin
        lens1 = SphericalLens(50e-3, -50e-3, 5e-3, 10e-3, 1.5)
        lens2 = SphericalLens(100e-3, -100e-3, 5e-3, 10e-3, 1.5)
        translate3d!(lens2, [0.0, 20e-3, 0.0])

        s1 = BMO.shape(lens1)
        s2 = BMO.shape(lens2)
        u_shape = s1 + s2

        p1 = BMO.position(s1) + [0.0, -2.5e-3, 0.0]
        p2 = BMO.position(s2) + [0.0, -2.5e-3, 0.0]

        @test BMO.active_constituent_sdf(u_shape, p1) isa BMO.AbstractSphericalSurfaceSDF
        @test BMO.active_constituent_sdf(u_shape, p2) isa BMO.AbstractSphericalSurfaceSDF

        n_glass = λ -> 1.5
        struct CompoundLens{T, C <: Tuple} <:
               BMO.AbstractRefractiveOptic{T, BMO.RefractiveIndex}
            shape::BMO.UnionSDF{T}
            n::Function
            coatings::C
            function CompoundLens(shape::BMO.UnionSDF{T}, n::Function,
                    coatings::C = ()) where {T, C <: Tuple}
                return new{T, C}(shape, n, coatings)
            end
        end
        BMO.shape_trait_of(::CompoundLens) = BMO.SingleShape()
        BMO.shape(c::CompoundLens) = c.shape
        BMO.refractive_index(c::CompoundLens, λ::Real) = c.n(λ)
        BMO._attach_coatings(c::CompoundLens, c_tuple; deepcopy_shape::Bool = false) = CompoundLens(
            deepcopy_shape ? deepcopy(c.shape) : c.shape, c.n, c_tuple)

        compound = CompoundLens(u_shape, n_glass)
        coated_compound = compound |> with_coatings(front = SimpleARCoating(0.02))

        sys = System([coated_compound])
        beam = Beam(Ray([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0], 1000e-9))
        solve_system!(sys, beam; retrace = false)

        model1 = BMO.get_coating_model_at_hit(coated_compound, rays(beam)[1])
        @test model1 isa SimpleARCoating
        @test model1.R ≈ 0.02
    end

    @testset "Unified Coating Transmittance and Reflectance API" begin
        ar = SimpleARCoating(0.01) # R=1%
        R = coating_reflectance(ar, 0.0, 1064e-9, 1.0, 1.5)
        T = coating_transmittance(ar, 0.0, 1064e-9, 1.0, 1.5)
        @test R≈0.01 atol=1e-6
        @test T≈0.99 atol=1e-6
        @test R + T≈1.0 atol=1e-6
    end

    @testset "Analytic Normals and Surface Tags" begin
        box = BMO.BoxSDF(10e-3, 20e-3, 30e-3)
        p_front = BMO.Point3(0.0, -10e-3, 0.0)
        n_analytic_front = BMO.normal3d(box, p_front)
        n_fd_front = BMO.normal_fd(box, p_front)
        @test n_analytic_front≈n_fd_front atol=1e-4
        @test BMO.surface_tag(box, p_front, n_analytic_front) === :front

        p_back = BMO.Point3(0.0, 10e-3, 0.0)
        n_analytic_back = BMO.normal3d(box, p_back)
        @test BMO.surface_tag(box, p_back, n_analytic_back) === :back

        cyl = BMO.CylinderSDF(5e-3, 10e-3)
        p_top = BMO.Point3(0.0, 5e-3, 0.0)
        n_top = BMO.normal3d(cyl, p_top)
        @test n_top≈[0.0, 1.0, 0.0] atol=1e-4
        @test BMO.surface_tag(cyl, p_top, n_top) === :top

        p_side = BMO.Point3(5e-3, 0.0, 0.0)
        n_side = BMO.normal3d(cyl, p_side)
        @test n_side≈[1.0, 0.0, 0.0] atol=1e-4
        @test BMO.surface_tag(cyl, p_side, n_side) === :side

        sph = BMO.SphereSDF(10e-3)
        p_sph_front = BMO.Point3(0.0, -10e-3, 0.0)
        n_sph_front = BMO.normal3d(sph, p_sph_front)
        @test n_sph_front≈[0.0, -1.0, 0.0] atol=1e-4
        @test BMO.surface_tag(sph, p_sph_front, n_sph_front) === :front

        p_sph_back = BMO.Point3(0.0, 10e-3, 0.0)
        n_sph_back = BMO.normal3d(sph, p_sph_back)
        @test BMO.surface_tag(sph, p_sph_back, n_sph_back) === :back
    end

    @testset "AbstractSurfaceModel Hierarchy" begin
        @test SimpleARCoating() isa AbstractCoatingModel
        @test SimpleARCoating() isa AbstractSurfaceModel
        @test Uncoated() isa AbstractSurfaceModel
        @test ThinFilmCoating([], []) isa AbstractSurfaceModel
    end

    @testset "Multi-Face Selector Syntax" begin
        lens = BMO.SphericalLens(25e-3, 50e-3, 50e-3)
        coated = BMO.with_coatings(lens, (:front, :back) => BMO.SimpleARCoating(0.01))

        p_front = BMO.Point3(0.0, -2e-3, 0.0)
        p_back = BMO.Point3(0.0, 2e-3, 0.0)
        n_front = BMO.Point3(0.0, -1.0, 0.0)
        n_back = BMO.Point3(0.0, 1.0, 0.0)

        @test BMO.get_matching_coating(
            BMO.coatings(coated), BMO.shape(coated), p_front, n_front) isa
              BMO.SimpleARCoating
        @test BMO.get_matching_coating(
            BMO.coatings(coated), BMO.shape(coated), p_back, n_back) isa BMO.SimpleARCoating
    end

    @testset "CompositeSurfaceModel" begin
        ar1 = BMO.SimpleARCoating(0.05) # R=0.05, T=0.95
        ar2 = BMO.SimpleARCoating(0.10) # R=0.10, T=0.90
        comp = BMO.CompositeSurfaceModel(ar1, ar2)
        @test BMO.coating_behavior(comp) isa BMO.Transmissive

        J = BMO.get_jones_matrix(comp, 0.0, 632.8e-9, 1.0, 1.5, false)
        T_val = BMO.unpolarized_transmittance(J)
        @test T_val≈(0.95 * 0.90) atol=1e-4
    end

    @testset "deepcopy_shape option in with_coatings" begin
        lens = BMO.SphericalLens(25e-3, 50e-3, 50e-3)
        c1 = BMO.with_coatings(
            lens, :front => BMO.SimpleARCoating(0.01); deepcopy_shape = false)
        c2 = BMO.with_coatings(
            lens, :front => BMO.SimpleARCoating(0.01); deepcopy_shape = true)

        @test c1.shape === lens.shape
        @test c2.shape !== lens.shape
    end

    @testset "GradedThinFilmCoating" begin
        tf = BMO.ThinFilmCoating([1.38], [150e-9])
        graded = BMO.GradedThinFilmCoating(p -> 1.0 + 0.1 * p[1], tf)
        @test BMO.coating_behavior(graded) isa BMO.Transmissive

        p1 = BMO.Point3(0.0, 0.0, 0.0)
        p2 = BMO.Point3(1.0, 0.0, 0.0)

        J1 = BMO.get_jones_matrix(graded, 0.0, 632.8e-9, 1.0, 1.5, false; local_p = p1)
        J2 = BMO.get_jones_matrix(graded, 0.0, 632.8e-9, 1.0, 1.5, false; local_p = p2)

        @test BMO.unpolarized_transmittance(J1) != BMO.unpolarized_transmittance(J2)
    end

    @testset "Coating Engine Bug Prevention & Surface Tag Contracts" begin
        # Test surface_tag 1-arg vs 3-arg contract on translated/rotated Lens SDFs
        lens = SphericalLens(25.0e-3, 50.0e-3, 5.0e-3, 10.0e-3) # thickness = 5mm, y ranges from 0 to 5mm
        translate3d!(lens, [1.0e-3, 2.0e-3, 3.0e-3])
        zrotate3d!(lens, π / 4)

        shape_lens = BMO.shape(lens)
        world_front_pt = position(lens) + orientation(lens) * BMO.Point3(0.0, -1.0e-3, 0.0) # y < 0
        world_back_pt = position(lens) + orientation(lens) * BMO.Point3(0.0, 6.0e-3, 0.0) # y > 5mm

        # 1-arg surface_tag takes world coordinates and returns correct tag
        @test BMO.surface_tag(shape_lens, world_front_pt) === :front
        @test BMO.surface_tag(shape_lens, world_back_pt) === :back

        # 3-arg surface_tag takes local coordinates and returns identical tag
        local_front_pt = BMO._world_to_sdf(shape_lens, world_front_pt)
        local_back_pt = BMO._world_to_sdf(shape_lens, world_back_pt)
        @test BMO.surface_tag(shape_lens, local_front_pt, BMO.Point3(0.0, -1.0, 0.0)) ===
              :front
        @test BMO.surface_tag(shape_lens, local_back_pt, BMO.Point3(0.0, 1.0, 0.0)) ===
              :back

        # Test spatial hit position resolution in resolve_coated_boundary
        # Front face has AR coating (R=0), back face has Beamsplitter coating (R=0.8)
        ar_coat = SimpleARCoating(0.0)
        bs_coat = SimpleBeamsplitterCoating(0.8, 0.8, 0.2, 0.2)
        coated_lens = lens |> with_coatings(front = ar_coat, back = bs_coat)

        sys = System([coated_lens])

        # Ray starting outside front face targeting lens (s-polarized along z-axis)
        ray_dir = orientation(lens) * BMO.Point3(0.0, 1.0, 0.0)
        ray_in = PolarizedRay(
            world_front_pt - 10.0e-3 * ray_dir, ray_dir, 632.8e-9, [0.0, 0.0, 1.0])
        beam_in = Beam(ray_in)
        solve_system!(sys, beam_in; retrace = false)
        ray_hit = BMO.rays(beam_in)[1]

        # Test boundary resolution at hit point (front face)
        _, coating_front = BMO.resolve_coated_boundary(sys, coated_lens, ray_hit)
        @test coating_front === ar_coat

        # Substrate object identity matching in interface context
        # Check that shape matching detects substrate correctly (6th return element is is_substrate)
        ctx = BMO._resolve_interface_context(sys, coated_lens, ray_hit)
        @test ctx[6] === true # is_substrate is true for coated_lens

        # Fabry-Pérot cavity tracing onto detector test
        m1 = SphericalLens(Inf, Inf, 0.5e-3, 5.0e-3) |>
             with_coatings(front = SimpleBeamsplitterCoating(0.8, 0.8, 0.2, 0.2),
            back = SimpleARCoating(0.0))
        m2 = SphericalLens(Inf, Inf, 0.5e-3, 5.0e-3) |>
             with_coatings(front = SimpleBeamsplitterCoating(0.8, 0.8, 0.2, 0.2),
            back = SimpleARCoating(0.0))
        translate3d!(m2, [0.0, 1.0e-3, 0.0])

        det = Detector(10.0e-3)
        translate3d!(det, [0.0, 5.0e-3, 0.0])

        fp_sys = System([m1, m2, det])
        fp_beam = Beam(PolarizedRay(
            [0.0, -1.0e-3, 0.0], [0.0, 1.0, 0.0], 632.8e-9, [0.0, 0.0, 1.0]))
        solve_system!(fp_sys, fp_beam; r_max = 50, depth_max = 12, retrace = false)

        @test !isnothing(BMO.hits(det)) && length(BMO.hits(det)) > 0
    end

    @testset "Audit Improvements: Layer TIR Robustness & AD Compatibility" begin
        # Real-index layer TIR in TMM without DomainError
        # Layer index nl = 1.2, incident n1 = 1.5, substrate n2 = 1.5.
        # At θi = 60° (sin 60° ≈ 0.866), n1 * sinθi = 1.299 > 1.2 -> TIR occurs inside layer!
        θ_tir = π / 3
        coat_low_index = ThinFilmCoating(1.2, 500e-9)
        rs, rp, ts, tp = BMO.fresnel_coefficients(coat_low_index, θ_tir, 1000e-9, 1.5, 1.5)
        @test rs isa Complex
        @test rp isa Complex
        @test ts isa Complex
        @test tp isa Complex
        # Energy conservation check
        R_val = coating_reflectance(coat_low_index, θ_tir, 1000e-9, 1.5, 1.5)
        T_val = coating_transmittance(coat_low_index, θ_tir, 1000e-9, 1.5, 1.5)
        @test R_val + T_val≈1.0 atol=1e-5

        # Parametric Coatings
        ar_f32 = SimpleARCoating(Float32(0.01))
        @test ar_f32.R isa Float32
        hr_f32 = SimpleHRCoating(Float32(0.99))
        @test hr_f32.R isa Float32
        bs_gen = SimpleBeamsplitterCoating(0.5, 0.5, 0.5, 0.5)
        @test bs_gen.rs isa ComplexF64

        # ForwardDiff / Automatic Differentiation through ThinFilmCoating
        using ForwardDiff
        # Compute gradient of transmittance with respect to layer thickness d
        λ_test = 632.8e-9
        calc_trans = d -> begin
            c = ThinFilmCoating(1.38, d)
            coating_transmittance(c, 0.0, λ_test, 1.0, 1.5)
        end
        # Quarter wave thickness: d0 = λ / (4 * 1.38)
        d0 = λ_test / (4 * 1.38)
        grad_d = ForwardDiff.derivative(calc_trans, d0)
        # At optimal AR thickness, transmittance is at a local extremum -> derivative is ~0
        @test abs(grad_d) < 1e-4
    end
end
