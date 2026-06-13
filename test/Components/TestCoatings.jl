using Test
using BeamletOptics
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

    dl_front_coated = CoatedLens(dl_front; back = coat_cement)

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
    coated_segmented = CoatedLens(
        base_lens_flat, Pair{Any, Any}[
            left_half => ar_model,
            right_half => hr_model
        ])

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
    coated_prism = CoatedLens(prism_base, Pair{Any, Any}[
        :hypotenuse => hr_model
    ])

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
        splitting_coating = ThinFilmCoating(Float64[], Float64[]; behavior = BMO.Splitting())
        
        # External Normal Reflection (Air to Glass)
        coated_lens_ext = CoatedLens(planar_lens, front = splitting_coating)
        system_ext = System([coated_lens_ext])
        
        gauss_ext = GaussianBeamlet([0.0, -5e-3, 0.0], [0.0, 1.0, 0.0], 1000e-9, 1e-3)
        solve_system!(system_ext, gauss_ext; retrace = false)
        
        @test length(gauss_ext.children) == 2
        t_ext = gauss_ext.children[1]
        r_ext = gauss_ext.children[2]
        @test electric_field(r_ext) / electric_field(gauss_ext) ≈ -0.2 atol=1e-8
        
        # Internal Normal Reflection (Glass to Air)
        coated_lens_int = CoatedLens(planar_lens, back = splitting_coating)
        system_int = System([coated_lens_int])
        
        gauss_int = GaussianBeamlet([0.0, -5e-3, 0.0], [0.0, 1.0, 0.0], 1000e-9, 1e-3)
        solve_system!(system_int, gauss_int; retrace = false)
        
        @test length(gauss_int.children) == 2
        t_int = gauss_int.children[1]
        r_int = gauss_int.children[2]
        @test electric_field(r_int) / electric_field(gauss_int) ≈ 0.2 atol=1e-8
        
        # Consistency of Polarization and Amplitude for AstigmaticGaussianBeamlet
        p_ray_ext = PolarizedRay([0.0, -5e-3, 0.0], [0.0, 1.0, 0.0], 1000e-9, [1.0, 0.0, 0.0])
        p_beam = Beam(p_ray_ext)
        solve_system!(system_ext, p_beam; retrace=false)
        
        @test length(p_beam.children) == 2
        p_r = first(rays(p_beam.children[2]))
        
        agb_ext = AstigmaticGaussianBeamlet([0.0, -5e-3, 0.0], [0.0, 1.0, 0.0], 1000e-9, 1e-3; E0=[1.0, 0.0, 0.0], support=[0.0, 0.0, 1.0])
        solve_system!(system_ext, agb_ext; retrace=false)
        
        @test length(agb_ext.children) == 2
        agb_r = agb_ext.children[2]
        agb_r_chief = first(rays(agb_r.c))
        @test agb_r_chief.E0 ≈ p_r.E0 atol=1e-8
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
end
