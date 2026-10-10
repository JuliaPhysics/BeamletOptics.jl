module TestGaussianBeamlet

using BeamletOptics
using Test

const BMO = BeamletOptics

const mm = 1e-3

@testset "Gaussian beamlet" begin
    @testset "Testing type definitions" begin
        @test isdefined(BMO, :GaussianBeamlet)
    end

    @testset "Testing analytical equations" begin
        λ = 500e-9
        w0 = 1mm
        M2 = 1
        zR = BMO.rayleigh_range(λ, w0, M2)
        # Test Rayleigh range and div. angle against Paschotta (https://www.rp-photonics.com/gaussian_beams.html)
        @test isapprox(zR, 6.28, atol = 1e-2)
        @test isapprox(BMO.beam_waist(zR, w0, zR), sqrt(2) * w0)
        @test isapprox(BMO.gouy_phase(zR, zR), -π / 4)
        @test isapprox(BMO.wavefront_curvature(zR, zR), 1 / (2 * zR))
        @test isapprox(BMO.divergence_angle(λ, w0, M2), 159e-6, atol = 1e-6)
    end

    @testset "Testing parameter correctness" begin
        # Gauss beam parameters
        y = -5:0.01:5       # m
        λ_1 = 500e-9        # m
        λ_2 = 1000e-9       # m
        P0 = 1              # W
        r = 0               # m
        w0_1 = 1mm         # m
        w0_2 = 2mm         # m
        M2_1 = 1mm         # m
        M2_2 = 2mm         # m
        E0_1 = BMO.electric_field(2 * P0 / (π * w0_1^2))
        E0_2 = BMO.electric_field(2 * P0 / (π * w0_2^2))
        gauss_1 = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0],
            λ_1,
            w0_1,
            M2 = M2_1,
            P0 = P0)
        gauss_2 = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0],
            λ_2,
            w0_2,
            M2 = M2_2,
            P0 = P0)
        # Calculate analytical values
        zr_1 = BMO.rayleigh_range(λ_1, w0_1, M2_1)
        zr_2 = BMO.rayleigh_range(λ_2, w0_2, M2_2)
        wa_1 = BMO.beam_waist.(y, w0_1, zr_1)
        wa_2 = BMO.beam_waist.(y, w0_2, zr_2)
        Ra_1 = BMO.wavefront_curvature.(y, zr_1)
        Ra_2 = BMO.wavefront_curvature.(y, zr_2)
        ψa_1 = BMO.gouy_phase.(y, zr_1)
        ψa_2 = BMO.gouy_phase.(y, zr_2)
        Ea_1 = BMO.electric_field.(r, y, E0_1, w0_1, λ_1, M2_1)
        Ea_2 = BMO.electric_field.(r, y, E0_2, w0_2, λ_2, M2_2)
        # Calculate numerical values
        wn_1, Rn_1, ψn_1, w0n_1 = BMO.gauss_parameters(gauss_1, y)
        wn_2, Rn_2, ψn_2, w0n_2 = BMO.gauss_parameters(gauss_2, y)
        En_1 = [BMO.electric_field(gauss_1, r, yi) for yi in y]
        En_2 = [BMO.electric_field(gauss_2, r, yi) for yi in y]
        # Compare beam diameter within 0.1 nm
        @test all(isapprox.(wa_1, wn_1, atol = 1e-10))
        @test all(isapprox.(wa_2, wn_2, atol = 1e-10))
        # Compare wavefront curvature
        @test all(isapprox.(Ra_1, Rn_1, atol = 5e-9))
        @test all(isapprox.(Ra_2, Rn_2, atol = 5e-9))
        # Compare Gouy phase within 0.1 μrad
        @test all(isapprox.(ψa_1, ψn_1, atol = 1e-7))
        @test all(isapprox.(ψa_2, ψn_2, atol = 1e-7))
        # Compare calculated waist radius
        @test all(isapprox.(w0_1, w0n_1))
        @test all(isapprox.(w0_2, w0n_2))
        # Compare calculate electric field at r
        @test all(isapprox.(Ea_1, En_1, atol = 1e-8))
        @test all(isapprox.(Ea_2, En_2, atol = 1e-7))
        # Compare beam power with original value
        @test P0 ≈ BMO.optical_power(gauss_1)
        @test P0 ≈ BMO.optical_power(gauss_2)
    end

    @testset "Support vector" begin
        pos = [0.0, 0, 0]
        dir = [0.0, 1, 0]
        λ = 1000e-9
        w0 = 1mm
        waist_start(gauss) = position(first(BMO.rays(gauss.waist)))
        ref = GaussianBeamlet(pos, dir, λ, w0; support = [1.0, 0, 0])
        @test waist_start(ref) ≈ [w0, 0, 0]
        # the length of an orthogonal support vector does not matter
        for scale in (1e-9, 5, 1e9)
            gauss = GaussianBeamlet(pos, dir, λ, w0; support = scale .* [1.0, 0, 0])
            @test waist_start(gauss) ≈ waist_start(ref)
            @test BMO.direction(first(BMO.rays(gauss.divergence))) ≈
                  BMO.direction(first(BMO.rays(ref.divergence)))
        end
        # nor does the length of the direction
        @test waist_start(GaussianBeamlet(pos, 3 .* dir, λ, w0; support = [2.0, 0, 0])) ≈ waist_start(ref)
        # deviations within the threshold pass, e.g. a support vector from a rotation
        @test GaussianBeamlet(pos, dir, λ, w0; support = [1.0, 1e-12, 0]) isa GaussianBeamlet
        # not orthogonal: the waist radius would be w0 * sin of the angle to the direction
        @test_throws ArgumentError GaussianBeamlet(pos, dir, λ, w0; support = [1.0, 1, 0])
        @test_throws ArgumentError GaussianBeamlet(pos, dir, λ, w0; support = [1.0, 1e-6, 0])
        @test_throws ArgumentError GaussianBeamlet(pos, dir, λ, w0; support = 1e9 .* [1.0, 1e-6, 0])
        # parallel and anti-parallel
        @test_throws ArgumentError GaussianBeamlet(pos, dir, λ, w0; support = [0.0, 1, 0])
        @test_throws ArgumentError GaussianBeamlet(pos, dir, λ, w0; support = [0.0, -2, 0])
        # no direction
        @test_throws ArgumentError GaussianBeamlet(pos, dir, λ, w0; support = [0.0, 0, 0])
        # the default is orthogonal for any direction
        for d in ([1.0, 0, 0], [0.0, 0, -1], [1.0, 2, 3], [-1e-3, 1, 1e-3])
            gauss = GaussianBeamlet(pos, d, λ, w0)
            s = waist_start(gauss) - pos
            @test isapprox(BMO.dot(s, BMO.normalize(d)), 0; atol = 1e-15)
            @test BMO.norm(s) ≈ w0
        end
    end

    @testset "Testing propagation correctness" begin
        # Analytical result using complex q factor
        q_ana(q0::Complex, M::Matrix) = (M[1] * q0 + M[3]) / (M[2] * q0 + M[4])
        R_ana(q::Complex) = real(1 / q)
        w_ana(q::Complex, λ, n = 1) = sqrt(-λ / (π * n * imag(1 / q)))
        propagate_ABCD(d) = [1 d; 0 1]
        lensmaker_ABCD(f) = [1 0; -1/f 1]
        # Beam parameters
        λ = 1000e-9
        w0 = 1mm
        M2 = 1
        zr = BMO.rayleigh_range(λ, w0, M2)
        # Lens parameters
        R1 = 1
        R2 = 1
        lens_y_location = 0.1
        nl = 1.5
        f = BMO.lensmakers_eq(R1, -R2, nl)
        # Stuff
        dy = 0.001
        ys = 0:dy:1.5
        w_analytical = Vector{Float64}(undef, length(ys))
        R_analytical = Vector{Float64}(undef, length(ys))
        # Propagate using ABCD formalism
        q0 = 0 + zr * im
        for i in 1:length(ys)
            w_analytical[i] = w_ana(q0, λ)
            R_analytical[i] = R_ana(q0)
            # catch first lens
            if i * dy == lens_y_location
                q0 = q_ana(q0, lensmaker_ABCD(f))
                continue
            end
            q0 = q_ana(q0, propagate_ABCD(dy))
        end

        # Numerical result
        tl = BMO.ThinLensSDF(R1, R2, 0.025)
        lens = Lens(tl, x -> nl)
        system = System(lens)
        translate3d!(lens, [0, lens_y_location, 0])
        # Create and solve beam, calculate beam parameters
        gauss = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0],
            λ,
            w0,
            support = [1, 0, 0],
            M2 = 1)
        solve_system!(system, gauss)
        w_numerical, R_numerical, ψ_numerical, w0_numerical = BMO.gauss_parameters(
            gauss,
            ys)
        # Compare beam radius to within 1 μm
        @test all(isapprox.(w_analytical, w_numerical, atol = 1e-6))
        # Compare if radius of curvature agreement to within 1 cm above 95% and any NaN
        temp = isapprox.(R_analytical, R_numerical, atol = 1e-2)
        @test sum(temp) / length(temp) > 0.95
        @test !any(isnan.(R_numerical))
        # Compare if Gouy phase zero at waists
        w0, i = findmin(w_analytical)
        @test isapprox(ψ_numerical[1], 0, atol = 1mm)
        @test isapprox(ψ_numerical[i], 0, atol = 1mm)
        # Compare calculated waist after lens
        @test isapprox(w0_numerical[i], w0, atol = 1e-7)

        @testset "Testing isparaxial and istilted" begin
            # Before lens rotation
            @test BMO.istilted(system, gauss) == false
            @test BMO.isparaxial(system, gauss) == true
            # Tilt lens, test again with 30° threshold for paraxial approx.
            zrotate3d!(lens, deg2rad(45))
            solve_system!(system, gauss)
            @test BMO.istilted(system, gauss) == true
            @test BMO.isparaxial(system, gauss, deg2rad(30)) == false
        end
    end
end

@testset "Field phase in and in front of a medium" begin
    # https://github.com/JuliaPhysics/BeamletOptics.jl/issues/126
    λ = 1e-6
    n = 1.5
    # glass from y = 50.0 mm to about 50.5 mm
    lens = ThinLens(50mm, 50mm, 10mm, n)
    translate3d!(lens, [0, 50mm, 0])
    system = System([lens])
    Δφ(a, b) = rad2deg(angle(a / b))
    args = ([0.0, 0, 0], [0.0, 1, 0], λ, 1mm)
    free = GaussianBeamlet(args...; support = [1.0, 0, 0])
    gb = GaussianBeamlet(args...; support = [1.0, 0, 0])
    agb = AstigmaticGaussianBeamlet(args...; support = [1.0, 0, 0])
    solve_system!(system, gb)
    solve_system!(system, agb)
    @test BMO.refractive_index(gb, 2) == n
    @test length(BMO.rays(gb.chief)[2]) > 0.4mm
    # the lens does not change the field in front of it
    @test electric_field(gb, 0.0, 10mm) ≈ electric_field(free, 0.0, 10mm)
    @test electric_field(gb, 0.3mm, 40mm) ≈ electric_field(free, 0.3mm, 40mm)
    # inside the glass the phase advances with n k0 along the axis ...
    Δz = 100.2e-6
    @test Δφ(electric_field(gb, 0.0, 50.1mm + Δz), electric_field(gb, 0.0, 50.1mm)) ≈ mod(360 * n * Δz / λ, 360) atol = 0.01
    # ... and across the beam, as in the astigmatic model
    r = 0.1mm
    across_gb = Δφ(electric_field(gb, r, 50.25mm), electric_field(gb, 0.0, 50.25mm))
    across_agb = Δφ(electric_field(agb, [r, 0, 0], 50.25mm), electric_field(agb, [0.0, 0, 0], 50.25mm))
    @test across_gb ≈ across_agb atol = 0.01
    @test abs(across_gb) > 15
    # a hit of the segment in the glass gives the field of the beam
    hit = BMO.GaussianBeamletHit(gb, 2)
    p = BMO.Point3(r, 50.25mm, 0.0)
    @test BMO.beamlet_hit_field(hit, p) ≈ electric_field(gb, r, 50.25mm)
end

@testset "Gouy phase is collected along the beam" begin
    # https://github.com/JuliaPhysics/BeamletOptics.jl/issues/125
    λ = 1e-6
    n = 1.5
    w0 = 1mm
    R = 50mm
    k = 2π / λ
    lens(y) = (l = ThinLens(R, R, 10mm, n); translate3d!(l, [0, y, 0]); l)   # f = 50 mm
    gb() = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; support = [1.0, 0, 0])
    agb() = AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; support = [1.0, 0, 0])
    # field on the axis of the detector without the phase of the optical path [deg]
    function axis_field(system, pd, beam)
        empty!(pd)
        solve_system!(system, beam)
        _, _, E = electric_field(pd; n = 1, x_min = 0, x_max = 0, z_min = 0, z_max = 0)
        return E[1] * cis(-k * BMO.optical_path_length(beam))
    end
    gouy(args...) = rad2deg(angle(axis_field(args...)))
    detector(y) = (pd = Detector(20mm); translate3d!(pd, [0, y, 0]); pd)

    @testset "Behind one lens" begin
        for y in (70mm, 95mm, 110mm, 200mm)
            pd = detector(y)
            system = System([lens(50mm), pd])
            beam = gb()
            ψ = gouy(system, pd, beam)
            # ABCD law for the traced thick lens: Gouy phase = -arg(A + B/q) with q = -i z_R at the waist
            l1, l2, l3 = length.(BMO.rays(beam.chief))
            M = [1 l3; 0 1] * [1 0; (n - 1) / -R n] * [1 l2; 0 1] * [1 0; (1 - n) / (n * R) 1 / n] * [1 l1; 0 1]
            ψ_abcd = rad2deg(-angle(M[1, 1] + M[1, 2] / (-im * π * w0^2 / λ)))
            @test rad2deg(angle(cis(deg2rad(ψ - ψ_abcd)))) ≈ 0 atol = 0.1
        end
        # the astigmatic model agrees away from the focus, where its parabasal rays see other aberrations
        pd = detector(70mm)
        system = System([lens(50mm), pd])
        @test gouy(system, pd, gb()) ≈ gouy(system, pd, agb()) atol = 0.1
        # the phase is continuous through the lens
        beam = gb()
        empty!(pd)
        solve_system!(system, beam)
        Δφ(a, b) = rad2deg(angle(a / b))
        step = Δφ(electric_field(beam, 0.0, 50.6mm), electric_field(beam, 0.0, 49.9mm))
        opl = BMO.optical_path_length(beam) - length(beam) + 0.7mm
        @test rad2deg(angle(cis(deg2rad(step) - k * opl))) ≈ 0 atol = 0.5
    end

    @testset "Behind two lenses" begin
        pd = detector(250mm)
        system = System([lens(50mm), lens(150mm), pd])
        @test gouy(system, pd, gb()) ≈ gouy(system, pd, agb()) atol = 0.1
    end

    @testset "Without a lens" begin
        pd = detector(70mm)
        system = System([pd])
        @test gouy(system, pd, gb()) ≈ rad2deg(-atan(70mm / (π * w0^2 / λ))) atol = 1e-6
        m = SquarePlanoMirror2D(25mm)
        translate3d!(m, [0, 100mm, 0])
        zrotate3d!(m, π / 4)
        pd = detector(0.0)
        translate3d!(pd, [50mm, 100mm, 0])
        zrotate3d!(pd, -π / 2)
        beam = gb()
        ψ = gouy(System([m, pd]), pd, beam)
        @test length(BMO.rays(beam.chief)) == 2
        @test ψ ≈ rad2deg(-atan(150mm / (π * w0^2 / λ))) atol = 1e-6
    end

    @testset "From a parent to its children" begin
        pd_t = detector(200mm)
        pd_r = detector(0.0)
        translate3d!(pd_r, [100mm, 100mm, 0])
        zrotate3d!(pd_r, -π / 2)
        bs = ThinBeamsplitter(25mm)
        translate3d!(bs, [0, 100mm, 0])
        zrotate3d!(bs, π / 4)
        direct = axis_field(System([lens(50mm), pd_t]), pd_t, gb())
        system = System([lens(50mm), bs, pd_t, pd_r])
        empty!(pd_r)
        beam = gb()
        transmitted = axis_field(system, pd_t, beam)
        @test BMO.hit_count(pd_t) == 1 && BMO.hit_count(pd_r) == 1
        @test transmitted ≈ direct / sqrt(2) rtol = 1e-6
        # same distance behind the splitter: the reflected field differs by its sign only
        _, _, E_r = electric_field(pd_r; n = 1, x_min = 0, x_max = 0, z_min = 0, z_max = 0)
        reflected = E_r[1] * cis(-k * BMO.optical_path_length(beam.children[2]))
        @test abs(reflected) ≈ abs(transmitted) rtol = 1e-6
        @test abs(sind(rad2deg(angle(reflected / transmitted)))) < 1e-6
    end
end

end # MODULE