module TestPolarizedRays

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics

@testset "Polarized rays" begin
    @testset "Polarization transforms" begin
        # Reflection matrix and lin. x-pol
        J = BMO.SPBasis(-1, 0, 0, 1)
        E0 = [1, 0, 0]

        @testset "90° reflection" begin
            in_dir = [0, 0, 1]
            out_dir = [1, 0, 0]
            nml = normalize([1, 0, -1])
            P90 = BMO._calculate_global_E0(in_dir, out_dir, nml, J)
            @test P90 * E0 ≈ [0, 0, -1]
            @test P90 * in_dir ≈ out_dir
        end

        @testset "0° reflection" begin
            in_dir = [0, 0, 1]
            out_dir = [0, 0, -1]
            nml = normalize([0, 0, -1])
            P00 = BMO._calculate_global_E0(in_dir, out_dir, nml, J)
            @test P00 * E0 ≈ [-1, 0, 0]
            @test P00 * in_dir ≈ out_dir
        end
    end

    @testset "Test error messages" begin
        # Test dir. error msg
        @test_throws ErrorException PolarizedRay(zeros(3), zeros(3))
        @test_throws ErrorException PolarizedRay(zeros(3), ones(3)*eps())
        # Test polarization orthogonal error msg
        @test_throws ErrorException PolarizedRay(zeros(3), [0,1,0], 1e-6, [0,1,1])
    end

    @testset "Mirror reflections" begin
        # Setup system as in https://opg.optica.org/ao/fulltext.cfm?uri=ao-50-18-2855&id=218813
        m1 = SquarePlanoMirror2D(1.0)
        m2 = SquarePlanoMirror2D(1.0)
        m3 = SquarePlanoMirror2D(1.0)
        translate3d!(m2, [2, 0, 0])
        translate3d!(m3, [2, 2, 0])
        zrotate3d!(m1, deg2rad(-90))
        yrotate3d!(m1, deg2rad(45))
        zrotate3d!(m2, deg2rad(45))
        xrotate3d!(m3, deg2rad(135))

        system = StaticSystem([m1, m2, m3])

        I0_1 = 1
        I0_2 = 5
        lin_x_pol = [I0_1, 0, 0]
        lin_y_pol = [0, I0_2, 0]

        # Beam of polarized rays
        ray = PolarizedRay([0.0, 0, -2], [0, 0, 1], 1000e-9, lin_x_pol)
        beam = Beam(ray)

        @testset "x-Polarization" begin
            BMO.polarization!(ray, lin_x_pol)
            # test tracing
            solve_system!(system, beam)
            @test BMO.polarization(beam.rays[1]) ≈ lin_x_pol
            @test BMO.polarization(beam.rays[2]) ≈ [0, 0, -I0_1]
            @test BMO.polarization(beam.rays[3]) ≈ [0, 0, I0_1]
            @test BMO.polarization(beam.rays[4]) ≈ [0, -I0_1, 0]
            @test length(beam) ≈ 6.0
        end

        @testset "y-Polarization" begin
            BMO.polarization!(ray, lin_y_pol)
            translate3d!(m3, [0, 2, 0])
            # test retracing
            solve_system!(system, beam)
            @test BMO.polarization(beam.rays[1]) ≈ lin_y_pol
            @test BMO.polarization(beam.rays[2]) ≈ [0, -I0_2, 0]
            @test BMO.polarization(beam.rays[3]) ≈ [I0_2, 0, 0]
            @test BMO.polarization(beam.rays[4]) ≈ [-I0_2, 0, 0]
            @test length(beam) ≈ 8.0
        end
    end

    @testset "Brewster windows" begin
        brewster_angle(n) = atan(n)
        # Testcase based on 5 successive Brewster windows
        n = 1.5
        θb = brewster_angle(n)
        d = 0.1
        # Calculate transmission efficiency
        rs, rp, ts, tp = BMO.fresnel_coefficients(θb, n)
        Ts = 1 - abs2(rs)
        Tp = 1 - abs2(rp)
        # Setup testcase
        s1 = BMO.CuboidMesh(1.0, d, 1.0)
        s2 = BMO.CuboidMesh(1.0, d, 1.0)
        s3 = BMO.CuboidMesh(1.0, d, 1.0)
        s4 = BMO.CuboidMesh(1.0, d, 1.0)
        s5 = BMO.CuboidMesh(1.0, d, 1.0)
        l1 = Lens(s1, x -> n)
        l2 = Lens(s2, x -> n)
        l3 = Lens(s3, x -> n)
        l4 = Lens(s4, x -> n)
        l5 = Lens(s5, x -> n)
        translate3d!.([l1, l2, l3, l4, l5], Ref([-0.5, -d / 2, -0.5]))
        BMO.set_new_origin3d!.(BMO.shape.([l1, l2, l3, l4, l5]))
        translate3d!(l2, [0, 0.5, -1d / 2])
        translate3d!(l3, [0, 1.0, -2d / 2])
        translate3d!(l4, [0, 1.5, -3d / 2])
        translate3d!(l5, [0, 2.0, -4d / 2])
        xrotate3d!.([l1, l2, l3, l4, l5], -θb)
        # Solve system of s- and p-polarized beams
        system = StaticSystem([l1, l2, l3, l4, l5])
        x_pol_ray = PolarizedRay(
            [-0.1, -1, 0], [0, 1.0, 0], 1000e-9, [BMO.electric_field(1), 0, 0])
        z_pol_ray = PolarizedRay(
            [+0.1, -1, 0], [0, 1.0, 0], 1000e-9, [0, 0, BMO.electric_field(1)])
        s_beam = Beam(x_pol_ray)
        p_beam = Beam(z_pol_ray)
        solve_system!(system, s_beam)
        solve_system!(system, p_beam)
        # Since system is non-focussing, calculate pseudo-intensity
        pseudo_Is = abs2(BMO.polarization(last(BMO.rays(s_beam)))[1]) /
                    (2 * BMO.Z_vacuum)
        pseudo_Ip = abs2(BMO.polarization(last(BMO.rays(p_beam)))[3]) /
                    (2 * BMO.Z_vacuum)
        # Test against m interfaces
        m = length(system.objects) * 2
        @test pseudo_Is ≈ Ts^m
        @test pseudo_Ip ≈ Tp^m
    end

    @testset "Fresnel rhomb" begin
        # Create Fresnel rhomb with n=1.5 and θ=53.3° for quarter-wave plate effect
        n = 1.5
        s1 = BMO.CuboidMesh(0.5, 1.25, 0.5, deg2rad(53.3))
        l1 = Lens(s1, x -> n)
        translate3d!(l1, [-0.25, 0, -0.25])
        BMO.set_new_origin3d!(s1)
        # Rotate prism to obtain 45° beam input polarization
        yrotate3d!(l1, deg2rad(135))
        # Solve system
        system = StaticSystem([l1])
        ray = PolarizedRay(
            [0, -1, 0], [0, 1.0, 0], 1000e-9, [0, 0, BMO.electric_field(1)])
        beam = Beam(ray)
        solve_system!(system, beam)
        # Assumes propagation along the y-axis after rhomb, calculate polarization state
        Ex = getindex.(BMO.polarization.(beam.rays), 1)
        Ey = getindex.(BMO.polarization.(beam.rays), 2)
        Ez = getindex.(BMO.polarization.(beam.rays), 3)
        # Test for circular polarization and Ey error
        phi = angle(last(Ez)) - angle(last(Ex))
        @test phi ≈ π / 2
        @test abs(last(Ey)) < 2e-14
    end

    @testset "Detector electric_field (coherent vector sum)" begin
        # Two plane waves crossing at the origin at ±θ from the y-axis. The OPL is chosen
        # so that each ray's phasor at the crossing point is 1+0im. Along the detector x-axis
        # they form fringes of period Λ whose contrast equals |E₊⋅E₋| / (|E₊||E₋|), which
        # a scalar sum cannot reproduce.
        λ = 1e-6
        θ = deg2rad(30)
        Λ = λ / (2 * sin(θ))
        d_plus = [sin(θ), cos(θ), 0]
        d_minus = [-sin(θ), cos(θ), 0]
        s_pol = [0, 0, 1.0]
        p_plus = [cos(θ), -sin(θ), 0]
        p_minus = [cos(θ), sin(θ), 0]

        function manual_hit(dir, E0; focus = zeros(3))
            L = λ   # k*L = 2π -> cis(k*opl) = 1 at the focus
            ray = PolarizedRay(focus .- L .* dir, dir, λ, E0)
            BMO.intersection!(ray, BMO.Intersection(L, BMO.Point3(-dir)))
            return BMO.PolarizedRayHit(ray, BMO.optical_path_length(ray))
        end

        function two_ray_detector(E_plus, E_minus)
            pd = Detector(1.0)
            push!(pd, manual_hit(d_plus, E_plus))
            push!(pd, manual_hit(d_minus, E_minus))
            return pd
        end

        # index 101 is x = 0 (crossing point), index 51 is x = -Λ/2 (half a fringe away)
        x_scan(pd) = electric_field(pd; n = 201, x_min = -Λ, x_max = Λ, z_min = 0.0, z_max = 0.0)
        visibility(I) = (maximum(I) - minimum(I)) / (maximum(I) + minimum(I))

        @testset "s-polarization: full-contrast fringes" begin
            xs, zs, E = x_scan(two_ray_detector(s_pol, s_pol))
            @test eltype(E) <: BMO.Point3{<:Complex}
            I = intensity.(E[:, 1])

            @test sum(abs2, E[101, 1]) ≈ 4
            @test I[51] / I[101] < 1e-12
            @test visibility(I) ≈ 1
        end

        @testset "p-polarization: fringe contrast |cos 2θ|" begin
            xs, zs, E = x_scan(two_ray_detector(p_plus, p_minus))
            E_focus = E[101, 1]
            I = intensity.(E[:, 1])

            # in-plane (y) components cancel by symmetry
            @test abs(E_focus[2]) < 1e-12
            @test sum(abs2, E_focus) ≈ 4 * cos(θ)^2
            @test intensity(E_focus) ≈ sum(abs2, E_focus) / (2 * BMO.Z_vacuum)
            @test I[51] / I[101] ≈ (1 - cos(2θ)) / (1 + cos(2θ))
            @test visibility(I) ≈ abs(cos(2θ))
        end

        @testset "orthogonal polarizations do not interfere" begin
            pd = two_ray_detector(s_pol, p_minus)
            xs, zs, I = intensity(pd; n = 101, x_min = -Λ, x_max = Λ, z_min = -Λ, z_max = Λ)

            @test maximum(I) - minimum(I) < 1e-12 * maximum(I)
            @test 2 * BMO.Z_vacuum * maximum(I) ≈ 2
        end

        @testset "tilted detector: no projection factor" begin
            # Documents that, unlike RayHit, the vector field is not weighted by the incidence angle
            tilt = deg2rad(40)
            probe(pd) = intensity(pd; n = 1, x_min = 0.0, x_max = 0.0, z_min = 0.0, z_max = 0.0)[3][1, 1]

            pd_pol = Detector(1.0)
            zrotate3d!(pd_pol, tilt)
            solve_system!(StaticSystem([pd_pol]), Beam([0, -1.0, 0], [0, 1.0, 0], λ, s_pol))

            pd_ray = Detector(1.0)
            zrotate3d!(pd_ray, tilt)
            solve_system!(StaticSystem([pd_ray]), Beam([0, -1.0, 0], [0, 1.0, 0], λ))

            @test 2 * BMO.Z_vacuum * probe(pd_pol) ≈ 1
            @test 2 * BMO.Z_vacuum * probe(pd_ray) ≈ cos(tilt)^2
        end

        @testset "aplanatic high-NA focus: vectorial PSF" begin
            # Converging x-polarized ray fan on a sphere around the focus (optical axis +y)
            # with the Richards–Wolf field and apodization of an aplanatic lens.
            function aplanatic_detector(NA; N = 60, pol = [1.0, 0, 0])
                pd = Detector(1.0)
                for i in 1:N, j in 1:(4N)
                    α = asin(NA) * (i - 0.5) / N
                    φ = 2π * (j - 0.5) / (4N)
                    dir = [sin(α) * cos(φ), cos(α), sin(α) * sin(φ)]
                    e_r = [cos(φ), 0, sin(φ)]
                    e_φ = [-sin(φ), 0, cos(φ)]
                    e_ρ = [cos(α) * cos(φ), -sin(α), cos(α) * sin(φ)]
                    E0 = sqrt(cos(α)) * sin(α) .* (dot(pol, e_r) .* e_ρ .+ dot(pol, e_φ) .* e_φ)
                    push!(pd, manual_hit(dir, E0))
                end
                return pd
            end

            function vector_psf(NA; n = 151)
                R = 1.5λ / NA
                xs, zs, E = electric_field(aplanatic_detector(NA); n, x_min = -R, x_max = R, z_min = -R, z_max = R)
                c = (n + 1) ÷ 2
                I = intensity.(E)
                Ex2 = map(e -> abs2(e[1]), E)
                Ey2 = map(e -> abs2(e[2]), E)
                fwhm_count(v) = count(>(maximum(v) / 2), v)
                return (
                    peak = argmax(I) == CartesianIndex(c, c),
                    fwhm_ratio = fwhm_count(I[:, c]) / fwhm_count(I[c, :]),
                    Ey_center = Ey2[c, c] / maximum(Ey2),
                    Ey_ratio = maximum(Ey2) / maximum(Ex2),
                    Ey_lobe_row = argmax(Ey2)[2] == c
                )
            end

            low = vector_psf(0.1)
            @test low.peak
            @test low.fwhm_ratio ≈ 1 atol = 0.1
            @test low.Ey_ratio < 0.01

            high = vector_psf(0.9)
            @test high.peak
            # PSF is stretched along the input polarization (x)
            @test high.fwhm_ratio > 1.2
            # longitudinal component vanishes on axis, but forms strong lobes along x
            @test high.Ey_center < 1e-12
            @test high.Ey_ratio > 0.1
            @test high.Ey_lobe_row
        end
    end

    @testset "Polarized point spread function (vector coherent sum)" begin
        # Same optical setup as the scalar Airy-disc test (see TestDetector.jl),
        # but tracing a linearly polarized ray fan built from the same
        # deterministic UniformDiscSource sampling, so the geometry matches exactly.
        mm = 1e-3
        R1 = 100mm
        R2 = Inf
        l = 1mm
        d = 25.4mm
        n = 1.5
        λ = 1e-6
        D = 15mm
        num_rays = 1000

        lens = SphericalLens(R1, R2, l, d, x -> n)

        x_shift = y_shift = -2mm
        psfd = Detector(10e-3)
        translate3d!(psfd, [x_shift, 200e-3 + 0.13e-3, y_shift])

        cs = UniformDiscSource([0, -10e-3, 0], [0, 1, 0], D, λ; num_rays)
        E0 = [1.0, 0, 0]
        pol_beams = [Beam(position(first(rays(b))), direction(first(rays(b))), λ, E0)
                     for b in BMO.beams(cs)]
        pcs = BMO.CollimatedSource(pol_beams, D)

        sys = System([lens, psfd])
        solve_system!(sys, pcs)

        x, y, I_num = intensity(psfd; n = 500, crop_factor = 5, center = MinMax())

        # walk from the peak to the first local minimum through the centre column
        ix_ctr, jx_ctr = Tuple(argmax(I_num))
        col = I_num[:, jx_ctr]
        i_min = ix_ctr
        while i_min < length(col) && col[i_min + 1] < col[i_min]
            i_min += 1
        end

        # compare relative to the peak, since the absolute offset (x_shift) dwarfs the Airy radius
        airy_radius = 1.22 * λ * 200e-3 / D
        @test x[i_min] - x[ix_ctr] ≈ airy_radius rtol = 2e-2

        # vector field result feeds the scalar intensity/power pipeline
        @test eltype(electric_field(psfd; n = 2)[3]) <: BMO.Point3{<:Complex}
        @test BMO.optical_power(psfd; n = 100, crop_factor = 5, center = MinMax()) > 0
    end
end

end # MODULE