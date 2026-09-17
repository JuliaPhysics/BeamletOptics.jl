module TestPolarizers

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics

const mm = 1e-3

@testset "Polarizing components" begin
    # Test fcts.
    pseudo_I(v) = norm(v)^2
    pseudo_I(ray::BMO.PolarizedRay) = pseudo_I(BMO.polarization(ray))
    malus_law(θ) = cosd(θ)^2
    function generate_linpol_angles(thetas)
        angles = zeros(length(thetas))
        # consider flip in angle quadrant above 90°  
        for i in eachindex(thetas)
            theta = thetas[i] - 1
            if theta > 90 && theta <= 270
                angles[i] = 180
            end
        end
        return angles
    end

    @testset "Jones polarization matrices" begin
        @test isdefined(BMO, :AbstractJonesMatrix)
        @test isdefined(BMO, :LocalJonesBasis)
        @test isdefined(BMO, :GlobalJonesBasis)
        # test matrix operations
        J = BMO.XYBasis(1, 2, 3, 4)
        R = [1 2 0; 3 4 0; 0 0 1]
        @test J == R
        @test J*2 == R*2
        @test transpose(J) == transpose(R)
        @test inv(J) ≈ inv(R)
    end

    @testset "Polarizing filter" begin
        filter = PolarizationFilter(5mm)
        system = System([filter])
        # rotation angle steps, intensity and angle comparison
        thetas = 1:10:360
        i_n = zeros(length(thetas))
        i_a = similar(i_n)
        angles_n = similar(i_n)
        angles_a = generate_linpol_angles(thetas)

        @testset "Pol. filter normal incidence - rotation around y-axis" begin 
            Rfil = orientation(filter)
            ray_dir = Rfil[:,2]
            ray_pos = position(filter) - 10mm * ray_dir
            local_x = Rfil[:,1]
            local_z = Rfil[:,3]
            pol_vec = local_x
            lambda = 1e-6
            beam = Beam(ray_pos, ray_dir, lambda, pol_vec)            
            # step ref. vector by angle increment 
            RotMat = BMO.rotate3d(ray_dir, deg2rad(step(thetas)))
            for i in eachindex(thetas)
                solve_system!(system, beam)
                E1 = BMO.polarization(BMO.rays(beam)[2])
                i_n[i] = pseudo_I(E1)
                i_a[i] = malus_law(thetas[i]-1)
                rotate3d!(filter, ray_dir, deg2rad(step(thetas)))
                # test polarization state
                @test all(BMO.islinear.(beam.rays))
                # test projected polarization direction
                v1 = real(E1)
                angles_n[i] = rad2deg(BMO.angle3d(v1, pol_vec))
                # update ref vector
                pol_vec = RotMat * pol_vec
            end    
            @test angles_n ≈ angles_a
        end

        # Move and rotate filter
        translate3d!(filter, [0,10mm,0])
        xrotate3d!(filter, deg2rad(45))
        zrotate3d!(filter, deg2rad(30))

        @testset "Pol. filter tilted incidence - rotation around ray optical axis" begin
            Rfil = orientation(filter)
            ray_dir = Rfil[:,2]
            ray_pos = position(filter) - 10mm * ray_dir
            local_x = Rfil[:,1]
            local_z = Rfil[:,3]
            pol_vec = local_x
            lambda = 1e-6
            beam = Beam(ray_pos, ray_dir, lambda, pol_vec)
            # tilt filter
            θ_tilt = 45
            rotate3d!(filter, local_x, deg2rad(θ_tilt))
            i_n_tilt = similar(i_n)
            for i in eachindex(thetas)
                solve_system!(system, beam)
                i_n_tilt[i] = pseudo_I(beam.rays[2].E0)
                rotate3d!(filter, ray_dir, deg2rad(step(thetas)))
                @test all(BMO.islinear.(beam.rays))
            end
            # only applicable for θ_tilt = 45°
            i_a_tilt = malus_law.(thetas .- 1) .* 0.75 .+ 0.25
            @test i_n_tilt ≈ i_a_tilt
        end
    end

    @testset "Polarizing filter - astigmatic Gaussian beamlet" begin
        # Default filter at origin: transmits along global x, blocks global z
        filter = PolarizationFilter(5mm)
        system = System([filter])
        λ = 633e-9
        w0 = 0.5mm
        pos = [0, -10mm, 0]
        dir = [0, 1, 0]

        @testset "Transmitted polarization" begin
            ref = AstigmaticGaussianBeamlet(pos, dir, λ, w0; E0 = [1, 0, 0], support = [1, 0, 0])
            agb = AstigmaticGaussianBeamlet(pos, dir, λ, w0; E0 = [1, 0, 0], support = [1, 0, 0])
            @test_nowarn solve_system!(system, agb)
            # All component beams pass the filter in sync
            lengths = map(b -> length(BMO.rays(b)), BMO._component_beams(agb))
            @test all(==(2), lengths)
            # Beam radii behind the ideal filter match free-space propagation
            z = 25mm
            w1, w2 = BMO.gauss_parameters(agb, z)
            w1_ref, w2_ref = BMO.gauss_parameters(ref, z)
            @test w1 ≈ w1_ref
            @test w2 ≈ w2_ref
            # Chief polarization is unchanged
            @test BMO.polarization(BMO.rays(agb.c)[2]) ≈ BMO.polarization(BMO.rays(agb.c)[1])
        end

        @testset "Blocked polarization" begin
            agb = AstigmaticGaussianBeamlet(pos, dir, λ, w0; E0 = [0, 0, 1], support = [1, 0, 0])
            @test_nowarn solve_system!(system, agb)
            # Auxiliary beams stay in sync with the chief beam
            lengths = map(b -> length(BMO.rays(b)), BMO._component_beams(agb))
            @test all(==(lengths[1]), lengths)
            # Beamlet is terminated at the filter
            @test all(==(1), lengths)
            @test BMO.object(BMO.intersection(last(BMO.rays(agb.c)))) === filter
        end
    end

    @testset "Round polarization filter" begin
        filter = RoundPolarizationFilter(10mm)
        system = System([filter])

        @testset "Hit" begin
            beam = Beam([4mm, -10mm, 0], [0, 1, 0], 1e-6, [1.0, 0, 0])
            solve_system!(system, beam)
            @test length(BMO.rays(beam)) == 2
            @test BMO.polarization(BMO.rays(beam)[2]) ≈ BMO.polarization(BMO.rays(beam)[1])
        end

        @testset "Miss" begin
            beam = Beam([6mm, -10mm, 0], [0, 1, 0], 1e-6, [1.0, 0, 0])
            solve_system!(system, beam)
            @test length(BMO.rays(beam)) == 1
        end
    end

    @testset "Transmission axis" begin
        @testset "PolarizationFilter" begin
            pf = PolarizationFilter(5mm)
            t = transmission_axis(pf)
            @test norm(t) ≈ 1
            @test abs(dot(t, [1, 0, 0])) ≈ 1 atol = 1e-12

            for θ in 0:15:180
                pf2 = PolarizationFilter(5mm)
                rotate3d!(pf2, [0, 1, 0], deg2rad(θ))
                t2 = transmission_axis(pf2)
                expected = [cosd(θ), 0, -sind(θ)]
                @test norm(t2) ≈ 1
                @test abs(dot(t2, expected)) ≈ 1 atol = 1e-12
            end

            pf3 = PolarizationFilter(5mm)
            xrotate3d!(pf3, deg2rad(30))
            t3 = transmission_axis(pf3)
            n3 = orientation(pf3)[:, 2]
            @test norm(t3) ≈ 1
            @test abs(dot(t3, n3)) ≈ 0 atol = 1e-12
            @test abs(dot(t3, [1, 0, 0])) ≈ 1 atol = 1e-12
        end

        @testset "LinearPolarizer" begin
            n = λ -> 1.5
            lp = RoundLinearPolarizer(25.4mm, 1.6mm, 1.6mm, n)
            t = transmission_axis(lp)
            @test norm(t) ≈ 1
            @test abs(dot(t, [1, 0, 0])) ≈ 1 atol = 1e-12

            for θ in 0:15:180
                lp2 = RoundLinearPolarizer(25.4mm, 1.6mm, 1.6mm, n)
                rotate3d!(lp2, [0, 1, 0], deg2rad(θ))
                t2 = transmission_axis(lp2)
                expected = [cosd(θ), 0, -sind(θ)]
                @test norm(t2) ≈ 1
                @test abs(dot(t2, expected)) ≈ 1 atol = 1e-12
            end

            lp3 = RoundLinearPolarizer(25.4mm, 1.6mm, 1.6mm, n)
            xrotate3d!(lp3, deg2rad(30))
            t3 = transmission_axis(lp3)
            n3 = orientation(lp3)[:, 2]
            @test norm(t3) ≈ 1
            @test abs(dot(t3, n3)) ≈ 0 atol = 1e-12
            @test abs(dot(t3, [1, 0, 0])) ≈ 1 atol = 1e-12
        end

        @testset "Malus consistency" begin
            filter = PolarizationFilter(5mm)
            t = transmission_axis(filter)
            system = System([filter])
            beam = Beam([0, -10mm, 0], [0, 1, 0], 1e-6, t)
            solve_system!(system, beam)
            E = BMO.polarization(last(BMO.rays(beam)))
            @test norm(E) ≈ 1 atol = 1e-9
        end
    end

    @testset "Linear polarizer - geometry and kinematics" begin
        n = λ -> 1.5
        tf = 1.6mm
        tb = 1.0mm

        lp = RoundLinearPolarizer(25.4mm, tf, tb, n)
        @test position(lp) ≈ zeros(3)
        @test position(lp.front) ≈ [0, -tf, 0]
        @test position(lp.back) ≈ [0, 0, 0]
        @test thickness(lp) ≈ tf + tb

        v = [1mm, 2mm, 3mm]
        translate3d!(lp, v)
        rotate3d!(lp, [0, 0, 1], π / 4)
        yp = orientation(lp)[:, 2]
        @test position(lp) ≈ v
        @test thickness(lp) ≈ tf + tb
        @test position(lp.front) ≈ v - tf * yp
        @test position(lp.back) ≈ v
    end

    @testset "Linear polarizer - Malus law (forward and reverse)" begin
        n = λ -> 1.5
        tf = 1.6mm
        tb = 1.0mm
        D = 25.4mm
        λ = 1e-6

        # Reference system: single plate of thickness tf + tb, no film in between
        ref = Prism(BMO.PlanoSurfaceSDF(tf + tb, D), n)
        translate3d!(ref, [0, -tf, 0])
        system_ref = System([ref])

        thetas = 0:10:80

        @testset "Forward" begin
            for θ in thetas
                lp = RoundLinearPolarizer(D, tf, tb, n)
                rotate3d!(lp, [0, 1, 0], deg2rad(θ))
                system = System([lp])

                beam = Beam([0, -10mm, 0], [0, 1, 0], λ, [1.0, 0, 0])
                solve_system!(system, beam)
                beam_ref = Beam([0, -10mm, 0], [0, 1, 0], λ, [1.0, 0, 0])
                solve_system!(system_ref, beam_ref)

                E = BMO.polarization(last(BMO.rays(beam)))
                Eref = BMO.polarization(last(BMO.rays(beam_ref)))
                ratio = norm(E)^2 / norm(Eref)^2
                @test ratio ≈ cosd(θ)^2 rtol = 1e-9
                @test length(BMO.rays(beam)) == 4
            end
        end

        @testset "Reverse" begin
            for θ in thetas
                lp = RoundLinearPolarizer(D, tf, tb, n)
                rotate3d!(lp, [0, 1, 0], deg2rad(θ))
                system = System([lp])

                beam = Beam([0, 10mm, 0], [0, -1, 0], λ, [1.0, 0, 0])
                solve_system!(system, beam)
                beam_ref = Beam([0, 10mm, 0], [0, -1, 0], λ, [1.0, 0, 0])
                solve_system!(system_ref, beam_ref)

                E = BMO.polarization(last(BMO.rays(beam)))
                Eref = BMO.polarization(last(BMO.rays(beam_ref)))
                ratio = norm(E)^2 / norm(Eref)^2
                @test ratio ≈ cosd(θ)^2 rtol = 1e-9
                @test length(BMO.rays(beam)) == 4
            end
        end

        @testset "Blocked orientation (θ = 90°)" begin
            @testset "Forward" begin
                lp = RoundLinearPolarizer(D, tf, tb, n)
                rotate3d!(lp, [0, 1, 0], deg2rad(90))
                system = System([lp])
                beam = Beam([0, -10mm, 0], [0, 1, 0], λ, [1.0, 0, 0])
                @test_nowarn solve_system!(system, beam)
                # The film blocks the beam at the cemented interface: start -> front prism, then terminated
                @test length(BMO.rays(beam)) == 2
                @test BMO.shape(BMO.intersection(last(BMO.rays(beam)))) === BMO.shape(lp.front)
            end

            @testset "Reverse" begin
                lp = RoundLinearPolarizer(D, tf, tb, n)
                rotate3d!(lp, [0, 1, 0], deg2rad(90))
                system = System([lp])
                beam = Beam([0, 10mm, 0], [0, -1, 0], λ, [1.0, 0, 0])
                @test_nowarn solve_system!(system, beam)
                # The film blocks the beam at the cemented interface: start -> back prism, then terminated
                @test length(BMO.rays(beam)) == 2
                @test BMO.shape(BMO.intersection(last(BMO.rays(beam)))) === BMO.shape(lp.back)
            end
        end
    end

    @testset "Linear polarizer - tilted plate offset" begin
        n = λ -> 1.5
        tf = 1.6mm
        tb = 1.0mm
        D = 25.4mm
        λ = 1e-6

        lp = RoundLinearPolarizer(D, tf, tb, n)
        xrotate3d!(lp, deg2rad(30))
        system = System([lp])

        p1 = [0.0, -10mm, 0.0]
        d1 = [0.0, 1.0, 0.0]
        beam = Beam(p1, d1, λ, [1.0, 0, 0])
        solve_system!(system, beam)

        lastray = last(BMO.rays(beam))
        dir_out = BMO.direction(lastray)
        pos_out = BMO.position(lastray)

        @test dir_out ≈ d1 atol = 1e-12

        # perpendicular distance between the incoming line (p1, d1) and the outgoing line
        delta = pos_out .- p1
        perp = delta .- dot(delta, d1) .* d1
        offset = norm(perp)

        θi = deg2rad(30)
        expected_offset = (tf + tb) * sin(θi) * (1 - cos(θi) / sqrt(1.5^2 - sin(θi)^2))
        @test offset ≈ expected_offset atol = 1e-9
    end

    @testset "Linear polarizer - beamlets" begin
        n = λ -> 1.5
        tf = 1.6mm
        tb = 1.0mm
        D = 25.4mm
        λ = 633e-9
        w0 = 0.5mm
        pos = [0, -10mm, 0]
        dir = [0, 1, 0]

        lp = RoundLinearPolarizer(D, tf, tb, n)
        system = System([lp])

        @testset "GaussianBeamlet" begin
            gauss = GaussianBeamlet(pos, dir, λ, w0)
            @test_nowarn solve_system!(system, gauss)
        end

        @testset "AstigmaticGaussianBeamlet - transmitted" begin
            agb = AstigmaticGaussianBeamlet(pos, dir, λ, w0; E0 = [1, 0, 0], support = [1, 0, 0])
            @test_nowarn solve_system!(system, agb)
            lengths = map(b -> length(BMO.rays(b)), BMO._component_beams(agb))
            @test all(==(4), lengths)
        end

        @testset "AstigmaticGaussianBeamlet - blocked" begin
            agb = AstigmaticGaussianBeamlet(pos, dir, λ, w0; E0 = [0, 0, 1], support = [1, 0, 0])
            @test_nowarn solve_system!(system, agb)
            lengths = map(b -> length(BMO.rays(b)), BMO._component_beams(agb))
            @test all(==(2), lengths)
            @test BMO.shape(BMO.intersection(last(BMO.rays(agb.c)))) === BMO.shape(lp.front)
        end
    end
end

end # MODULE