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
end

end # MODULE