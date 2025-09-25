module TestMichelson

using BeamletOptics
using Test

const BMO = BeamletOptics

const mm = 1e-3

@testset "Michelson interferometer" begin
    @testset "Testing full Michelson interferometer model" begin
        # setup Michelson Interferometer
        l_0 = 0.1
        pd_size = BMO.inch / 5
        pd_resolution = 100
        m1 = SquarePlanoMirror2D(BMO.inch)
        m2 = SquarePlanoMirror2D(BMO.inch)
        bs = ThinBeamsplitter(BMO.inch, reflectance = 0.5)
        pd = Detector(pd_size)
        translate3d!(m1, [l_0, 0, 0])
        translate3d!(m2, [0, l_0, 0])
        translate3d!(pd, [-l_0, 0, 0])
        zrotate3d!(bs, deg2rad(45))
        zrotate3d!(m1, deg2rad(90))
        zrotate3d!(pd, deg2rad(90))
        
        system = System([m1, m2, bs, pd])
        
        # Test correct values for reflectivity/transmission
        @test isvalid(bs)
        
        @testset "Equal armlength MI - integrated power" begin
            # setup 635 nm laser with 0.1 mm waist for fast divergence
            λ = 635e-9
            P_0 = 5e-3
            beam = GaussianBeamlet([0, -l_0, 0], [0, 1.0, 0], λ, 1e-4, P0 = P_0)
            
            # Shift mirror #2 by -λ to +λ
            lambdas = LinRange(-λ, λ, 200)
            
            path_length_numerical = zeros(length(lambdas))
            optical_pwr_numerical = zeros(length(lambdas))
            
            for (i, lambda) in enumerate(lambdas)
                translate_to3d!(m2, [0, l_0, 0] + [0, lambda, 0])
                empty!(pd)
                solve_system!(system, beam)
                
                # Moving mirror path length
                path_length_numerical[i] = length(beam.children[1].children[2])
                optical_pwr_numerical[i] = BMO.optical_power(
                    pd;
                    # restore pre-v0.11 behavior, e.g. no autolims
                    n=pd_resolution,
                    x_min=-pd_size/2,
                    x_max=pd_size/2,
                    z_min=-pd_size/2,
                    z_max=pd_size/2
                )
            end
            path_length_analytical = @. 2 * lambdas + 4l_0
            optical_pwr_analytical = @. P_0 * (1 / 2 * cos(2π * (2lambdas / λ) + π) + 1 / 2)
            
            # Compare correct PD signal and λ shift in moving arm
            @test all(isapprox.(optical_pwr_analytical, optical_pwr_numerical, atol = 5e-6))
            @test all(isapprox.(path_length_analytical, path_length_numerical))
        end
        
        @testset "Unequal armlength MI - electrical field" begin
            λ = 635e-9
            w0 = 1e-4
            P0 = 1e-3
            M2 = 1
            I0 = 2 * P0 / (π * w0^2)
            E0 = BMO.electric_field(I0) * 1 / sqrt(2)^2
            zR = BMO.rayleigh_range(λ, w0, M2)
            
            beam = GaussianBeamlet([0, -l_0, 0], [0, 1.0, 0], λ, w0, P0 = P0, M2 = M2)
            
            # arm length diff
            Δl = 1 * l_0
            translate_to3d!(m2, [0, l_0 + Δl, 0])
            
            # numerical solution
            empty!(pd)
            solve_system!(system, beam)
            
            xs, ys, Efield = electric_field(pd; n=pd_resolution)
            
            # analytical solution
            short_arm = 4l_0
            long_arm = short_arm + 2Δl
            screen = zeros(ComplexF64, length(xs), length(ys))
            for (j, y) in enumerate(ys)
                for (i, x) in enumerate(xs)
                    r = sqrt(x^2 + y^2)
                    screen[i, j] += BMO.electric_field(
                        r, short_arm, E0, w0, λ, M2)
                    screen[i, j] += BMO.electric_field(
                        r, long_arm, E0, w0, λ, M2) * exp(im * pi)
                end
            end
            
            Re_analytical = real.(screen)
            Im_analytical = imag.(screen)
            
            Re_numerical = real(Efield)
            Im_numerical = imag(Efield)
            
            # Compare solutions, units V/m
            @test all(isapprox.(Re_analytical, Re_numerical, atol = 5e-2))
            @test all(isapprox.(Im_analytical, Im_numerical, atol = 5e-2))
        end
    end

    @testset "Testing power conservation after beamsplitting" begin
        # variables
        P0 = 0.5 # W
        l0 = 0.1 # m
        w0 = 0.5e-3
        λ = 1064e-9
        
        bs = ThinBeamsplitter(10e-3)
        pd_resolution = 100
        pd_1 = Detector(10e-3)
        pd_2 = Detector(10e-3)
        
        zrotate3d!(bs, deg2rad(45))
        translate3d!(pd_1, [0, l0, 0])
        zrotate3d!(pd_1, deg2rad(180))
        
        translate3d!(pd_2, [l0, 0, 0])
        zrotate3d!(pd_2, deg2rad(90))
        
        # add BS and PD orientation error
        zrotate3d!(bs, deg2rad(0.017))
        zrotate3d!(pd_1, deg2rad(10))
        xrotate3d!(pd_1, deg2rad(15))
        
        # define system and beams -> solve
        system = System([bs, pd_1, pd_2])
        
        phis = LinRange(0, 2pi, 25)
        p1 = similar(phis)
        p2 = similar(phis)
        
        l1 = GaussianBeamlet([0, -l0, 0], [0, 1.0, 0], λ, w0; P0)
        l2 = GaussianBeamlet([-l0, 0, 0], [1.0, 0, 0], λ, w0; P0)
        
        E0_buffer = l1.E0
        
        for (i, phi) in enumerate(phis)
            # Iterate over relative phase shifts, use retracing
            l1.E0 = E0_buffer * exp(im * phi)
            empty!(pd_1)
            empty!(pd_2)
            solve_system!(system, l1)
            solve_system!(system, l2)
            p1[i] = BMO.optical_power(pd_1)
            p2[i] = BMO.optical_power(pd_2)
            # Test power conservation
            @test p1[i] + p2[i] - 2P0 < 1e-4 # W
        end
    end
end
    
end