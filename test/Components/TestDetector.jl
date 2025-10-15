module TestDetector

using BeamletOptics
using GeometryBasics
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

        # get numerical first zero of the Airy-disk
        ci = argmax(I_num)
        ix_ctr, jx_ctr = Tuple(ci)
        num_min = x[argmin(I_num[:, jx_ctr])]   # first zero through the centre column

        # theoretical Airy-disk 1st zero
        airy_min = 1.22*λ*200e-3/D - x_shift

        @test abs(num_min) ≈ airy_min rtol=1e-2
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
    @testset "Pre-Beamsplitter tests with seperate beams" begin
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

@testset "Issue #42" begin
    # Tests regressions of https://github.com/JuliaPhysics/BeamletOptics.jl/issues/42
    
    @testset "Testing RayHit mutability / retracing" begin
        # Setup system
        m1 = RoundPlanoMirror(25mm, 5mm)
        pd = Detector(50mm)
        sys = System([m1, pd])
        translate_to3d!(m1, [0,50mm,0])
        translate_to3d!(pd, [0,50mm,50mm])
        xrotate3d!(m1, deg2rad(-40))
        xrotate3d!(pd, deg2rad(90))
        beam = Beam([0,0,0], [0, 1, 0])
        # get hits at every step and once after all (hits_ref)
        hits_std = Vector{Point2{Float64}}()
        for i = 1:11
            solve_system!(sys, beam)
            # get point directly
            push!(hits_std, last(spot_diagram(pd)))
            # rotate mirror
            if i != 11
                xrotate3d!(m1, deg2rad(-1))
            end
        end
        hits_ref = spot_diagram(pd)
        # test correct number of unique hits
        @test length(hits_std) == length(hits_ref)
        @test length(hits_std) == length(unique(hits_std))
        @test length(hits_ref) == length(unique(hits_ref))
    end
end

end # TESTSET

end # MODULE