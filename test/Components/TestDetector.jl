module TestDetector

using BeamletOptics
using Test

const BMO = BeamletOptics

const mm = 1e-3

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

@testset "Point-spread functions" begin
    # parameters for an almost thin-lens
    l = 1e-3
    R1 = 100e-3
    R2 = Inf
    d = 25.4e-3
    n = 1.5
    λ = 1e-6    
    D = 15e-3
    num_rays = 1000
    
    # plane wave source
    cs = UniformDiscSource([0, -10e-3, 0], [0, 1, 0], D, λ; num_rays)
    
    # test lens
    lens = SphericalLens(R1, R2, l, d, x -> n)
    
    # PSF detector
    psfd = PSFDetector(10e-3)
    translate3d!(psfd, [0, 200e-3 + 0.13e-3, 0])

    @testset "Airy-disc test" begin
        # build system and solve it
        sys = System([lens, psfd])
        solve_system!(sys, cs)

        # get PSF
        x, y, I_num = intensity(psfd; n=500, crop_factor=5, center=:bbox)

        # get numerical first zero of the Airy-disk
        ci = argmax(I_num)
        ix_ctr, jx_ctr = Tuple(ci)
        num_min = x[argmin(I_num[:, jx_ctr])]   # first zero through the centre column

        # theoretical Airy-disk 1st zero
        airy_min = 1.22*λ*200e-3/D

        @test abs(num_min) ≈ airy_min rtol=1e-2
    end

    @testset "PSFDetector reset" begin
        # Test if data in psfd
        @test length(psfd.data) == num_rays
        # Empty psfd and test
        empty!(psfd)
        @test isempty(psfd.data)
    end
end


end