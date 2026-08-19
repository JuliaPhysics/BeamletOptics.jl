module TestMachZehnder

using BeamletOptics
using Test
using LinearAlgebra
using AbstractTrees

const BMO = BeamletOptics

@testset "Mach-Zehnder interferometer" begin
    # setup MZI
    m1 = SquarePlanoMirror2D(BMO.inch)
    m2 = SquarePlanoMirror2D(BMO.inch)
    b1 = ThinBeamsplitter(BMO.inch, reflectance = 0.5)
    b2 = ThinBeamsplitter(BMO.inch, reflectance = 0.5)
    
    system = StaticSystem([m1, m2, b1, b2])
    
    translate3d!(b1, [0 * BMO.inch, 0 * BMO.inch, 0])
    translate3d!(b2, [2 * BMO.inch, 2 * BMO.inch, 0])
    translate3d!(m1, [0 * BMO.inch, 2 * BMO.inch, 0])
    translate3d!(m2, [2 * BMO.inch, 0 * BMO.inch, 0])
    
    # Rotate with consideration to mirror/bs normal
    zrotate3d!(b1, deg2rad(360 - 135))
    zrotate3d!(b2, deg2rad(45))
    zrotate3d!(m1, deg2rad(360 - 135))
    zrotate3d!(m2, deg2rad(45))
    
    ray = PolarizedRay([0, -0.1, 0], [0.0, 1.0, 0], 1000e-9, [0, 0, 1])
    beam = Beam(ray)
    
    @testset "z-polarized ray along y-axis" begin
        # Solve with z-polarized ray along y-axis
        BMO.polarization!(ray, [0, 0, 1])
        solve_system!(system, beam)
        
        # Extract E0s: t - transmitted, r - reflected
        t = BMO.polarization(beam.children[1].rays[1])
        r = BMO.polarization(beam.children[2].rays[1])
        tr = BMO.polarization(beam.children[1].rays[2])
        rr = BMO.polarization(beam.children[2].rays[2])
        trt = BMO.polarization(beam.children[1].children[1].rays[1])
        trr = BMO.polarization(beam.children[1].children[2].rays[1])
        rrt = BMO.polarization(beam.children[2].children[1].rays[1])
        rrr = BMO.polarization(beam.children[2].children[2].rays[1])
        
        # Test phase flips
        @test t[3] ≈ sqrt(2) / 2
        @test r[3] ≈ -sqrt(2) / 2
        @test tr[3] ≈ -sqrt(2) / 2
        @test rr[3] ≈ sqrt(2) / 2
        @test trt ≈ -rrr
        @test trr ≈ rrt
    end
    
    @testset "x-polarized ray along y-axis" begin
        # Test num. of leaves before retracing
        @test length(collect(Leaves(beam))) == 4
        # Retrace with z-polarized ray along y-axis
        BMO.polarization!(ray, [1, 0, 0])
        solve_system!(system, beam)
        
        # Extract E0s: t - transmitted, r - reflected
        t = BMO.polarization(beam.children[1].rays[1])
        r = BMO.polarization(beam.children[2].rays[1])
        tr = BMO.polarization(beam.children[1].rays[2])
        rr = BMO.polarization(beam.children[2].rays[2])
        trt = BMO.polarization(beam.children[1].children[1].rays[1])
        trr = BMO.polarization(beam.children[1].children[2].rays[1])
        rrt = BMO.polarization(beam.children[2].children[1].rays[1])
        rrr = BMO.polarization(beam.children[2].children[2].rays[1])
        
        # Test phase flips
        @test t[1] ≈ sqrt(2) / 2
        @test r[2] ≈ -sqrt(2) / 2
        @test tr ≈ r
        @test rr ≈ t
        @test trt ≈ -rrr
        @test trr ≈ rrt
    end
end

end # MODULE