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
    
    # Power at each output port of the MZI: the fields of both arms that leave in its direction
    function port_powers(beam)
        ports = Dict{Vector{Float64}, Vector{ComplexF64}}()
        for leaf in Leaves(beam)
            E, d = BMO.polarization(first(BMO.rays(leaf))), BMO.direction(first(BMO.rays(leaf)))
            k = round.(collect(d), digits = 6)
            ports[k] = get(ports, k, zeros(ComplexF64, 3)) + E
        end
        return sort([norm(E)^2 for E in values(ports)])
    end

    @testset "z-polarized ray along y-axis" begin
        # Solve with z-polarized ray along y-axis
        BMO.polarization!(ray, [0, 0, 1])
        solve_system!(system, beam)

        # Extract E0s: t - transmitted, r - reflected
        t = BMO.polarization(beam.children[1].rays[1])
        r = BMO.polarization(beam.children[2].rays[1])

        # Lossless 50:50 split
        @test norm(t)^2 ≈ 0.5
        @test norm(r)^2 ≈ 0.5
        # The balanced arms interfere: all power leaves one port (r' = -r, Stokes relations)
        @test port_powers(beam) ≈ [0, 1] atol = 1e-12
    end

    @testset "x-polarized ray along y-axis" begin
        # Test num. of leaves before solving again
        @test length(collect(Leaves(beam))) == 4
        # Retrace with x-polarized ray along y-axis
        BMO.polarization!(ray, [1, 0, 0])
        solve_system!(system, beam)

        t = BMO.polarization(beam.children[1].rays[1])
        r = BMO.polarization(beam.children[2].rays[1])
        @test norm(t)^2 ≈ 0.5
        @test norm(r)^2 ≈ 0.5
        @test port_powers(beam) ≈ [0, 1] atol = 1e-12
    end
end

end # MODULE