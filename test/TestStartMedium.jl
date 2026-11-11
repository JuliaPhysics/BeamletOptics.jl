module TestStartMedium

using BeamletOptics
using Test
using LinearAlgebra

const BMO = BeamletOptics
const mm = 1e-3
const λ0 = 1e-6

@testset "isinside" begin
    # A signed distance function: a lens
    lens = SphericalLens(50mm, -50mm, 5mm, 25mm, ConstantRefractiveIndex(1.5))
    @test BMO.isinside(BMO.shape(lens), [0, 2.5mm, 0])
    @test !BMO.isinside(BMO.shape(lens), [0, -1mm, 0])
    @test !BMO.isinside(BMO.shape(lens), [20mm, 2.5mm, 0])
    # A closed mesh, also on its axes and diagonals, where a parity ray could hit an edge
    cuboid = BMO.CuboidMesh(10mm, 10mm, 10mm)
    for p in ([5mm, 5mm, 5mm], [1mm, 1mm, 1mm], [9mm, 5mm, 1mm], [5mm, 5mm, 9.99mm])
        @test BMO.isinside(cuboid, p)
    end
    for p in ([-1mm, 5mm, 5mm], [5mm, 11mm, 5mm], [15mm, 15mm, 15mm])
        @test !BMO.isinside(cuboid, p)
    end
    # A flat mesh has no volume, a shape without a method contains no point
    @test !BMO.isinside(BMO.RectangularFlatMesh(10mm, 10mm), [0, 0, 0])
    @test BMO.isinside(BMO.SingleBoundingSphere([0, 0, 0], 1.0), [0.5, 0, 0])
    @test BMO.isinside(BMO.NoBoundingSphere(), [1e3, 0, 0])
end

# A block of the index n, its faces at y = 0 and y = thickness
block(n; size = 40mm, thickness = 10mm) =
    RectangularCompensatorPlate(size, size, thickness, ConstantRefractiveIndex(n))

@testset "A ray that starts in glass" begin
    n = 1.5
    glass = block(n)
    system = System([glass])
    # Leaving through the back face it refracts from the glass into air
    θ = deg2rad(20)
    beam = Beam([0, 5mm, 0], [sin(θ), cos(θ), 0], λ0)
    solve_system!(system, beam)
    @test BMO.refractive_index(first(rays(beam))) == n
    @test asin(BMO.direction(last(rays(beam)))[1]) ≈ asin(n * sin(θ))
    # The same holds for a block that is turned and moved
    zrotate3d!(glass, deg2rad(30))
    translate3d!(glass, [3mm, -2mm, 1mm])
    d = orientation(glass) * [sin(θ), cos(θ), 0]
    beam = Beam(position(glass) + orientation(glass) * [0, 5mm, 0], d, λ0)
    solve_system!(system, beam)
    @test dot(BMO.direction(last(rays(beam))), orientation(glass)[:, 1]) ≈ n * sin(θ)
    # Moved out of the glass, the source starts in the medium of the system again
    beam = Beam([0, -50mm, 0], [0, 1, 0], λ0)
    BMO.refractive_index!(first(rays(beam)), n)
    solve_system!(system, beam)
    @test BMO.refractive_index(first(rays(beam))) == 1
end

@testset "A ray that starts on a surface" begin
    # On the front vertex of a lens: into the lens it starts in the glass, away from it in air
    lens = SphericalLens(50mm, -50mm, 5mm, 25mm, ConstantRefractiveIndex(1.5))
    for (d, n) in (([0, 1, 0], 1.5), ([0, -1, 0], 1.0))
        beam = Beam([0, 0, 0], d, λ0)
        solve_system!(System([lens]), beam)
        @test BMO.refractive_index(first(rays(beam))) == n
    end
    # On the face of a mesh plate
    plate = block(1.5)
    for (d, n) in (([0, 1, 0], 1.5), ([0, -1, 0], 1.0))
        beam = Beam([0, 0, 0], d, λ0)
        solve_system!(System([plate]), beam)
        @test BMO.refractive_index(first(rays(beam))) == n
    end
end

@testset "A ray that starts in a lens within water" begin
    n_glass, n_water = 1.5, 1.33
    # Water from y = -40 mm to 0, the lens from y = -20 mm to -16 mm within it
    water = block(n_water; size = 80mm, thickness = 40mm)
    translate3d!(water, [0, -40mm, 0])
    lens = SphericalLens(Inf, Inf, 4mm, 20mm, ConstantRefractiveIndex(n_glass))
    translate3d!(lens, [0, -20mm, 0])
    θ = deg2rad(15)
    beam = Beam([0, -18mm, 0], [sin(θ), cos(θ), 0], λ0)
    solve_system!(System([water, lens]), beam)
    @test BMO.refractive_index(first(rays(beam))) == n_glass
    # Glass into water, then water into air
    @test BMO.refractive_index(rays(beam)[2]) == n_water
    @test BMO.direction(rays(beam)[2])[1] ≈ n_glass * sin(θ) / n_water
    @test BMO.direction(last(rays(beam)))[1] ≈ n_glass * sin(θ)
end

@testset "Beamlets and sources start in glass" begin
    n = 1.5
    system = System([block(n)])
    gauss = GaussianBeamlet([0, 5mm, 0], [0, 1, 0], λ0, 0.5mm)
    solve_system!(system, gauss)
    @test all(b -> BMO.refractive_index(first(rays(b))) == n, BMO._component_beams(gauss))
    source = CollimatedSource([0, 5mm, 0], [0, 1, 0], 10mm, λ0; num_rings = 2)
    solve_system!(system, source)
    @test all(b -> BMO.refractive_index(first(rays(b))) == n, BMO.beams(source))
end

end # module
