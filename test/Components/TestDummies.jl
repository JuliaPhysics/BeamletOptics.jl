module TestDummies

using BeamletOptics
using Test

const BMO = BeamletOptics

@testset "Dummy objects" begin
    # Setup dummy cube and test beam
    cube_shape = BMO.CubeMesh(1)
    translate3d!(cube_shape, -[0.5, 0, 0.5])
    translate3d!(cube_shape, [0, 5, 0])
    @testset "IntersectableObject" begin
        beam = Beam([0, 0, 0], [0, 1, 0], 1e-6)
        intersectable = IntersectableObject(cube_shape)
        system = System([intersectable])
        solve_system!(system, beam)
        # Test nothing interaction
        @test length(BMO.rays(beam)) == 1
        @test BMO.object(BMO.intersection(last(BMO.rays(beam)))) == intersectable
        @test BMO.shape(BMO.intersection(last(BMO.rays(beam)))) == cube_shape
        @test isnothing(BMO.interact3d(system, intersectable, beam, first(BMO.rays(beam))))
    end

    @testset "NonInteractableObject" begin
        beam = Beam([0, 0, 0], [0, 1, 0], 1e-6)
        noninteract = NonInteractableObject(cube_shape)
        system = System([noninteract])
        solve_system!(system, beam)
        # Test nothing interaction and intersection
        @test length(BMO.rays(beam)) == 1
        @test isnothing(BMO.intersection(last(BMO.rays(beam))))
        @test isnothing(BMO.interact3d(system, noninteract, beam, first(BMO.rays(beam))))
    end

    @testset "KM100CPMount" begin
        mount = BMO.KM100CPMount()
        @test mount isa NonInteractableObject
        @test BMO.position(mount) == zeros(3)
        # post base 81.8 mm below the mirror center
        zmin, zmax = extrema(BMO.vertices(BMO.shape(mount))[:, 3])
        @test zmin ≈ -0.0818 atol=1e-4
        @test zmax ≈ 0.0264 atol=1e-4
    end
end

end