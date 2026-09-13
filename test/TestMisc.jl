module TestMisc

using BeamletOptics
using Test
using LinearAlgebra: I

const BMO = BeamletOptics

@testset "Render" begin
    axis = nothing
    cube = BMO.CubeMesh(1)
    @test_throws BMO.MissingBackendError render!(axis, cube)
    @test_throws BMO.MissingBackendError BMO.get_view(axis)
    @test_throws BMO.MissingBackendError BMO.set_view(axis, [1 1; 0 0])
    @test_throws BMO.MissingBackendError BMO.set_view(axis, [0,0,0], [0,1,0], [0,0,1])
    @test_throws BMO.MissingBackendError BMO.set_orthographic(axis)
    @test_throws BMO.MissingBackendError BMO.hide_axis(axis, true)
    @test_throws BMO.MissingBackendError BMO.arrow!(axis, [0,0,0], [1,0,0])
    @test_throws BMO.MissingBackendError BMO.render_lcs!(axis)
    @test_throws BMO.MissingBackendError BMO.render_lcs!(axis, [0,0,0], Matrix{Float64}(I, 3, 3))
    @test_throws BMO.MissingBackendError BMO.render_lcs!(axis, RoundPlanoMirror(25e-3, 5e-3))
    @test_throws BMO.MissingBackendError BMO.look_at!(axis, [0,0,0], [1,0,0])
end

@testset "Aqua" begin
    using Aqua
    Aqua.test_all(BMO)
end
   
end