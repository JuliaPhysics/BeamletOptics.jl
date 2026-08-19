module TestMisc

using BeamletOptics
using Test

const BMO = BeamletOptics

@testset "Render" begin
    axis = nothing
    cube = BMO.CubeMesh(1)
    @test_throws BMO.MissingBackendError render!(axis, cube)
    @test_throws BMO.MissingBackendError BMO.get_view(axis)
    @test_throws BMO.MissingBackendError BMO.set_view(axis, [1 1; 0 0])
    @test_throws BMO.MissingBackendError BMO.hide_axis(axis, true)
end

@testset "Aqua" begin
    if Base.find_package("Aqua") !== nothing
        @eval using Aqua
        @eval Aqua.test_all($BMO; ambiguities = false)
    end
end
   
end