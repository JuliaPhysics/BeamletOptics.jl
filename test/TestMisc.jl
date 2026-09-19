module TestMisc

using BeamletOptics
using Aqua
using Test

@testset "Aqua" begin
    Aqua.test_all(BeamletOptics)
end
   
end