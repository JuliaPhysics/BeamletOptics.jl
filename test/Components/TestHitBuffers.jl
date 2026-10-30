module TestHitBuffers

using BeamletOptics
using Test

const BMO = BeamletOptics

@testset "Hit buffers" begin
    # Hit of a traced ray
    pd = Detector(10e-3)
    translate3d!(pd, [0, 0.1, 0])
    solve_system!(System([pd]), Beam([0.0, 0, 0], [0.0, 1, 0]))
    hit = only(BMO.hits(pd))

    @testset "push! returns the detector" begin
        empty!(pd)
        @test push!(pd, hit) === pd
        @test BMO.hits(pd) == [hit]
        # within `with_hit_buffers` the hit reaches the detector when the buffers are flushed
        empty!(pd)
        BMO.with_hit_buffers() do
            @test push!(pd, hit) === pd
        end
        @test BMO.hits(pd) == [hit]
    end
end

end
