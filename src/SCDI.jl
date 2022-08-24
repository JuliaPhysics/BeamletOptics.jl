module SCDI

using LinearAlgebra, FileIO
using MakieCore
using SnoopPrecompile

# Do not change order of inclusion!
include("Ray Tracing/Utils.jl")
include("Ray Tracing/Rays.jl")
include("Ray Tracing/Mesh.jl")
include("Ray Tracing/Sphere.jl")
include("Ray Tracing/Intersections.jl")
include("Ray Tracing/System.jl")
include("Ray Tracing/Interactions.jl")
include("Ray Tracing/Render.jl")

if get(ENV, "CI", "false") == "false"
    @precompile_setup begin
        # setup dummy workload
        vertices = [
            1 1 0
            1 -1 0
            -1 -1 0
            -1 1 0
        ]
        faces = [
            1 2 3
            3 4 1
        ]
        pos = [0, 0, 0]
        dir = Matrix{Int}(I, 3, 3)
        scale = 1
        @precompile_all_calls begin
            # execute workload
            plane = Mirror{Float64}(SCDI.Mesh{Float64}(
                vertices,
                faces,
                dir,
                pos,
                scale
            ))

            SCDI.translate3d!(plane, [0, 0, 1])
            SCDI.xrotate3d!(plane, π / 2)
            SCDI.yrotate3d!(plane, π / 4)
            SCDI.zrotate3d!(plane, π / 2)
            SCDI.scale3d!(plane, 2)
        end
    end
end

end # module
