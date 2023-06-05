module SCDI

using LinearAlgebra, FileIO, UUIDs, MarchingCubes
using PrecompileTools
using MakieCore

import Base: length

# Do not change order of inclusion!
include("Ray Tracing/Utils.jl")
include("Ray Tracing/Types.jl")
include("Ray Tracing/Rays.jl")
include("Ray Tracing/Gaussian.jl")
include("Ray Tracing/Mesh.jl")
include("Ray Tracing/SDF.jl")
include("Ray Tracing/Sphere.jl")
include("Ray Tracing/Surfaces.jl")
include("Ray Tracing/Intersections.jl")
include("Ray Tracing/System.jl")
include("Ray Tracing/Interactions.jl")
include("Ray Tracing/Render.jl")

if get(ENV, "CI", "false") == "false"
    @setup_workload begin
        # setup dummy workload
        _vertices = [
            1 1 0
            1 -1 0
            -1 -1 0
            -1 1 0
        ]
        _faces = [
            1 2 3
            3 4 1
        ]
        _pos = [0, 0, 0]
        _dir = Matrix{Int}(I, 3, 3)
        _scale = 1
        @compile_workload begin
            # execute workload
            plane = SCDI.Mirror(
                uuid4(),
                SCDI.Mesh{Float64}(
                    uuid4(),
                    _vertices,
                    _faces,
                    _dir,
                    _pos,
                    _scale
            ))
            SCDI.translate3d!(shape(plane), [0, 0, 1])
            SCDI.xrotate3d!(shape(plane), π / 2)
            SCDI.yrotate3d!(shape(plane), π / 4)
            SCDI.zrotate3d!(shape(plane), π / 2)
            SCDI.scale3d!(shape(plane), 2)
        end
    end
end

end # module
