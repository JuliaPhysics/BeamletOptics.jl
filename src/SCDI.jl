module SCDI

using LinearAlgebra: norm, normalize, normalize!, dot, cross, I
using MarchingCubes: MC, march
using Trapz: trapz
using PrecompileTools: @setup_workload, @compile_workload
using StaticArrays: @SArray, @SArray, SMatrix, SArray
using GeometryBasics: Point3, Point2, Mat
using AbstractTrees: AbstractTrees, parent, children, NodeType, nodetype, nodevalue, print_tree, HasNodeType, Leaves, StatelessBFS, PostOrderDFS, PreOrderDFS

import Base: length

# Do not change order of inclusion!
include("Ray Tracing/Constants.jl")
include("Ray Tracing/Utils.jl")
include("Ray Tracing/AbstractTypes/AbstractTypes.jl")
include("Ray Tracing/Rays.jl")
include("Ray Tracing/PolarizedRays.jl")
include("Ray Tracing/Beam.jl")
include("Ray Tracing/Mesh.jl")
include("Ray Tracing/SDFs/SDF.jl")
include("Ray Tracing/Gaussian.jl")
include("Ray Tracing/System.jl")
include("Ray Tracing/Interactions.jl")
include("Ray Tracing/Groups.jl")
include("Ray Tracing/Render.jl")
include("Ray Tracing/Components.jl")
include("Ray Tracing/Doublets.jl")

if get(ENV, "CI", "false") == "false"
    @setup_workload begin
        # setup dummy workload
        _vertices = [1 1 0
            1 -1 0
            -1 -1 0
            -1 1 0]
        _faces = [1 2 3
            3 4 1]
        _pos = [0, 0, 0]
        _dir = Matrix{Int}(I, 3, 3)
        _scale = 1
        @compile_workload begin
            # execute workload
            plane = SCDI.Mirror(
                SCDI.Mesh{Float64}(
                    _vertices,
                    _faces,
                    _dir,
                    _pos,
                    _scale))
            SCDI.translate3d!(plane, [0, 0, 1])
            SCDI.xrotate3d!(plane, π / 2)
            SCDI.yrotate3d!(plane, π / 4)
            SCDI.zrotate3d!(plane, π / 2)
        end
    end
end

end # module
