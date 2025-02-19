module SCDI

using LinearAlgebra: norm, normalize, normalize!, dot, cross, I, tr
using MarchingCubes: MC, march
using Trapz: trapz
using PrecompileTools: @setup_workload, @compile_workload
using StaticArrays: @SArray, @SArray, SMatrix, SArray
using GeometryBasics: Point3, Point2, Mat
using AbstractTrees: AbstractTrees, parent, children, NodeType, nodetype, nodevalue, print_tree, HasNodeType, Leaves, StatelessBFS, PostOrderDFS, PreOrderDFS, TreeIterator
using InteractiveUtils: subtypes
using FileIO: load
using MeshIO

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
include("Ray Tracing/OpticalComponents/Components.jl")
include("Ray Tracing/Groups.jl")
include("Ray Tracing/Render.jl")

if get(ENV, "CI", "false") == "false"
    include("CompileWorkload.jl")
end

end # module
