module BeamletOptics

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
using ForwardDiff: gradient

import Base: length, push!, empty!

# Do not change order of inclusion!
include("Constants.jl")
include("Utils.jl")
include("AbstractTypes/AbstractTypes.jl")
include("Rays.jl")
include("PolarizedRays.jl")
include("Beam.jl")
include("BeamGroups.jl")
include("Mesh.jl")
include("SDFs/SDF.jl")
include("Gaussian.jl")
include("System.jl")
include("OpticalComponents/Components.jl")
include("Groups.jl")
include("Render.jl")
include("Exports.jl")

include("Workloads/precompile.jl")

end # module
