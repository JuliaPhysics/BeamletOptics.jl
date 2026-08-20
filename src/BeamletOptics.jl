module BeamletOptics

using LinearAlgebra: norm, normalize, normalize!, dot, cross, I, tr, eigen, Symmetric
using MarchingCubes: MC, march
using Trapz: trapz
using PrecompileTools: @setup_workload, @compile_workload
using StaticArrays: @SArray, @SVector, SMatrix, SArray, SVector
using GeometryBasics: Point3, Point2, Mat
using AbstractTrees: AbstractTrees, parent, children, NodeType, nodetype, nodevalue,
                     print_tree, HasNodeType, Leaves, StatelessBFS, PostOrderDFS,
                     PreOrderDFS, TreeIterator
using InteractiveUtils: subtypes
using FileIO: load
using MeshIO
using ForwardDiff: gradient
using Random

import Base: length, push!, empty!, position

# Do not change order of inclusion!
include("Constants.jl")
include("Config.jl")
using .Config: get_invariant_threshold, set_invariant_threshold!,
             get_sdf_surface_threshold, get_sdf_raymarch_eps, get_sdf_inside_step,
             get_internal_reflection_threshold, get_line_plane_intersection_threshold,
             get_orthogonality_threshold, get_default_r_max, get_default_depth_max,
             get_default_wavelength, get_default_waist, get_default_power,
             get_coincident_boundary_tolerance
include("Utils/Utils.jl")
include("AbstractTypes/AbstractTypes.jl")
include("Rays.jl")
include("PolarizedRays.jl")
include("Transitions.jl")
include("SurfaceInteractions.jl")
include("Beam.jl")
include("Gaussian.jl")
include("AstigmaticGaussian.jl")
include("BeamGroups.jl")
include("Mesh.jl")
include("SDFs/SDF.jl")
include("System.jl")
include("OpticalComponents/Components.jl")
include("ObjectGroups.jl")
include("Render.jl")
include("Exports.jl")

include("Workloads/precompile.jl")

end # module
