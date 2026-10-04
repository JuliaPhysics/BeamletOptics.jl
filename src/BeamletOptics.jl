"""
    BeamletOptics

Non-sequential 3D ray and Gaussian beamlet tracing for optical setups ("BMO").
Documentation: https://juliaphysics.github.io/BeamletOptics.jl/stable/

# AI coding assistants

The package ships an agent skill that teaches assistants such as Claude Code how to use BMO.
Install the copy matching this package version into a project with
[`BeamletOptics.install_agent_skill`](@ref).
"""
module BeamletOptics

using LinearAlgebra: norm, normalize, normalize!, dot, cross, I, eigen, Symmetric, Hermitian, svd
using MarchingCubes: MC, march
using Trapz: trapz
using PrecompileTools: @setup_workload, @compile_workload
using StaticArrays: @SArray, @SVector, SMatrix, SArray, SVector
using GeometryBasics: Point3, Point2, Mat
using AbstractTrees: AbstractTrees, parent, children, isroot, NodeType, nodetype, nodevalue,
                     print_tree, HasNodeType, Leaves, StatelessBFS, PostOrderDFS,
                     PreOrderDFS, TreeIterator
using InteractiveUtils: subtypes
using FileIO: load
using MeshIO
using ForwardDiff: gradient
using Random
using ProgressMeter: Progress, update!, finish!, cancel

import Base: length, push!, empty!, position

# Do not change order of inclusion!
include("Constants.jl")
include("Config.jl")
using .Config: get_invariant_threshold, set_invariant_threshold!,
             get_sdf_surface_threshold, get_sdf_raymarch_eps, get_sdf_inside_step,
             get_internal_reflection_threshold, get_line_plane_intersection_threshold,
             get_orthogonality_threshold, get_default_r_max, get_default_depth_max,
             get_default_wavelength, get_default_waist, get_default_power,
             get_progress_threshold, set_progress_threshold!
include("Utils/Utils.jl")
include("AbstractTypes/AbstractTypes.jl")
include("Rays.jl")
include("PolarizedRays.jl")
include("Beam.jl")
include("Gaussian.jl")
include("AstigmaticGaussian.jl")
include("BeamGroups/BeamGroups.jl")
include("Mesh.jl")
include("SDFs/SDF.jl")
include("System.jl")
include("OpticalComponents/Components.jl")
include("ObjectGroups.jl")
include("Properties.jl")
include("Render.jl")
include("AgentSkill.jl")
include("Exports.jl")

include("Workloads/precompile.jl")

end # module
