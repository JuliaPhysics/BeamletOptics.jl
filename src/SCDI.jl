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
    cm = 1e-2
    splitter_origin = [18.81cm, 23.5cm,0]

    @setup_workload begin
        # setup dummy workload
        cm = 1e-2
        splitter_origin = [18.81cm, 23.5cm,0]
        @compile_workload begin
            # execute workload
            NBK7 = SCDI.DiscreteRefractiveIndex([632.8e-9], [1.51509])

            # Mirror
            rpm = SCDI.RightAnglePrismMirror(25e-3, 25e-3)
            SCDI.zrotate3d!(rpm, deg2rad(45))
            SCDI.translate3d!(rpm, [0,33.5cm,0])

            mirror_assembly = SCDI.ObjectGroup([rpm])

            # Beamsplitter
            cbs = SCDI.CubeBeamsplitter(SCDI.inch, NBK7)
            SCDI.zrotate3d!(cbs, deg2rad(-90))

            splitter_assembly = SCDI.ObjectGroup([cbs])

            # Arms
            m1 = SCDI.RoundPlanoMirror(SCDI.inch, 5e-3)
            SCDI.zrotate3d!(m1, deg2rad(-90))
            SCDI.translate3d!(m1, [22cm,0,0])
            m2 = SCDI.RoundPlanoMirror(SCDI.inch, 5e-3)
            SCDI.zrotate3d!(m2, deg2rad(-90))
            SCDI.translate3d!(m2, [12cm,0,0])

            arm_1 = SCDI.ObjectGroup([m1])
            arm_2 = SCDI.ObjectGroup([m2])

            # PD
            pd = SCDI.Photodetector(8e-3, 200)
            SCDI.translate3d!(pd, [0, -12cm, 0])

            pd_assembly = SCDI.ObjectGroup([pd])

            system = SCDI.System(
                [
                    mirror_assembly,
                    splitter_assembly,
                    arm_1,
                    arm_2,
                    pd_assembly
                ]
            )

            ##
            SCDI.translate_to3d!(mirror_assembly, [0,-10cm,0])
            SCDI.translate_to3d!(splitter_assembly, [18.81cm,23.5cm,0])

            SCDI.translate_to3d!(arm_1, splitter_origin)
            SCDI.translate_to3d!(arm_2, splitter_origin)
            SCDI.translate3d!(arm_1, [3.81cm/2, 0, 0])
            SCDI.translate3d!(arm_2, [0, 3.81cm/2, 0])
            SCDI.zrotate3d!(arm_2, deg2rad(90))

            SCDI.translate_to3d!(pd_assembly, splitter_origin)
            SCDI.translate3d!(pd_assembly, [0, -3.81cm/2, 0])

            ##
            beam = SCDI.GaussianBeamlet([0.,0,0], [0.,1,0], 632.8e-9, 5e-4, M2=2)
            SCDI.solve_system!(system, beam)

            ##
            SCDI.reset_detector!(pd)
            SCDI.translate3d!(m1, [5e-9, 0, 0])
            SCDI.solve_system!(system, beam)
            p = SCDI.optical_power(pd)
        end
    end
end

end # module
