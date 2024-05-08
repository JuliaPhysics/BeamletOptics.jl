#=
Types:
    AbstractShape
        AbstractSDF
        AbstractMesh
    AbstractObject
        AbstractReflectiveOptic
        AbstractRefractiveOptic
        AbstractDetector
        AbstractBeamSplitter
        AbstractObjectGroup
    AbstractRay
    AbstractBeam
    AbstractSystem
    Intersection
    Interaction

Core Functions:
    intersect3d(AbstractObject, AbstractRay)
    interact3d(AbstractSystem, AbstractObject, AbstractBeam, AbstractRay)
    rotate3d!(AbstractObject, axis, angle)
    translate3d!(AbstractObject, offset)
    trace_system!(system, Beam)
    trace_system!(system, GaussianBeamlet)
    retrace_system!(system, Beam)
    retrace_system!(system, GaussianBeamlet)
=#

# Order of inclusion matters!
include("AbstractShape.jl")
include("AbstractObject.jl")
include("AbstractRay.jl")
include("AbstractBeam.jl")
include("AbstractSystem.jl")
include("AbstractUtils.jl")