#=
Types:
    AbstractShape
        AbstractSDF
        AbstractMesh
    AbstractObject
        AbstractReflectiveOptic
        AbstractRefractiveOptic
        AbstractDetector
        AbstractBeamsplitter
        AbstractObjectGroup
    AbstractRay
    AbstractBeam
    AbstractSystem
    Intersection
    Interaction
    Hint

Core Functions:
    intersect3d(AbstractObject, AbstractRay)
    interact3d(AbstractSystem, AbstractObject, AbstractBeam, AbstractRay)
    rotate3d!(AbstractObject, axis, angle)
    translate3d!(AbstractObject, offset)
    trace_system!(system, Beam)
    trace_system!(system, GaussianBeamlet)
=#

# A call of `intersect3d` on an abstract shape or object is dispatched at runtime. Inference must not
# compile the few generic methods it matches for abstract arguments: their code would call `dot`,
# `normalize` etc. on values of unknown type, and loading a package that adds such methods, e.g.
# Makie, would invalidate the precompiled tracing code. Calls with concrete arguments match one method.
Base.Experimental.@max_methods 1 function intersect3d end

# Order of inclusion matters!
include("AbstractKinematicTrait.jl")
include("AbstractShape.jl")
include("AbstractObject.jl")
include("AbstractShapeTrait.jl")
include("AbstractRay.jl")
include("AbstractBeam.jl")
include("AbstractBeamGroup.jl")
include("AbstractGaussian.jl")
include("AbstractSystem.jl")
include("AbstractBoundingSphere.jl")
include("AbstractUtils.jl")
