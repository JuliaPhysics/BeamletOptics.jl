#=
Types:
    AbstractSDF
        UnionSDF
        SphereSDF
        CylinderSDF
        CutSphereSDF
        RingSDF
        AbstractRotationallySymmetricLensSDF
            AbstractSphericalSurfaceSDF
                ConvexSphericalSurfaceSDF
                ConcaveSphericalSurfaceSDF
            AbstractAsphericalSurfaceSDF
                ConvexAsphericalSurfaceSDF
                ConcaveAsphericalSurfaceSDF
=#

# Order of inclusion matters!
include("AbstractSDF.jl")
include("UnionSDF.jl")
include("PrimitiveSDF.jl")
include("SphericalLensSDF.jl")
include("MeniscusLensSDF.jl")
include("AsphericalLensSDF.jl")