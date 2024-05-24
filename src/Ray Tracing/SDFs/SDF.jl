#=
Types:
    AbstractSDF
        UnionSDF
        SphereSDF
        CylinderSDF
        CutSphereSDF
        RingSDF
        AbstractRotationallySymmetricLensSDF
            AbstractSphericalLensSDF

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
include("AsphericalLensSDF.jl")