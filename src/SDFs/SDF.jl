#=
Generate a hierarchy of all AbstractSDFs via:

    list_subtypes(BeamletOptics.AbstractSDF)

=#

# Order of inclusion matters!
include("AbstractSDF.jl")
include("UnionSDF.jl")
include("PrimitiveSDF.jl")
include("SphericalLensSDF.jl")
include("MeniscusLensSDF.jl")
include("AsphericalLensSDF.jl")