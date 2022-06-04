module SCDI

using LinearAlgebra, FileIO

# include("ComplexRayTracing.jl")

include("Ray Tracing\\Geometry.jl")
include("Ray Tracing\\Rays.jl")
include("Ray Tracing\\Tracing.jl")
include("Ray Tracing\\Utils.jl")
include("Ray Tracing\\Optics.jl")

end # module
