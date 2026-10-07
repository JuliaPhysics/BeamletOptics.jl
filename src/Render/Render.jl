#=
Entry points of the render pipeline with their docstrings. Without a loaded `Makie` backend they
throw a `MissingBackendError`; the `Makie` extension (ext/) implements them, with one file per file
below where the names match.
=#

# Order of inclusion matters!
include("RenderCore.jl")
include("RenderWavelength.jl")
include("RenderCamera.jl")
include("RenderBoundingSphere.jl")
include("RenderHandles.jl")
include("RenderLive.jl")
include("RenderLook.jl")
