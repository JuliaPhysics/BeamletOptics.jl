# Coatings main entrypoint. Includes subfiles in the correct dependency order.

include("Coatings/Models.jl")
include("Coatings/BoundaryPhysics.jl")
include("Coatings/StandaloneCoating.jl")
include("Coatings/CoatedComponents.jl")
