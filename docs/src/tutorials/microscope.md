```@setup plots
using GLMakie, BeamletOptics

const BMO = BeamletOptics

tutorial_dir = joinpath(@__DIR__, "..", "assets", "ms_assets")

include(joinpath(tutorial_dir, "miniscope_showcase.jl"))
include(joinpath(tutorial_dir, "miniscope_plots.jl"))
```

# Miniature microscope

The UCLA Miniscope is a lightweight microscope that utilizes 2-photon fluorescence imaging to record neural activity in awake, freely moving mice [Madruga:2024](@cite). This beginner tutorial aims to reproduce the optical path of the imaging system from the data provided in the [UCLA 2P Miniscope repository](https://github.com/golshanilab/UCLA_2P_Miniscope). You will learn how to:

1. Define optical components (e.g. a [`SphericalDoubletLens`](@ref))
2. Position the mentioned components using the kinematic API
3. Define an optical [`System`](@ref)
4. Add a source of [`Ray`](@ref)s 
5. Trace a [`Beam`](@ref) through the optical system
6. Visualize the results 

All specifications for the optical system and CAD files are taken from this [repository](https://github.com/golshanilab/UCLA_2P_Miniscope) under GNU GPL-3.0.

![Test](test.png)

## Thingy