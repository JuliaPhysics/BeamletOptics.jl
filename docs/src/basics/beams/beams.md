```@setup beams
include(joinpath(@__DIR__, "..", "..", "assets", "cond_save.jl"))

beam_showcase_dir = joinpath(@__DIR__, "..", "..", "assets", "beam_renders")

conditional_include(joinpath(beam_showcase_dir, "beam_showcase.jl"))
```

# Beams

As mentioned in the [Rays](@ref) section, a beam within the context of this package serves as a data structure for storing collections of rays, forming the backbone of the simulation framework. Beams are intended to be designed as [AbstractTrees](https://github.com/JuliaCollections/AbstractTrees.jl) to allow for ray bifurcations, e.g. in the case of optical elements such as beamsplitters. The [`solve_system!`](@ref) function relies on this data structure to perform ray tracing computations within optical systems. 

To ensure compatibility and extensibility, beam types must adhere to the [`BeamletOptics.AbstractBeam`](@ref) interface. Refer to its documentation for more information.

## Basic beam

A minimal implementation of the [`BeamletOptics.AbstractBeam`](@ref) type is provided by the [`Beam`](@ref). It can be used to store a light path through an optical system. If the beam is split, its children will be recursively traced until all paths are solved.

```@docs; canonical=false
Beam
```

A ray tracing example through an arbitrary system using a [`Beam`](@ref) is shown below. Individual [`Ray`](@ref) segments are marked by their starting position and direction. The [Beam expander](@ref) and [Miniature microscope](@ref) tutorial covers the use of the [`Beam`](@ref) in more detail. 

![Beam structure](beam_showcase.png)