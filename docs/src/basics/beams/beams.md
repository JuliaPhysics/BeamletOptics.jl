```@setup beams
beam_showcase_dir = joinpath(@__DIR__, "..", "..", "assets", "beam_renders")

Main.DocUtils.conditional_include(joinpath(beam_showcase_dir, "beam_showcase.jl"))
```

# Basic beam

A minimal implementation of the [`BeamletOptics.AbstractBeam`](@ref) type is provided by the [`Beam`](@ref). It can be used to store a light path through an optical system. If the beam is split, its children will be recursively traced until all paths are solved.

```@docs; canonical=false
Beam
```

A ray tracing example through an arbitrary system using a [`Beam`](@ref) is shown below. Individual [`Ray`](@ref) segments are marked by their starting position and direction. The [Laser alignment](@ref) and [Miniature microscope](@ref) tutorial covers the use of the [`Beam`](@ref) in more detail. 

![Beam structure](beam_showcase.png)

## Inspecting a traced beam

The rays of a [`Beam`](@ref) are available via `rays`, its child beams via `beam.children`, and `point_on_beam` returns the point at a given distance along one beam. To analyse or animate the whole tree, [`path_segments`](@ref) flattens a traced [`Beam`](@ref), [`GaussianBeamlet`](@ref) (its chief ray), [`BeamletOptics.AstigmaticGaussianBeamlet`](@ref) or beam group into a vector of segments with start and end point, the accumulated geometric and optical path length, the wavelength, the position in the beam tree and the object hit at the end of the segment. A final ray without intersection has no length of its own and is drawn `flen` long.

```julia
solve_system!(system, beam)
segs = path_segments(beam; flen = 0.1)

# path length of the last intersection of each branch
for seg in segs
    seg.final && println("branch ends after ", seg.s_start, " m (OPL: ", seg.opl_start, " m)")
end
```

```@docs; canonical=false
path_segments
BeamletOptics.PathSegment
```
