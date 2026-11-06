```@setup beams
beam_showcase_dir = joinpath(@__DIR__, "..", "..", "assets", "beam_renders")

Main.DocUtils.conditional_include(joinpath(beam_showcase_dir, "beam_showcase.jl"))
```

# Basic beam

A minimal implementation of the [`BeamletOptics.AbstractBeam`](@ref) type is provided by the [`Beam`](@ref). It can be used to store a light path through an optical system. If the beam is split, its children will be recursively traced until all paths are solved.

```@docs; canonical=false
Beam
```

A [`Beam`](@ref) of [`Ray`](@ref)s is spawned from a position and a direction, a [`Beam`](@ref) of [`PolarizedRay`](@ref)s additionally from the electric field vector:

```@docs; canonical=false
Beam(::AbstractArray, ::AbstractArray, ::Any)
Beam(::AbstractArray, ::AbstractArray, ::Any, ::Any)
```

A ray tracing example through an arbitrary system using a [`Beam`](@ref) is shown below. Individual [`Ray`](@ref) segments are marked by their starting position and direction. The [Laser alignment](@ref) and [Miniature microscope](@ref) tutorial covers the use of the [`Beam`](@ref) in more detail. 

![Beam structure](beam_showcase.png)

## Inspecting a traced beam

The rays of a [`Beam`](@ref) are available via [`rays`](@ref) (for a beamlet, those of its chief ray), its child beams via `children(beam)` of AbstractTrees, and [`point_on_beam`](@ref) returns the point at a given distance along one beam. To analyse or animate the whole tree, [`path_segments`](@ref) flattens a traced [`Beam`](@ref), [`GaussianBeamlet`](@ref) (its chief ray), [`AstigmaticGaussianBeamlet`](@ref) or beam group into a vector of [`BeamletOptics.PathSegment`](@ref)s with start and end point, the accumulated geometric and optical path length, the wavelength, the position in the beam tree and the object hit at the end of the segment. A final ray without intersection has no length of its own and is drawn `flen` long. The last segment of each branch has `leaf = true`, whether it leaves the system as a final ray or ends on an object, e.g. a detector.

```julia
solve_system!(system, beam)
segs = path_segments(beam; flen = 0.1)

# path length at the last intersection of each branch
for seg in segs
    seg.leaf || continue
    s, opl = seg.final ? (seg.s_start, seg.opl_start) : (seg.s_stop, seg.opl_stop)
    println("branch ends after ", s, " m (OPL: ", opl, " m)")
end
```

```@docs; canonical=false
path_segments
BeamletOptics.PathSegment
```
