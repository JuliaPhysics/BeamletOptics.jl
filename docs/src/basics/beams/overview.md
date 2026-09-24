# Beams

As mentioned in the [Rays](@ref) section, a beam within the context of this package serves as a data structure for storing collections of rays, forming the backbone of the simulation framework. Beams are intended to be designed as [AbstractTrees](https://github.com/JuliaCollections/AbstractTrees.jl) to allow for ray bifurcations, e.g. in the case of optical elements such as beamsplitters. The [`solve_system!`](@ref) function relies on this data structure to perform ray tracing computations within optical systems. 

To ensure compatibility and extensibility, beam types must adhere to the [`BeamletOptics.AbstractBeam`](@ref) interface. Refer to its documentation for more information. For detailed documentation of the provided beam types and sources, refer to the following table of contents.

## Beam overview

```@contents
Pages = ["beams.md", "stigmatic_beam.md", "astigmatic_beam.md", "beam_groups.md"]
Depth = 2
```

## Moving beams

Rays, beams and beam groups can be moved with the same kinematic API as optical elements, see [Moving sources](@ref) in the [Kinematics](@ref) section.

## Listing available beams

When using this package in the REPL, a tree view of all implemented [`BeamletOptics.AbstractBeam`](@ref)s and [`BeamletOptics.AbstractBeamGroup`](@ref)s can be generated via the [`BeamletOptics.list_subtypes`](@ref) helper function. Note that this function is not able to determine all available constructors.

```@repl
using BeamletOptics # hide
BeamletOptics.list_subtypes(BeamletOptics.AbstractBeam);
BeamletOptics.list_subtypes(BeamletOptics.AbstractBeamGroup);
```
