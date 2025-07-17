# Optical components

A collection of basic optical elements is provided with this package as is. They are tested for the correctness of their optical interactions and are verified to work with reasonable fidelity. For detailed documentation, refer to the following table of contents. 

## Component overview

```@contents
Pages = ["mirrors.md", "lenses.md", "beamsplitters.md", "polarization_optics.md", "detectors.md"]
Depth = 2
```
## Listing available components

When using this package in the REPL, a tree view of all implemented [`BeamletOptics.AbstractObject`](@ref)s can be generated via the [`BeamletOptics.list_subtypes`](@ref) helper function. Note that this function is not able to determine all available constructors.

```@repl
using BeamletOptics # hide
BeamletOptics.list_subtypes(BeamletOptics.AbstractObject);
```