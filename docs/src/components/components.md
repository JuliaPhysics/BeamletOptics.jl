# Optical components

A collection of basic optical elements is provided with this package as is. They are tested for the correctness of their optical interactions and are verified to work with reasonable fidelity. For detailed documentation, refer to the following sections ([Component overview](@ref)). 

## Component overview

```@contents
Pages = ["mirrors.md", "lenses.md", "beamsplitters.md", "photodetectors.md"]
Depth = 2
```
## Component list

When using this package in the REPL, a tree view of all implemented [`SCDI.AbstractObject`](@ref)s can be generated via the [`SCDI.list_subtypes`](@ref) helper function. Note that this function is not able to determine all available constructors.

```@repl
using SCDI # hide
SCDI.list_subtypes(SCDI.AbstractObject)
```