# API design

This part of the documentation is intended for users that want to change the internals of BMO. You can develop the package locally by typing `] dev BeamletOptics` in the REPL. An overview of all relevant topics is provided below. The following sections assume that you are in principle familiar with the concept of rays and beams, as well as optical elements -- e.g. mirrors -- in the context of this package. If not, it is recommended that you read the [Rays](@ref), [Beams](@ref) and [Optical components](@ref) sections first.

!!! warning
    The developer section is constantly changing as we iterate on the pre-1.0 BMO releases. If you find any outdated information, please open an issue!

```@contents
Pages = ["conventions.md", "core.md", "geometry.md", "meshes.md", "sdfs.md", "kinematics_api.md"]
Depth = 2
```