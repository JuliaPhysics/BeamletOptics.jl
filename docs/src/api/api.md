# API design

This part of the documentation is intended for users that want to change the internals of BMO. You can develop the package locally by typing `] dev BeamletOptics` in the REPL. An overview of all relevant topics is provided below. The following sections assume that you are in principle familiar with the concept of rays and beams, as well as optical elements -- e.g. mirrors -- in the context of this package. If not, it is recommended that you read the [Rays](@ref), [Beams](@ref) and [Optical components](@ref) sections first.

!!! warning
    The developer section is constantly changing as we iterate on the pre-1.0 BMO releases. If you find any outdated information, please open an issue!

```@contents
Pages = ["conventions.md", "core.md", "geometry.md", "meshes.md", "sdfs.md", "kinematics_api.md"]
Depth = 2
```

## Public API for dependent packages

Besides the exported names, a few names are declared `public` in `src/Exports.jl` without being exported, e.g. the progress interface `BeamletOptics.ProgressSink`, the sampling of the sources `BeamletOptics.AbstractSourceSampling` and the samplings of a detector `BeamletOptics.AbstractDetectorSampling` with `BeamletOptics.sample_points` and its result `BeamletOptics.DetectorSamples` (see [Sampling](@ref)), the optional `BeamletOptics.bounding_sphere_of` of a shape (see [Bounding sphere pretest](@ref)), `BeamletOptics.AbstractObjectGroup` with `BeamletOptics.check_kinematic_members` for own assembly types of separate parts (see [Assemblies: group or multi-shape object](@ref)), or `BeamletOptics.relaunch!` and the beamlet hit records (`BeamletOptics.AstigmaticGaussianBeamletHit(agb, id)`, `electric_field(hit, p)`) for components that re-emit light, e.g. solvers coupled through a field (see [Re-emitting components](@ref)). Packages built on BMO, such as [BeamletOpticsGUI](https://github.com/StackEnjoyer/BeamletOpticsGUI.jl), may rely on them like on the exported API; other unexported names are internal and may change in any release. Whether a name is public can be checked with `Base.ispublic(BeamletOptics, name)`. Their docstrings are listed in the [Reference](@ref), the render handle protocol is described on the [Live rendering](@ref "Render handle protocol") page. A new `AbstractObject` subtype can add a method of [`BeamletOptics.properties`](@ref) to list its own parameters (`properties` and `default_properties` are public, not exported: write `BeamletOptics.properties`), see [Inspecting components](@ref).
