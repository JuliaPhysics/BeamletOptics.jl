```@setup beam_catalog
beam_catalog_dir = joinpath(@__DIR__, "..", "..", "assets", "beam_renders")

Main.DocUtils.conditional_include(joinpath(beam_catalog_dir, "beam_catalog_showcase.jl"))
```

# Beams

As mentioned in the [Rays](@ref) section, a beam within the context of this package serves as a data structure for storing collections of rays, forming the backbone of the simulation framework. A beam splits into child beams, e.g. at a beamsplitter, and [`solve_system!`](@ref) traces it together with all its children. Browse the catalog below, or refer to the table of contents for detailed documentation of the provided beam types and sources.

## Beam catalog

Beams, beamlets and beam groups are traced through a bi-convex lens onto a detector near its focus.

```@raw html
<ComponentCatalog catalog="beams" />
```

## Beam overview

```@contents
Pages = ["beams.md", "stigmatic_beam.md", "astigmatic_beam.md", "beam_groups.md"]
Depth = 2
```

## Beam models

The beamlets are traced as several [`Beam`](@ref)s of rays side by side (blue: beam types, purple: rays):

```mermaid
flowchart TB
    BM["<b>Beam</b><br/><i>beam type</i>"]
    GB["<b>GaussianBeamlet</b><br/><i>beam type</i>"]
    AGB["<b>AstigmaticGaussianBeamlet</b><br/><i>beam type</i>"]
    RAY["<b>Ray, PolarizedRay</b><br/><i>rays</i>"]

    GB -- "3 beams:<br/>chief, waist, divergence" --> BM
    AGB -- "9 beams:<br/>chief and 8 auxiliary" --> BM
    BM -- rays --> RAY

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    class BM,GB,AGB blue
    class RAY purple
```

A beamlet stops when its rays hit different shapes or its invariant check fails. The tracing loop and the interface of a beam type are described in [Tracing a beamlet](@ref) and [Beam interface](@ref).

## Moving beams

Rays, beams and beam groups can be moved with the same kinematic API as optical elements, see [Moving sources](@ref) in the [Kinematics](@ref) section.
