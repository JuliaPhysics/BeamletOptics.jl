```@setup beam_catalog
beam_catalog_dir = joinpath(@__DIR__, "..", "..", "assets", "beam_renders")

Main.DocUtils.conditional_include(joinpath(beam_catalog_dir, "beam_catalog_showcase.jl"))
```

# Beams

As mentioned in the [Rays](@ref) section, a beam within the context of this package serves as a data structure for storing collections of rays, forming the backbone of the simulation framework. Beams are intended to be designed as [AbstractTrees](https://github.com/JuliaCollections/AbstractTrees.jl) to allow for ray bifurcations, e.g. in the case of optical elements such as beamsplitters. The [`solve_system!`](@ref) function relies on this data structure to perform ray tracing computations within optical systems. 

To ensure compatibility and extensibility, beam types must adhere to the [`BeamletOptics.AbstractBeam`](@ref) interface. Refer to its documentation for more information. Browse the catalog below, or refer to the table of contents for detailed documentation of the provided beam types and sources.

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

A [`Beam`](@ref) stores a vector of rays, and the beamlets are built from it: a [`GaussianBeamlet`](@ref) traces a chief, a waist and a divergence `Beam` side by side, an [`AstigmaticGaussianBeamlet`](@ref) a chief `Beam` of [`PolarizedRay`](@ref)s and eight auxiliary ones. The diagram shows the beam types, the beams and rays they are made of, and the interface a new beam type implements (green):

```mermaid
flowchart TB
    ABS["<b>AbstractBeam</b><br/><i>abstract type</i>"]
    IF["<b>beam interface</b><br/>trace_system!, push!, and<br/>every beam: _component_beams<br/>ray-storing beam: also rays, empty!"]
    TREE["<b>parent, children</b><br/><i>beam tree, split by children!</i>"]
    BM["<b>Beam</b><br/><i>beam type</i>"]
    GB["<b>GaussianBeamlet</b><br/><i>beam type</i>"]
    AGB["<b>AstigmaticGaussianBeamlet</b><br/><i>beam type</i>"]
    RAY["<b>AbstractRay</b><br/><i>Ray, PolarizedRay</i>"]

    ABS -- "subtype implements" --> IF
    ABS -- fields --> TREE
    ABS -.-> BM
    ABS -.-> GB
    ABS -.-> AGB
    GB -- "chief, waist, divergence" --> BM
    AGB -- "chief c and 8 auxiliary<br/>(waist, divergence in x, y)" --> BM
    BM -- rays --> RAY

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    class ABS,BM,GB,AGB blue
    class IF green
    class TREE,RAY purple
```

### Interactions of beamlets

Components define [`BeamletOptics.interact3d`](@ref) for a `Beam` and one of its rays; a beamlet hands this call on to each of its component beams, unless the component defines a method for the beamlet itself. The diagram shows the tracing loop of a `GaussianBeamlet` and where the methods of a component author (green) are called; an `AstigmaticGaussianBeamlet` runs the same loop with its nine beams and also checks its complex invariant:

```mermaid
flowchart TB
    S(["<b>trace_system!</b><br/>system, gauss"])
    T["<b>tracing_step!</b><br/>last ray of chief,<br/>then waist and divergence"]
    C{"all rays hit<br/>the same shape?"}
    I{"<b>interact3d</b><br/>(system, object, gauss, id)<br/><i>dispatch on object</i>"}
    OWN["<b>method for the beamlet</b><br/>e.g. ThinBeamsplitter: children!<br/>Detector: push! hit"]
    GEN["<b>generic beamlet method</b><br/>one call per component beam"]
    CB["<b>interact3d</b><br/>(system, object, beam, ray)<br/><i>method for a Beam</i>"]
    GBI["<b>GaussianBeamletInteraction</b><br/>chief, waist, divergence"]
    P["<b>push!(gauss, interaction)</b><br/>push! to each component Beam"]
    K{"check_invariant:<br/>Lagrange invariant within<br/>get_lagrange_tolerance?"}
    E([stop])

    S --> T
    T -- "a ray hits nothing" --> E
    T -- hit --> C
    C -- no --> E
    C -- yes --> I
    I -- "method for<br/>the object" --> OWN
    I -- otherwise --> GEN
    GEN --> CB
    CB -- "3 BeamInteractions" --> GBI
    CB -- "any nothing" --> E
    OWN -- interaction --> GBI
    OWN -- nothing --> E
    GBI --> P
    P --> K
    K -- "no: warning" --> E
    K -- "yes: next segment<br/>(fewer than r_max rays)" --> T

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class T,C,P,K blue
    class I,GEN,GBI purple
    class OWN,CB green
    class S,E terminal
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
