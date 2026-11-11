# Beams and sources

This page explains how the beams and sources described in the [Beams](@ref) section, including the beam groups, work internally, and what a developer has to provide so that a new beam type, beam group or source sampling works with the solver.

## Beam interface

A [`Beam`](@ref) stores a vector of rays, and the beamlets are built from it: a [`GaussianBeamlet`](@ref) traces a chief, a waist and a divergence `Beam` side by side, an [`AstigmaticGaussianBeamlet`](@ref) a chief `Beam` of [`PolarizedRay`](@ref)s and eight auxiliary ones. The diagram shows the beam types, the beams and rays they are made of, and the interface a new beam type implements. All three beam types are subtypes of `AbstractBeam` (frame); the beamlets store their rays in `Beam`s, so in the end every beam type holds its rays in the `Vector` of a `Beam`:

```mermaid
flowchart TB
    ABS["<b>AbstractBeam{T, R}</b><br/><i>abstract type, R #lt;: AbstractRay</i><br/>every subtype: parent, children<br/>of its own type, split by children!"]
    IF["<b>beam interface</b><br/>trace_system!, push!, and<br/>every beam: _component_beams<br/>ray-storing beam: also rays, empty!"]
    subgraph SUB ["subtypes of AbstractBeam"]
        GB["<b>GaussianBeamlet</b><br/><i>R = Ray</i>"]
        AGB["<b>AstigmaticGaussianBeamlet</b><br/><i>R = PolarizedRay</i>"]
        BM["<b>Beam{T, R}</b><br/>rays::Vector{R}"]
    end
    RAY["<b>AbstractRay</b><br/><i>Ray, PolarizedRay</i>"]

    ABS -.-> SUB
    ABS -- "subtype implements" --> IF
    GB -- "chief, waist, divergence:<br/>Beam of Rays" --> BM
    AGB -- "c: Beam of PolarizedRays<br/>wxp … dym: 8 Beams of Rays" --> BM
    BM -- "rays" --> RAY

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    class ABS,BM,GB,AGB blue
    class IF green
    class RAY purple
    style SUB fill:transparent,stroke:#4063D8,stroke-width:1px,stroke-dasharray:4 3
```

## Tracing a beamlet

Components define [`BeamletOptics.interact3d`](@ref) for a `Beam` and one of its rays; a beamlet hands this call on to each of its component beams, unless the component defines a method for the beamlet itself or is a thin layer, whose surface model handles the whole beamlet, see [Transitions](@ref). The diagram shows the tracing loop of a `GaussianBeamlet` and where the methods of a component author are called; an `AstigmaticGaussianBeamlet` runs the same loop with its nine beams and also checks its complex invariant:

```mermaid
flowchart TB
    S(["<b>trace_system!</b><br/>system, gauss"])
    T["<b>tracing_step!</b><br/>last ray of chief,<br/>then waist and divergence"]
    C["<b>check hits</b><br/>every ray hit,<br/>all the same shape?"]
    E1(["stop"])
    I["<b>interact3d(system, object, gauss, id)</b><br/><i>dispatch on object</i>"]
    OWN["<b>method for the beamlet</b><br/>e.g. Detector: push! hit"]
    MOD["<b>surface model of a thin layer</b><br/>e.g. BeamsplitterCoating: children!"]
    GEN["<b>generic beamlet method</b><br/>one call per component beam"]
    CB["<b>interact3d(system, object, beam, ray)</b><br/><i>method for a Beam</i>"]
    E2(["stop:<br/>an interaction is nothing"])
    GBI["<b>GaussianBeamletInteraction</b><br/>chief, waist, divergence"]
    P["<b>push!(gauss, interaction)</b><br/>push! to each component Beam"]
    K["<b>check_lagrange_invariants</b><br/>Lagrange invariant within<br/>get_lagrange_tolerance?"]
    E3(["stop with warning"])
    RM{"fewer than<br/>r_max rays?"}
    E4(["stop:<br/>r_max reached"])
    AGAIN(["next loop pass:<br/><b>tracing_step!</b>"])

    S --> T --> C
    C -- no --> E1
    C -- yes --> I
    I -- "method for<br/>the object" --> OWN
    I -- "thin layer" --> MOD
    I -- otherwise --> GEN
    GEN --> CB
    CB -- nothing --> E2
    OWN -- nothing --> E2
    CB -- "3 BeamInteractions" --> GBI
    OWN -- interaction --> GBI
    MOD -- nothing --> E2
    MOD -- interaction --> GBI
    GBI --> P --> K
    K -- no --> E3
    K -- yes --> RM
    RM -- yes --> AGAIN
    RM -- no --> E4

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class T,C,P,K,RM blue
    class I,GEN,GBI purple
    class OWN,CB,MOD green
    class S,E1,E2,E3,E4,AGAIN terminal
```

## Interface contracts

The docstrings of the abstract types list what a subtype implements: fields, getters and functions.

```@docs; canonical=false
BeamletOptics.AbstractBeam
BeamletOptics.AbstractBeamGroup
BeamletOptics.AbstractSourceSampling
```

## Source sampling

A [`CollimatedSource`](@ref) or [`PointSource`](@ref) stores how its rays were sampled as a [`BeamletOptics.AbstractSourceSampling`](@ref), such that [`num_rays!`](@ref) can regenerate them with another number of rays in the current pose of the source. The diagram shows how the constructors and `num_rays!` generate the beams with [`BeamletOptics.source_beams`](@ref): dotted arrows are its dispatch branches on the sampling type.

```mermaid
flowchart TB
    C(["<b>catalog constructor</b><br/>e.g. CollimatedSource(pos, dir, diameter, λ)"])
    GC(["<b>generic constructor</b><br/>CollimatedSource(..., sampling, λ; num_rays)<br/>PointSource(..., sampling, λ; num_rays)"])
    SB["<b>sampling_basis(dir, basis, T)</b><br/>reference vector b1"]
    N(["<b>num_rays!(src, n)</b>"])
    SBM["<b>source_beams(sampling, pos, dir, b1, ...)</b><br/>one method per sampling type"]
    SAMP["<b>samplings of BMO</b><br/>CollimatedSource: DiscRings,<br/>DiscSunflower, DiscLine<br/>PointSource: ConeRings,<br/>ConeSunflower, ConeFan"]
    OWN["<b>own AbstractSourceSampling</b><br/>source_beams, min_num_rays (optional),<br/>numerical_aperture (PointSource)"]
    ERR(["ArgumentError"])
    BM["<b>Vector{Beam}</b><br/>first beam: center ray"]
    W(["<b>given beams</b><br/>e.g. CollimatedSource(beams, diameter, pos, dir)"])
    G(["<b>CollimatedSource / PointSource</b><br/>beams, center, dir, sampling<br/>PointSource: NA from<br/>numerical_aperture(sampling)"])

    C -- "its sampling, e.g. DiscRings(num_rings)" --> GC
    GC --> SB --> SBM
    N -- "src.sampling, orientation(src), n,<br/>emission point of the old rays" --> SBM
    SBM -.-> SAMP
    SBM -.-> OWN
    SBM -. NoSampling .-> ERR
    SAMP -- "num_rays below<br/>min_num_rays" --> ERR
    SAMP --> BM
    OWN --> BM
    BM -- "new source, or num_rays!<br/>replaces src.beams" --> G
    W -- "sampling = NoSampling()" --> G

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class SB,SBM,BM blue
    class OWN green
    class SAMP purple
    class C,GC,N,W,ERR,G terminal
```

An own sampling is a subtype of [`BeamletOptics.AbstractSourceSampling`](@ref) with a method of `source_beams`; the generic constructors and [`num_rays!`](@ref) then work for it without further code. A sampling for a `PointSource` also implements `numerical_aperture(s)`. A sampling with a lower bound on the number of rays implements `min_num_rays(s)`, and its `source_beams` throws an `ArgumentError` for a `num_rays` below it, as the samplings of BMO do.

## Solving a beam group internally

The diagram shows the steps of `solve_system!(system, bg)`. A new [`BeamletOptics.AbstractBeamGroup`](@ref) subtype provides the getter of its beams, as a `beams` field or a method.

```mermaid
flowchart TB
    S(["<b>solve_system!(system, bg)</b>"])
    EM["<b>empty!(bg)</b><br/>every beam to its start state"]
    INI{"initialize?"}
    IN["<b>initialize!(system)</b><br/>once, before the first beam"]
    BS["<b>with_bounding_spheres(system)</b><br/>one table for all beams"]
    BG["<b>beams(bg)</b><br/>root beams of the group,<br/>split into chunks: at most<br/>4 per thread"]
    PR["<b>_with_progress</b><br/>progress bar, if progress = true"]
    CH["<b>Threads.@threads :greedy</b><br/>over the chunks"]
    BT["<b>with_bounding_spheres(table)</b><br/>the table set again in each task"]
    HB["<b>with_hit_buffers</b><br/>detector hits of a chunk<br/>stored at its end"]
    SB["<b>solve_system!(system, beam; ...)</b><br/>for each beam of the chunk"]
    E(["return nothing"])

    S --> EM --> INI
    INI -- yes --> IN --> BS
    INI -- no --> BS
    BS --> BG --> PR --> CH --> BT --> HB --> SB --> E

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class EM,IN,BS,PR,CH,BT,HB,SB blue
    class BG green
    class S,E terminal
```
