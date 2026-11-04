# Beams and sources

This page explains how the beams and sources described in the [Beams](@ref) section, including the beam groups, work internally, and what a developer has to provide so that a new beam type, beam group or source sampling works with the solver.

## Beam interface

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

## Tracing a beamlet

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

## Source sampling

A [`CollimatedSource`](@ref) or [`PointSource`](@ref) stores how its rays were sampled as a [`BeamletOptics.AbstractSourceSampling`](@ref), such that [`num_rays!`](@ref) can regenerate them with another number of rays in the current pose of the source. The diagram shows how the constructors and `num_rays!` generate the beams with [`BeamletOptics.source_beams`](@ref): dotted arrows are its dispatch branches on the sampling type, the green box is what an own sampling type implements.

```mermaid
flowchart TB
    C(["<b>catalog constructor</b><br/>e.g. CollimatedSource(pos, dir, diameter, λ)"])
    GC(["<b>generic constructor</b><br/>CollimatedSource(pos, dir, diameter, sampling, λ; num_rays)<br/>PointSource(pos, dir, sampling, λ; num_rays)"])
    N(["<b>num_rays!</b>(src, n)"])
    W(["<b>given beams</b><br/>e.g. CollimatedSource(beams, diameter, pos, dir)"])
    SB["<b>sampling_basis</b>(dir, basis, T)<br/>reference vector b1"]
    SBM["<b>source_beams</b>(sampling, pos, dir, b1,<br/>args..., λ, num_rays, T)"]
    DISC["<b>DiscRings, DiscSunflower, DiscLine</b><br/><i>CollimatedSource</i>"]
    CONE["<b>ConeRings, ConeSunflower, ConeFan</b><br/><i>PointSource</i>"]
    OWN["<b>own AbstractSourceSampling</b><br/>methods of source_beams, min_num_rays (optional)<br/>and numerical_aperture (PointSource)"]
    NO["<b>NoSampling</b>"]
    ERR(["ArgumentError"])
    BM["<b>Vector{Beam}</b><br/>first beam: center ray"]
    NA["<b>numerical_aperture</b>(sampling)<br/>NA of a PointSource"]
    G["<b>CollimatedSource / PointSource</b><br/>beams, center, dir, sampling"]
    MIN["<b>min_num_rays</b>(src.sampling)<br/>20 num_rings, 1 or nothing"]

    C -- "its sampling, e.g. DiscRings(num_rings)" --> GC
    GC --> SB
    SB --> SBM
    N -- "src.sampling, orientation(src), n" --> SBM
    SBM -.-> DISC
    SBM -.-> CONE
    SBM -.-> OWN
    SBM -.-> NO
    DISC --> BM
    CONE --> BM
    OWN --> BM
    DISC -- "num_rays < min_num_rays" --> ERR
    CONE -- "num_rays < min_num_rays" --> ERR
    NO --> ERR
    BM -- "replaces src.beams" --> G
    BM -- "new source" --> NA --> G
    W -- "sampling = NoSampling()" --> G
    G -- "min_num_rays(src)" --> MIN

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class SB,SBM,MIN,NA blue
    class OWN green
    class DISC,CONE,NO,BM,G purple
    class C,GC,N,W,ERR terminal
```

An own sampling is a subtype of [`BeamletOptics.AbstractSourceSampling`](@ref) with a method of `source_beams`; the generic constructors and [`num_rays!`](@ref) then work for it without further code. A sampling for a `PointSource` also implements `numerical_aperture(s)`. A sampling with a lower bound on the number of rays implements `min_num_rays(s)`, and its `source_beams` throws an `ArgumentError` for a `num_rays` below it, as the samplings of BMO do.

## Solving a beam group internally

The diagram shows the steps of `solve_system!(system, bg)`, see [Solving a beam group](@ref); the green box is the getter of the beams that a new [`BeamletOptics.AbstractBeamGroup`](@ref) subtype provides, as a `beams` field or a method.

```mermaid
flowchart TB
    S(["<b>solve_system!</b>(system, bg)"])
    BG["<b>beams</b>(bg)<br/>root beams of the group"]
    EM["<b>empty!</b>(bg)<br/>every beam to its start state"]
    INI{"initialize = true?"}
    IN["<b>initialize!</b>(system)<br/>once, before the first beam"]
    BS["<b>with_bounding_spheres</b>(system)<br/>one table for all beams"]
    CH["<b>Threads.@threads :greedy</b><br/>over chunks of beams(bg)"]
    HB["<b>with_hit_buffers</b><br/>detector hits of a chunk<br/>stored at its end"]
    SB["<b>solve_system!</b>(system, beam; kwargs...)<br/>for each beam of the chunk"]
    E([done])

    S --> EM --> INI
    INI -- yes --> IN --> BS
    INI -- no --> BS
    BS --> CH --> HB --> SB --> E
    BG -.-> EM
    BG -.-> CH

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class EM,IN,BS,CH,HB,SB blue
    class BG green
    class S,E terminal
```
