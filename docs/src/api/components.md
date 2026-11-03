# Implementing components

This page explains which functions a solve calls on an optical component and what a developer implements for a new [`BeamletOptics.AbstractObject`](@ref) subtype. The geometry of a component lives in its shape, its optical effect in its [`BeamletOptics.interact3d`](@ref) method; the solver connects both in the [Intersect-Interact-Repeat-Loop](@ref). For the components themselves, see [Optical components](@ref).

## Hooks called during a solve

The following diagram shows which functions a solve calls on an object: a new `AbstractObject` subtype implements the green ones, `interact3d` always and the others only where the defaults or the existing shapes do not fit, and may override the generic methods in purple:

```mermaid
flowchart TB
    S(["<b>solve_system!</b><br/>system, beam"])
    INI["<b>initialize!</b>(object)<br/><i>default: nothing to reset</i>"]
    BSO["<b>bounding_sphere_of</b>(shape)<br/><i>default: NoBoundingSphere</i>"]
    TR["<b>trace_system!</b><br/>system, beam"]
    STEP["<b>tracing_step!</b><br/>system, ray, hint"]
    IX["<b>intersect3d</b>(object, ray)<br/><i>generic method</i>"]
    SHP["<b>shape</b>(object)<br/>SingleShape (default): object.shape<br/>MultiShape (opt in via shape_trait_of):<br/>tuple of parts"]
    ISH["<b>intersect3d</b>(shape, ray)<br/><i>after the bounding sphere test</i>"]
    GB["<b>interact3d</b><br/>(system, object, beamlet, id)<br/><i>generic beamlet method</i>"]
    INT["<b>interact3d</b><br/>(system, object, beam, ray)"]
    E([beam done])

    S -- "initialize = true" --> INI
    S -- bounding_spheres --> BSO
    S -- "root beam,<br/>then its children" --> TR
    TR --> STEP
    TR -- "r_max reached, or<br/>beamlet:<br/>invariant check fails" --> E
    STEP -- "nothing" --> IX
    STEP -- "Hint" --> ISH
    IX -- shape_trait_of --> SHP
    SHP -- "shape" --> ISH
    SHP -- "object part" --> IX
    ISH -- "Intersection<br/>(Beam)" --> INT
    ISH -- "Intersection<br/>(beamlet)" --> GB
    GB -- "each component beam" --> INT
    ISH -- "hinted shape<br/>missed" --> IX
    ISH -- "no object hit" --> E
    ISH -- "beamlet rays hit<br/>different shapes" --> E
    INT -- "BeamInteraction:<br/>push!, hint" --> STEP
    INT -- "nothing" --> E

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class TR,STEP blue
    class INI,BSO,SHP,ISH,INT green
    class IX,GB purple
    class S,E terminal
```

## Interaction of the component families

The `interact3d` methods of the component families decide what follows a hit: a new ray, child beams, a recorded hit, a hand-over to the part of a composite that was hit, or the end of the beam:

```mermaid
flowchart TB
    INT(["<b>interact3d</b><br/>system, object, beam, ray"])
    CMP["<b>MultiShape composites</b><br/>DoubletLens, TripletLens,<br/>CubeBeamsplitter,<br/>AbstractPlateBeamsplitter,<br/>LinearPolarizer"]
    DET["<b>Detector</b><br/>push!(detector, hit)"]
    RAY["<b>new ray</b><br/>AbstractRefractiveOptic: refraction3d<br/>AbstractReflectiveOptic: reflection3d<br/>PolarizationFilter: Jones matrix"]
    SPL["<b>ThinBeamsplitter</b><br/>children!: transmitted, reflected"]
    BI(["<b>BeamInteraction</b><br/>next ray, optional Hint"])
    HT["<b>IntersectableObject</b><br/>hard target"]
    NO(["<b>nothing</b><br/>this beam ends"])

    INT --> CMP
    INT --> DET
    INT --> RAY
    INT --> SPL
    INT --> HT
    INT -- "no method:<br/>warning" --> NO
    CMP -- "part that was hit,<br/>optional Hint to the composite" --> RAY
    CMP -- coating --> SPL
    RAY --> BI
    RAY -- "PolarizationFilter:<br/>blocked" --> NO
    SPL -- "children are<br/>traced next" --> NO
    DET -- "stop = false:<br/>pass_through, same ray<br/>from the hit point" --> BI
    DET -- "stop = true" --> NO
    HT --> NO

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class CMP,DET,RAY,SPL,HT green
    class INT,BI,NO terminal
```
