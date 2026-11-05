# Implementing components

This page explains which functions a solve calls on an optical component and what a developer implements for a new [`BeamletOptics.AbstractObject`](@ref) subtype. The geometry of a component lives in its shape, its optical effect in its [`BeamletOptics.interact3d`](@ref) method; the solver connects both in the [Intersect-Interact-Repeat-Loop](@ref). For the components themselves, see [Optical components](@ref).

## Hooks called during a solve

The following diagram shows which functions a solve calls on an object: a new `AbstractObject` subtype implements the green ones, `interact3d` always and the others only where the defaults or the existing shapes do not fit, and may override the generic methods in purple; the loop itself is blue. How the loop chooses the next hit, including a `Hint`, is described in the [Intersect-Interact-Repeat-Loop](@ref):

```mermaid
flowchart TB
    S(["<b>solve_system!</b>"])
    INI["<b>initialize!(object)</b>"]
    BSO["<b>bounding_sphere_of(shape)</b>"]
    LOOP["<b>Intersect-Interact-Repeat-Loop</b>"]
    IX["<b>intersect3d(object, ray)</b><br/><i>generic method</i>"]
    GB["<b>interact3d(sys, obj, beamlet)</b><br/><i>generic beamlet method</i>"]
    SHP["<b>shape(object), intersect3d(shape, ray)</b><br/>SingleShape (default): object.shape<br/>MultiShape (opt in via shape_trait_of):<br/>tuple of parts, tested the same way"]
    INT["<b>interact3d</b>"]
    BACK(["<b>Repeat loop</b>"])

    S --> INI --> BSO --> LOOP
    LOOP --> IX
    LOOP -- "Beamlet" --> GB
    LOOP -- "Beam" --> INT
    IX -- shape_trait_of --> SHP
    GB -- "each component Beam" --> INT
    INT --> BACK

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class LOOP blue
    class INI,BSO,SHP,INT green
    class IX,GB purple
    class S,BACK terminal
```

## Interaction of the component families

The `interact3d` methods of the component families decide what follows a hit: a new ray, child beams, a recorded hit, a hand-over to the part of a composite that was hit, or the end of the beam. The families are grouped by their result (green: methods of the component families, purple: hand-over of a composite, red: entry point and results):

```mermaid
flowchart TB
    INT(["<b>interact3d</b><br/>system, object, beam, ray"])
    CMP["<b>MultiShape composites</b><br/>Custom interact3d code e.g.DoubletLens"]
    RAY["<b>Spawns new ray</b><br/>E.g. via refraction3d or reflection3d"]
    FIN["<b>Ends current beam</b><br/>E.g. for ThinBeamsplitter: <br/> attach children and terminate current"]
    AGAIN(["Repeat solver and <b>interact3d</b><br/>until part that was hit<br/>is a <b>SingleShape</b>"])
    BI(["Return a <b>BeamInteraction</b><br/> for next ray, optional <b>Hint</b>"])
    NO(["Return <b>nothing</b>"])

    INT -- isa <b>MultiShape</b> --> CMP
    INT -- isa <b>SingleShape</b> --> RAY
    INT -- isa <b>SingleShape</b> --> FIN
    CMP -- "optional Hint<br/>to the composite" --> AGAIN
    RAY --> BI
    RAY -- "Blocking components: <br/>e.g. PolarizationFilter" --> NO
    FIN --> NO

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class RAY,FIN green
    class CMP purple
    class INT,AGAIN,BI,NO terminal
```
