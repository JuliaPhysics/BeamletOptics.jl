# Implementing components

This page explains which functions a solve calls on an optical component and what a developer implements for a new [`BeamletOptics.AbstractObject`](@ref) subtype. The geometry of a component lives in its shape, its optical effect in its [`BeamletOptics.interact3d`](@ref) method; the solver connects both in the [Intersect-Interact-Repeat-Loop](@ref). For the components themselves, see [Optical components](@ref).

## Hooks called during a solve

The following diagram shows which functions a solve calls on an object. A new `AbstractObject` subtype implements `interact3d` always, and `initialize!`, `bounding_sphere_of` and `shape` only where the defaults or the existing shapes do not fit; it may also override the generic methods. `initialize!` runs only if the solve is called with `initialize = true`, for a beam group once before the first beam. How the loop chooses the next hit, including a `Hint`, is described in the [Intersect-Interact-Repeat-Loop](@ref):

```mermaid
flowchart TB
    S(["<b>solve_system!</b>"])
    INI["<b>initialize!(object)</b>"]
    BSO["<b>bounding_sphere_of(object)</b><br/>MultiShape: once per part"]
    LOOP["<b>Intersect-Interact-Repeat-Loop</b>"]
    IX["<b>intersect3d(object, ray)</b><br/><i>generic method</i>"]
    GB["<b>interact3d(sys, obj, beamlet, id)</b><br/><i>generic beamlet method</i>"]
    SHP["<b>shape(object), intersect3d(shape, ray)</b><br/>SingleShape (default): object.shape<br/>MultiShape (opt in via shape_trait_of):<br/>tuple of parts, tested the same way"]
    INT["<b>interact3d(sys, obj, beam, ray)</b>"]
    BACK(["<b>Repeat loop</b>"])

    S -- "initialize = true" --> INI --> BSO
    S -- "otherwise" --> BSO
    BSO --> LOOP
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

The `interact3d` methods of the component families decide what follows a hit: a new ray, child beams, a recorded hit, a hand-over to the part of a composite that was hit, or the end of the beam. The solver calls the method for a beam with the hit `ray`; for a beamlet, the generic method calls it for each component beam. The families are grouped by their result:

```mermaid
flowchart TB
    INT(["<b>interact3d</b><br/>system, object, beam, ray<br/><i>dispatch on the object type</i>"])
    CMP["<b>MultiShape composites</b><br/>E.g. DoubletLens, TripletLens:<br/>call interact3d of the part that was hit"]
    RAY["<b>Spawns new ray</b><br/>E.g. via refraction3d or reflection3d"]
    FIN["<b>Ends current beam</b><br/>E.g. for ThinBeamsplitter: <br/> attach children and terminate current"]
    DET["<b>Detector</b><br/>records a hit"]
    BI(["Return a <b>BeamInteraction</b><br/>for the next ray, optional <b>Hint</b>"])
    NO(["Return <b>nothing</b>"])

    INT -- "composite" --> CMP
    INT -- "single shape" --> RAY
    INT -- "single shape" --> FIN
    INT -- "single shape" --> DET
    INT -- "IntersectableObject" --> NO
    CMP -- "BeamInteraction(Hint(composite,<br/>next part), ray)" --> BI
    CMP -- "part returns nothing" --> NO
    RAY --> BI
    RAY -- "Blocking components: <br/>e.g. PolarizationFilter" --> NO
    FIN --> NO
    DET -- "stop = false:<br/>new ray from the hit point" --> BI
    DET -- "stop = true" --> NO

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class RAY,FIN,DET green
    class CMP purple
    class INT,BI,NO terminal
```

A composite calls the `interact3d` method of the part that was hit directly and does not run the solver again. It returns the `BeamInteraction` of that part with a `Hint` to the part to be tested next, or `nothing` if the part returns `nothing`. Why the hint is needed at coincident boundaries is explained in [Hints](@ref).
