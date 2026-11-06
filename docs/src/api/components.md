# Implementing components

This page explains which functions a solve calls on an optical component and what a developer implements for a new [`BeamletOptics.AbstractObject`](@ref) subtype. The geometry of a component lives in its shape, its optical effect at a boundary in its [`BeamletOptics.medium`](@ref) (a body) or [`BeamletOptics.surface_model`](@ref) (a thin layer), or in its own [`BeamletOptics.interact3d`](@ref) method; the solver connects both in the [Intersect-Interact-Repeat-Loop](@ref). For the components themselves, see [Optical components](@ref).

## Hooks called during a solve

The following diagram shows which functions a solve calls on an object. A new `AbstractObject` subtype implements `interact3d`, unless it is a body or a thin layer, which declares its `medium` or `surface_model` instead (see [Transitions](@ref)), and `initialize!`, `bounding_sphere_of` and `shape` only where the defaults or the existing shapes do not fit; it may also override the generic methods. `initialize!` runs only if the solve is called with `initialize = true`, for a beam group once before the first beam. How the loop chooses the next hit, including a `Hint`, is described in the [Intersect-Interact-Repeat-Loop](@ref):

```mermaid
flowchart TB
    S(["<b>solve_system!</b>"])
    INI["<b>initialize!(object)</b>"]
    BSO["<b>bounding_sphere_of(object)</b><br/>MultiShape: once per part"]
    LOOP["<b>Intersect-Interact-Repeat-Loop</b>"]
    IX["<b>intersect3d(object, ray)</b><br/><i>generic method</i>"]
    GB["<b>interact3d(sys, obj, beamlet, id)</b><br/><i>generic beamlet method</i><br/>thin layer: its surface model"]
    SHP["<b>shape(object), intersect3d(shape, ray)</b><br/>SingleShape (default): object.shape<br/>MultiShape (opt in via shape_trait_of):<br/>tuple of parts, tested the same way"]
    INT["<b>interact3d(sys, obj, beam, ray)</b><br/>or <b>medium</b> / <b>surface_model</b>"]
    BACK(["<b>Repeat loop</b>"])

    S -- "initialize = true" --> INI --> BSO
    S -- "otherwise" --> BSO
    BSO --> LOOP
    LOOP --> IX
    LOOP -- "Beamlet" --> GB
    LOOP -- "Beam" --> INT
    IX -- shape_trait_of --> SHP
    GB -- "otherwise: each component Beam" --> INT
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

The `interact3d` methods of the component families decide what follows a hit: a new ray, child beams, a recorded hit, a hand-over to the part of a composite that was hit, or the end of the beam. A body or a thin layer needs no own method: the fallback crosses the boundary through its [`BeamletOptics.surface_model`](@ref). The solver calls the method for a beam with the hit `ray`; for a beamlet, the generic method calls it for each component beam. The families are grouped by their result:

```mermaid
flowchart TB
    INT(["<b>interact3d</b><br/>system, object, beam, ray<br/><i>dispatch on the object type</i>"])
    CMP["<b>MultiShape composites</b><br/>E.g. DoubletLens, CubeBeamsplitter:<br/>interact3d of the part that was hit"]
    TR["<b>Body or thin layer</b><br/>declares medium or surface_model:<br/>transition, surface_model(t)"]
    RAY["<b>Spawns new ray</b><br/>E.g. a mirror via reflection3d"]
    FIN["<b>Ends current beam</b><br/>E.g. BeamsplitterCoating: <br/> attach children and terminate current"]
    DET["<b>Detector</b><br/>records a hit"]
    BI(["Return a <b>BeamInteraction</b><br/>for the next ray, optional <b>Hint</b>"])
    NO(["Return <b>nothing</b>"])

    INT -- "composite" --> CMP
    INT -- "no own method" --> TR
    INT -- "own method" --> RAY
    INT -- "own method" --> DET
    INT -- "IntersectableObject" --> NO
    CMP -- "result of the part" --> BI
    CMP -- "part returns nothing" --> NO
    TR -- "Uncoated: refraction, TIR" --> BI
    TR -- "a model that splits" --> FIN
    TR -- "blocked, e.g. JonesFilm" --> NO
    RAY --> BI
    FIN --> NO
    DET -- "stop = false: passes as<br/>a thin layer without a model" --> BI
    DET -- "stop = true" --> NO

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class RAY,DET green
    class TR,FIN blue
    class CMP purple
    class INT,BI,NO terminal
```

A composite calls the `interact3d` method of the part that was hit directly and does not run the solver again; it returns the result of that part. Parts that touch need no `Hint`: the solver finds the media on both sides of a shared face, see [Coincident hits](@ref) and [Transitions](@ref).
