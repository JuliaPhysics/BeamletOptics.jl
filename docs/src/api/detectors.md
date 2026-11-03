# Detector internals

This page explains how the [`Detector`](@ref) described in the basics is implemented internally, and what a developer has to provide so that a new hit type or sampling works with the readouts. What a readout computes is defined in [Reading out a detector](@ref) and [Field distributions](@ref).

## Hit recording

The following diagram shows how a solve records one hit type per beam model in a `Detector` and what empties it again; green boxes are the methods a new [`BeamletOptics.AbstractDetector`](@ref) implements.

```mermaid
flowchart TB
    S(["<b>solve_system!</b><br/>beam or beam group"]) --> IA["<b>interact3d</b>(system, d, beam, ray or id)<br/>one method per beam model"]
    IA -- Ray --> RH["<b>RayHit</b>"]
    IA -- PolarizedRay --> PH["<b>PolarizedRayHit</b>"]
    IA -- GaussianBeamlet --> GH["<b>GaussianBeamletHit</b>"]
    IA -- AstigmaticGaussianBeamlet --> AH["<b>AstigmaticGaussianBeamletHit</b>"]
    RH --> P["<b>push!</b>(d, hit)"]
    PH --> P
    GH --> P
    AH --> P
    P -- "solve of a beam group" --> B["<b>with_hit_buffers</b><br/>one buffer per chunk of beams and detector"]
    P -- "otherwise, under d.lock" --> H["<b>hits</b>(d)<br/>nothing or a Vector of one hit type,<br/>counted by length(d)"]
    B -- "append!(d, buffer)<br/>at the end of the chunk" --> H
    I(["<b>initialize!</b>(system)<br/>or initialize = true"]) -- "initialize!(d)" --> E["<b>empty!</b>(d)"]
    K(["<b>translate3d!</b> / <b>rotate3d!</b><br/>of d, also within an ObjectGroup"]) --> E
    E -- "hits = nothing" --> H

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class IA,E green
    class RH,PH,GH,AH purple
    class P,B,H blue
    class S,I,K terminal
```

## Readout pipeline

The following diagram shows how a readout turns a sampling into points and evaluates the stored hits there; green boxes are the methods a new sampling or hit type implements. A new hit type must also be added to the `hits` field of [`Detector`](@ref).

```mermaid
flowchart TB
    R(["<b>electric_field</b> / <b>intensity</b> / <b>optical_power</b><br/>(pd, sampling = RegularGrid())"]) --> SP["<b>sample_points</b>(pd, sampling)"]
    SP -- RegularGrid --> RG["<b>RegularGrid</b><br/>window around the hits via calc_local_lims,<br/>active surface without hits"]
    SP -- RectilinearGrid --> RL["<b>RectilinearGrid</b><br/>points xs × zs,<br/>trapezoidal weights"]
    SP -- "own AbstractDetectorSampling" --> OWN["own <b>sample_points</b> method"]
    RG -- "own hit type, unless an<br/>AbstractRayHit or AbstractBeamletHit" --> LIM["own <b>calc_local_lims</b>(pd, hits, crop_factor) methods"]
    RG -- "RectilinearGrid(xs, zs)" --> RL
    RL --> DS["<b>DetectorSamples</b><br/>coordinates, points, weights"]
    OWN --> DS
    DS -- points --> PT["<b>electric_field</b>(pd, points)<br/><b>intensity</b>(pd, points)<br/>threaded over the points"]
    PT -- "electric_field" --> SUM["<b>electric_field</b>(hits(pd), p)<br/>coherent sum over the hits"]
    PT -- "intensity" --> IH["<b>intensity</b>(hits(pd), p, n̂)<br/>Poynting flux through the detector"]
    IH -- "PolarizedRayHit,<br/>AstigmaticGaussianBeamletHit" --> IV["vector form<br/>½ Re[(E × H*) · n̂]"]
    IH -- "RayHit,<br/>GaussianBeamletHit" --> IS["scalar form<br/>½ Re[E* F]"]
    IH -- "own hit type" --> OWNI["own <b>intensity</b>(hits, p, n̂) method"]
    PT -- "no hits: hits(pd) = nothing" --> ZERO["zero field and flux<br/>at every point"]
    SUM --> EF["<b>electric_field</b>(hit, p)<br/>direction(hit), refractive_index(hit)"]
    IV --> EF
    IS --> EF
    PT -- "values" --> RET(["electric_field, intensity:<br/>(coordinates..., values)"])
    PT -- "intensity values I" --> OP(["optical_power:<br/>sum of weights .* I"])

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class SP,RG,RL,DS purple
    class PT,SUM,IH,IV,IS,ZERO blue
    class OWN,LIM,OWNI,EF green
    class R,RET,OP terminal
```

```@docs; canonical=false
BeamletOptics.AbstractDetectorSampling
BeamletOptics.sample_points
BeamletOptics.DetectorSamples
```

## OpticsBase extension internals

The diagram shows how the constructors share one sampling of the hits and how the return path leaves the `PlaneField` again. Blue boxes are BMO's detector internals that the extension reuses, green boxes depend on the hit type, purple boxes are the rest of the extension.

```mermaid
flowchart TB
    DET(["<b>PlaneField</b>(detector; size, spacing)"]) --> DH["<b>hits</b>(detector)<br/>ArgumentError if the detector has no hits"]
    DET -- "default axes" --> DA["<b>default_plane_axes</b>(detector)<br/>u from <b>local_axes</b>, n = −<b>direction</b>"]
    BEAM(["<b>PlaneField</b>(beamlet, id; size, spacing,<br/>origin, axes)"]) --> OH["one <b>GaussianBeamletHit</b> or<br/><b>AstigmaticGaussianBeamletHit</b><br/>(beamlet, id)"]
    DH --> PH
    DA -- "axes" --> PH
    OH --> PH
    PH["<b>PlaneField</b>(hits; size, spacing, origin, axes)<br/>beamlet hits only, otherwise ArgumentError"]
    PH --> CW["<b>common_wavelength_index</b>(hits)<br/>all hits share λ and n"]
    CW --> HFD["<b>hit_field_direction</b>(hit, u, v)<br/>per hit type"]
    HFD -- AstigmaticGaussianBeamletHit --> HA["nothing: <b>electric_field</b>(hit, p)<br/>is a 3D vector"]
    HFD -- GaussianBeamletHit --> HG["u made transverse to d<br/>(v if u ∥ d)"]
    HA --> SM["<b>_sample</b>(hits, origin, u, v, points)<br/>threaded over points, shared with the readouts"]
    HG --> SM
    SM -- "per point p" --> HF["<b>hit_field</b>(hit, p, pol)<br/>per hit type"]
    HF --> EH["E = Σᵢ Eᵢ projected on u, v<br/>H = Σᵢ (n/Z₀) dᵢ × Eᵢ projected on u, v"]
    EH --> PF(["<b>OpticsBase.PlaneField</b>"])
    PF -- "<b>forward</b>(f), n = 1" --> FF["<b>forward_field</b>(f)<br/>and <b>plane_frame</b>(f)"]
    FF --> WD(["<b>WavefrontBeamletDecomposition</b>(f)"])
    FF --> GD(["<b>GaussianModeDecomposition</b>(f)"])
    WD -- "position f.origin" --> BG["AstigmaticBeamGroup"]
    GD -- "starts at the centroid of the field" --> AB["AstigmaticGaussianBeamlet"]

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class DH,SM blue
    class OH,HFD,HA,HG,HF green
    class DA,PH,CW,EH,FF purple
    class DET,BEAM,PF,WD,GD,BG,AB terminal
```
