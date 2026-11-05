# Detector internals

This page explains how the [`Detector`](@ref) described in the basics is implemented internally, and what a developer has to provide so that a new hit type or sampling works with the readouts. What a readout computes is defined in [Reading out a detector](@ref) and [Field distributions](@ref).

## Hit recording

The following diagram shows how a solve records one hit type per beam model in a `Detector` and what empties it again. A new [`BeamletOptics.AbstractDetector`](@ref) implements `interact3d` and `empty!`.

```mermaid
flowchart TB
    S(["<b>solve_system!</b><br/>"])
    IA["<b>interact3d</b><br/>One method per beam model"]
    HT["<b>AbstractDetectorHit</b> <br/> One type per beam model, <br/>e.g. <b>Ray</b> generates a <b>RayHit</b>"]
    P["<b>push!(d, hit)</b>"]
    IK(["<b>initialize!(d)</b> or<br/><b>translate3d!</b> or<br/><b>rotate3d!</b><br/>Also within an ObjectGroup"])
    B["<b>with_hit_buffers</b><br/>one buffer per chunk<br/>of beams and detector"]
    E["<b>empty!(d)</b>"]
    H["<b>hits(d)</b><br/>Returns <b>nothing</b> or </br>a <b>Vector</b> of all hits"]

    S -- beam hits detector --> IA -- generates --> HT --> P
    P -- "beam group" --> B
    P -- "otherwise,<br/>under d.lock" --> H
    B -- "append!(d, buffer)<br/>at the end of the chunk" --> H
    IK --> E
    E -- "resets: hits = nothing" --> H

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class IA,E green
    class HT purple
    class P,B,H blue
    class S,IK terminal
```

## Readout pipeline

The following diagram shows how a readout turns a sampling into points and evaluates the stored hits there. A new sampling type implements `sample_points`; a new hit type implements `electric_field(hit, p)`, `direction(hit)`, `refractive_index(hit)` and an `intensity` method for a vector of its hits. Only the [`RegularGrid`](@ref) calls `calc_local_lims`, to find its window around the hits; a new hit type therefore implements it for a vector of its hits with a `crop_factor` and for `crop_factor = nothing`. A new hit type must also be added to the `hits` field of [`Detector`](@ref).

```mermaid
flowchart TB
    R(["<b>electric_field</b> / <b>intensity</b> / <b>optical_power</b><br/>(pd, sampling = RegularGrid())"])
    SP["<b>sample_points(pd, sampling)</b><br/><i>dispatch on the sampling</i>"]
    BMO["<b>BMO stock samplings</b><br/>"]
    OWN["own <b>sample_points</b> method"]
    LIM["own <b>calc_local_lims(pd, hits, crop_factor)</b><br/>and the method for crop_factor = nothing,<br/>unless an AbstractRayHit<br/>or AbstractBeamletHit"]
    DS["<b>DetectorSamples</b><br/>coordinates, points, weights"]
    PT["<b>electric_field(pd, points)</b><br/><b>intensity(pd, points)</b><br/>threaded over the points,<br/>no hits: zero at every point"]
    SUM["<b>electric_field(hits(pd), p)</b><br/>coherent sum over the hits"]
    IH["<b>intensity(hits(pd), p, n̂)</b><br/>Poynting flux through the detector<br/>vector form ½ Re[(E × H*) · n̂]:<br/>PolarizedRayHit, AstigmaticGaussianBeamletHit<br/>scalar form ½ Re[E* F]:<br/>RayHit, GaussianBeamletHit"]
    RET(["electric_field, intensity:<br/>(coordinates..., values)<br/>optical_power:<br/>sum of weights .* I"])
    OWNI["own <b>intensity(hits, p, n)</b> method"]
    EF["<b>electric_field(hit, p)</b><br/>direction(hit), refractive_index(hit)"]

    R --> SP
    SP -- "RegularGrid,<br/>RectilinearGrid" --> BMO
    SP -- "own AbstractDetectorSampling" --> OWN
    BMO -- "RegularGrid only,<br/>own hit type" --> LIM
    BMO --> DS
    OWN --> DS
    DS -- points --> PT
    PT -- electric_field --> SUM
    PT -- intensity --> IH
    PT -- values --> RET
    IH -- "own hit type" --> OWNI
    SUM --> EF
    IH --> EF

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class SP,BMO,DS purple
    class PT,SUM,IH blue
    class OWN,LIM,OWNI,EF green
    class R,RET terminal
```

```@docs; canonical=false
BeamletOptics.AbstractDetectorSampling
BeamletOptics.sample_points
BeamletOptics.DetectorSamples
```

## OpticsBase extension internals

The diagram shows how the constructors share one sampling of the hits and how the return path leaves the `PlaneField` again.

```mermaid
flowchart TB
    DET(["<b>PlaneField(detector; size, spacing)</b>"])
    BEAM(["<b>PlaneField(beamlet, id; ...)</b><br/>size, spacing, origin, axes"])
    DH["<b>hits(detector)</b><br/>ArgumentError if the<br/>detector has no hits"]
    DA["<b>default_plane_axes(detector)</b><br/>u from local_axes,<br/>n = −direction"]
    OH["one <b>GaussianBeamletHit</b> or<br/><b>AstigmaticGaussianBeamletHit</b><br/>(beamlet, id)"]
    PH["<b>PlaneField(hits; ...)</b><br/>size, spacing, origin, axes<br/>beamlet hits only,<br/>otherwise ArgumentError"]
    CW["<b>common_wavelength_index(hits)</b><br/>all hits share λ and n"]
    HFD["<b>hit_field_direction(hit, u, v)</b><br/>AstigmaticGaussianBeamletHit: nothing,<br/>electric_field(hit, p) is a 3D vector<br/>GaussianBeamletHit: u made<br/>transverse to d (v if u ∥ d)"]
    SM["<b>_sample(hits, origin, u, v, points)</b><br/>threaded over points,<br/>shared with the readouts"]
    HF["<b>hit_field(hit, p, pol)</b><br/>per hit type"]
    EH["<b>sum over the hits</b><br/>E = Σᵢ Eᵢ projected on u, v<br/>H = Σᵢ (n/Z₀) dᵢ × Eᵢ projected on u, v"]
    PF(["<b>OpticsBase.PlaneField</b>"])
    FF["<b>forward_field(f)</b><br/>and <b>plane_frame(f)</b>"]
    WD(["<b>WavefrontBeamletDecomposition(f)</b><br/>AstigmaticBeamGroup at f.origin"])
    GD(["<b>GaussianModeDecomposition(f)</b><br/>AstigmaticGaussianBeamlet<br/>at the centroid of the field"])

    DET --> DH
    DET -- "default axes" --> DA
    BEAM --> OH
    DH --> PH
    DA -- axes --> PH
    OH --> PH
    PH --> CW --> HFD --> SM
    SM -- "per point p" --> HF --> EH --> PF
    PF -- "forward(f), n = 1" --> FF
    FF --> WD
    FF --> GD

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class DH,SM blue
    class OH,HFD,HF green
    class DA,PH,CW,EH,FF purple
    class DET,BEAM,PF,WD,GD terminal
```
