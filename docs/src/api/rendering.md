# Rendering internals

This page explains how [Live rendering](@ref) is implemented, and what a developer has to provide so that a new component, shape or handle type can be drawn and kept up to date. The user-facing functions are described in the Basics section.

## Live rendering dispatch

The diagram shows which `live_render!` method draws a thing and which handle it returns. BMO provides `render!(ax, obj)` and `render!(ax, shape)` as defaults, which a component or shape replaces with an own method to be drawn differently:

```mermaid
flowchart TB
    S1["<b>live_render!(ax, system)</b><br/>each top-level object,<br/>groups per object"]
    S2["<b>live_render!(h, obj)</b><br/>each object of obj, added to h;<br/>ArgumentError if rendered,<br/>plots deleted if it fails"]
    B1["<b>live_render!(ax, beam)</b><br/>ray, beam, beam group or beamlet:<br/>plots of render!, live = true,<br/>coarser defaults for beamlets"]
    O1["<b>live_render!(ax, obj)</b><br/>object or object group"]
    SH["<b>SystemRenderHandle</b><br/><i>#lt;: AbstractSystemRenderHandle</i>"]
    BH["<b>BeamRenderHandle</b><br/><i>#lt;: AbstractBeamRenderHandle</i>"]
    D["<b>live_render!(draw, ax, x)</b><br/>keeps the plots that draw() adds"]
    R["<b>render!(ax, obj)</b><br/>own method, e.g. PolarizationFilter<br/>object group: each object"]
    OH["<b>ObjectRenderHandle</b><br/><i>#lt;: AbstractObjectRenderHandle</i>"]
    TR["<b>render!(ax, shape_trait_of(obj), obj)</b><br/>MultiShape: render! of each part"]
    SHP["<b>render!(ax, shape)</b><br/>or the analytic mesh of the shape"]

    S1 -- "each object" --> O1
    S2 -- "each object" --> O1
    S1 --> SH
    S2 -. "append! to h.handles" .-> SH
    B1 --> BH
    O1 -- "draw = () -> render!(ax, obj)" --> D
    D -- "draw()" --> R
    D --> OH
    SH -. "children(h)" .-> OH
    R -- "AbstractObject default" --> TR
    TR -- "SingleShape: shape(obj)" --> SHP

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    class S1,S2,O1,D,B1,TR blue
    class R,SHP green
    class SH,OH,BH purple
```

## Updating a handle

The diagram shows what `update_render!` and `rerender!` do for each handle type. `rerender!` skips
the checks and always draws an object again, e.g. after a change of its shape:

```mermaid
flowchart TB
    U(["<b>update_render!(h)</b>"])
    RE(["<b>rerender!(h)</b>"])
    RG["<b>_is_rigid(h)</b><br/>parts unchanged and<br/>moved rigidly?"]
    OB["<b>refill the observables</b><br/>obs[] = geometry()"]
    PC["<b>pose of the object</b><br/>unchanged?"]
    MM["<b>translate!</b>, <b>rotate!</b><br/>model matrix of the plots"]
    RR["<b>draw again</b><br/>delete the plots, draw() at the<br/>current pose, which becomes P0, R0"]
    E(["<b>return h</b>"])

    U -- "ObjectRenderHandle, or<br/>update_render!(oh) for each oh in<br/>children(h) of a system handle" --> RG
    U -- "BeamRenderHandle" --> OB
    RE -- "ObjectRenderHandle, or<br/>rerender!(oh) for each oh in<br/>children(h) of a system handle" --> RR
    RE -- "BeamRenderHandle" --> OB
    RG -- yes --> PC
    RG -- no --> RR
    PC -- yes --> E
    PC -- no --> MM

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class MM,RR,OB blue
    class RG,PC purple
    class U,RE,E terminal
```

## Render handle protocol

Every handle returned by `live_render!` is a [`BeamletOptics.AbstractRenderHandle`](@ref). Packages
built on BeamletOptics, e.g. [BeamletOpticsGUI](https://github.com/StackEnjoyer/BeamletOpticsGUI.jl),
use the handles through a protocol instead of the concrete handle types of the Makie extension. There
are three abstract handle types: [`BeamletOptics.AbstractObjectRenderHandle`](@ref) for a component,
[`BeamletOptics.AbstractSystemRenderHandle`](@ref) for a system and
[`BeamletOptics.AbstractBeamRenderHandle`](@ref) for a ray, beam or beam group. They share the
accessors below.

```@docs; canonical=false
BeamletOptics.rendered
BeamletOptics.render_plots
```

A system handle groups the handles of its objects. It can be searched and changed at runtime:

```@docs; canonical=false
BeamletOptics.children(::BeamletOptics.AbstractSystemRenderHandle)
BeamletOptics.render_parent
Base.push!(::BeamletOptics.AbstractSystemRenderHandle, ::BeamletOptics.AbstractObjectRenderHandle)
Base.delete!(::BeamletOptics.AbstractSystemRenderHandle, ::BeamletOptics.AbstractObjectRenderHandle)
```

[`update_render!`](@ref), [`rerender!`](@ref), [`remove_render!`](@ref) and [`pick_object`](@ref) work on any
system handle through these accessors, so a handle type of your own only needs to implement them.

The diagram shows which accessors and hooks the methods of a system handle call, and where they pass on to the object handles:

```mermaid
flowchart TB
    UR["<b>update_render!(h)</b><br/><b>rerender!(h)</b><br/><b>render_plots(h)</b>"]
    REM0["<b>remove_render!(h)</b>"]
    REM["<b>remove_render!(h, obj)</b>"]
    PCK["<b>pick_object(h, plot)</b>"]
    CH["<b>AbstractTrees.children(h)</b>"]
    PAR["<b>render_parent(h, obj)</b>"]
    DEL["<b>delete!(h, oh)</b><br/>after remove_render!(oh)"]
    OHM["<b>each object handle oh</b><br/>rendered, update_render!,<br/>rerender!, render_plots,<br/>remove_render!, pick_object"]
    PP["<b>pickable_plots(x, plots)</b><br/><i>method per type of x</i>"]

    UR --> CH
    REM0 --> CH
    REM --> CH
    REM -- "obj in a group:<br/>ArgumentError" --> PAR
    REM --> DEL
    PCK --> CH
    PCK -- "up to the<br/>outermost group" --> PAR
    CH -- "each oh; for remove_render!(h, obj)<br/>only if rendered(oh) in obj" --> OHM
    OHM -- "pick_object(oh, plot)" --> PP

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    class UR,REM0,REM,PCK blue
    class CH,PAR,DEL,PP green
    class OHM purple
```

[`pick_object`](@ref) uses `pickable_plots`; the colors of the material classes of the active look
are available via `look_colors`:

```@docs; canonical=false
BeamletOptics.pickable_plots
BeamletOptics.look_colors
```
