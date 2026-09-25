# Geometry representation

In BMO, the distinction between an *object* and its geometric representation (*shape*) is a central design principle. This separation is intended to ensure flexibility and modularity in modeling optical components. Unlike the surface based representation of optical elements in other tools, they are **treated as volumetric bodies** in this simulation framework.

## Separation of geometry and optical interactions

Each card lists the fields and functions that a subtype must provide. Solid arrows point from a field to the type it stores, dotted arrows point to subtypes.

```mermaid
flowchart TB
    SYS["<b>AbstractSystem</b><br/><i>abstract type</i><hr/>objects, n<br/>refractive_index(system, λ)"]
    GRP["<b>AbstractObjectGroup</b><br/><i>abstract type #lt;: AbstractObject</i><hr/>objects"]
    OBJ["<b>AbstractObject</b><br/><i>abstract type</i><hr/>interact3d(system, object, beam, ray)<br/>shape_trait_of(object)"]
    TRT["<b>AbstractShapeTrait</b><br/><i>trait</i>"]
    SGL["<b>SingleShape</b><br/><i>trait</i><hr/>object.shape"]
    MLT["<b>MultiShape</b><br/><i>trait</i><hr/>shape(object)"]
    SHP["<b>AbstractShape</b><br/><i>abstract type</i><hr/>pos, dir<br/>intersect3d(shape, ray)"]

    SYS -- objects --> OBJ
    GRP -- objects --> OBJ
    OBJ -- shape_trait_of --> TRT
    TRT -.-> SGL
    TRT -.-> MLT
    OBJ -- geometry --> SHP

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    class SHP blue
    class SYS,OBJ,GRP green
    class TRT,SGL,MLT purple
```

Light sources are structured analogously, but carry no shape:

```mermaid
flowchart TB
    BG["<b>AbstractBeamGroup</b><br/><i>abstract type</i><hr/>beams, center, orientation<br/>wavelength(group)"]
    BM["<b>AbstractBeam</b><br/><i>abstract type</i><hr/>parent, children<br/>first_ray(beam)<br/>empty!(beam)<br/>_modify_beam_head!<br/>_last_beam_intersection"]
    RY["<b>AbstractRay</b><br/><i>abstract type</i><hr/>pos, dir, intersection, λ, n<br/>empty!(ray)"]

    BG -- beams --> BM
    BM -- rays --> RY

    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    class BG,BM,RY green
```

The geometry, represented by a concrete subtype of the [`BeamletOptics.AbstractShape`](@ref), defines the physical boundaries of the element. Shapes can be represented in various forms, such as [Meshes](@ref) or [Signed Distance Functions (SDFs)](@ref). The main goal for this design choice is to allow for the possibility to switch out geometry representations for more advanced methods in the future, e.g. [NURBS.jl](https://github.com/HoBeZwe/NURBS.jl).

```@docs; canonical=false
BeamletOptics.AbstractShape
```

On the other hand, the optical behavior — how light interacts with the element — is defined by the [`BeamletOptics.AbstractObject`](@ref) type. This decoupling allows for independent development and extension of geometry representations and optical interaction models within the [Intersect-Interact-Repeat-Loop](@ref).

```@docs; canonical=false
BeamletOptics.AbstractObject
```

## Single and multi-shaped objects

An `AbstractObject` can consist of multiple `AbstractShape`s or even multiple subsidiary `AbstractObject`s, facilitating the creation of composite optical elements. For example, a lens with an anti-reflective coating could be represented as the substrate and a seperate model for the coating, each with its own geometric and optical properties. In general, an `object` can have the `SingleShape` or a `MultiShape` trait. The [`BeamletOptics.AbstractShapeTrait`] allows the solver and kinematic API to apply different implementations of basis functions via multiple dispatch. 

```@docs; canonical=false
BeamletOptics.SingleShape
BeamletOptics.MultiShape
```
