# Geometry representation

In BMO, the distinction between an *object* and its geometric representation (*shape*) is a central design principle. This separation is intended to ensure flexibility and modularity in modeling optical components. Unlike the surface based representation of optical elements in other tools, they are **treated as volumetric bodies** in this simulation framework.

## Separation of geometry and optical interactions

An object does not handle its geometry itself. Everything geometric, e.g. `intersect3d`, `position` or `orientation`, is forwarded to its shape trait, which returns the shape or shapes of the object. Moving an object works the same way: objects are [`BeamletOptics.Movable`](@ref) by default, and the movable methods of an object forward `translate3d!` and `rotate3d!` to the shape trait, which then moves every shape of the object. Solid arrows are labeled with the function or field that leads from one type to the next, dotted arrows point to subtypes. What a subtype has to implement is listed in the docstring of each type.

```mermaid
flowchart TB
    OBJ["<b>AbstractObject</b><br/><i>abstract type</i>"]
    MV["<b>Movable</b><br/><i>kinematic trait</i>"]
    TRT["<b>AbstractShapeTrait</b><br/><i>shape trait</i>"]
    SGL["<b>SingleShape</b><br/><i>shape trait</i>"]
    MLT["<b>MultiShape</b><br/><i>shape trait</i>"]
    SHP["<b>AbstractShape</b><br/><i>abstract type</i>"]

    OBJ -- kinematic_trait_of --> MV
    MV -- "translate3d!<br/>rotate3d!" --> TRT
    OBJ -- shape_trait_of --> TRT
    TRT -.-> SGL
    TRT -.-> MLT
    SGL -- object.shape --> SHP
    MLT -- "shape(object)" --> SHP

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    class SHP blue
    class OBJ green
    class MV,TRT,SGL,MLT purple
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

An `AbstractObject` can consist of multiple `AbstractShape`s or even multiple subsidiary `AbstractObject`s, facilitating the creation of composite optical elements. For example, a lens with an anti-reflective coating could be represented as the substrate and a seperate model for the coating, each with its own geometric and optical properties. In general, an `object` can have the `SingleShape` or a `MultiShape` trait. The [`BeamletOptics.AbstractShapeTrait`](@ref) allows the solver and kinematic API to apply different implementations of basis functions via multiple dispatch. 

```@docs; canonical=false
BeamletOptics.SingleShape
BeamletOptics.MultiShape
```
