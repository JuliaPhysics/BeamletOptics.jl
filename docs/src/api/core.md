# Core design

The BMO package is intended to provide optical simulation capabilites with as much "out-of-the-box" comfort as possible, meaning that users should not have to worry about e.g. providing the correct optical sequence and exact alignment of objects. The following five design principles are core assumptions of the underlying API:

1. Optical interactions are decoupled from the underlying geometry representation
2. Optical elements are closed volumes or must mimic as such (exceptions apply, e.g. coatings)
3. Elements should be easily moveable and have working interactions for most angles of incidence
4. Without additional knowledge, tracing is performed non-sequentially
5. With additional knowledge, tracing is performed sequentially

## Intersect-Interact-Repeat-Loop

The first two principles will be elaborated upon in more detail in the [Geometry representation](@ref) section. For the latter two design decisions, the following high-level solver schematic can be used to abstract the steps that are performed when calling [`solve_system!`](@ref) with an input system and beam:

```mermaid
flowchart TB
    S(["<b>solve_system!(system, beam)</b>"])
    I["<b>1. Intersect</b><br/>find the closest surface<br/>that the ray hits, and all<br/>surfaces at that point"]
    Q1{"hit?"}
    F["<b>first_to_interact</b><br/>thin layer, then entered body,<br/>then left body"]
    X["<b>2. Interact</b><br/>the hit object computes<br/>the next ray"]
    T["<b>transition(system, ray)</b><br/>media on both sides,<br/>thin layer"]
    Q2{"next ray?"}
    R["<b>3. Repeat</b><br/>append the next ray or children,<br/>assign their enclosure,<br/>keep the hint of the interaction"]
    B(["back to<br/><b>1. Intersect</b>"])
    E(["done"])

    S --> I --> Q1
    Q1 -- yes --> F --> X
    X -- "body or thin layer" --> T --> Q2
    X -- "own method,<br/>e.g. mirror" --> Q2
    Q1 -- "no, or r_max rays" --> E
    Q2 -- yes --> R --> B
    Q2 -- "no, e.g. blocked" --> E

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class I,Q1,F,Q2,R purple
    class X green
    class T blue
    class S,E,B terminal
```

This scheme is loosely referred to as the **Intersect-Interact-Repeat-Loop** and consists of the following steps:

1. Calculate the closest **Intersection** between a ray/beam and the objects within the system. All surfaces at the closest point are kept as [Coincident hits](@ref), e.g. the faces of a cemented doublet, and the one that interacts comes first.
2. Calculate the optical **Interaction**: `interact3d` is dispatched on the object that was hit. An object with its own method computes the interaction itself, e.g. a mirror or a detector. The default methods of refracting bodies and thin layers read the [Transitions](@ref) at the hit, the media on both sides, and hand the ray to the surface model of the boundary or to the thin layer.
3. Attach the next part of the ray chain or beam tree. The solver assigns the bodies that enclose the new rays.
4. Use the new information to repeat 1.

Once this procedure has been completed, the alignment of the system or other time-dependent optical properties (e.g. the phase of a Gaussian beamlet) can be updated. A new call of [`solve_system!`](@ref) solves the beam again from its start. This is described in more detail in the section: [Tracing systems](@ref).

The next sections will focus on the **Intersection** and **Interaction** steps.

## Intersections

Calculating intersections between straight lines, i.e. rays, and surfaces is a central challenge for every geometrical optics simulator. This must be done with high numerical precision, since many optical effects are sensitive on the order of the wavelength of the light under consideration with respect to position and direction [Hecht:2018; p. 265 ff](@cite). In order to define this mathematically or algorithmically, many different methods exist [Hanrahan:1989](@cite). The first question is, how is the geometry of the problem defined. This topic is treated in the [Geometry representation](@ref) section. The second question concerns then the algorithm or equation that allows to calculate the point of intersection between a ray and the surface of the element. This function is called `intersect3d` and is, at its core, defined for each `shape` and `ray`:

```@docs; canonical=false
BeamletOptics.intersect3d(::BeamletOptics.AbstractShape, ::BeamletOptics.AbstractRay)
```

Regardless of the underlying concrete implementation, each call of `intersect3d` must return `nothing` or the following type:

```@docs; canonical=false
BeamletOptics.Intersection
```

Since an optical element can consist of multiple joint shapes, the return type must store which specific part of the object was hit: the `shape` of the `Intersection` is that part, and the solver compares and looks up the hit part by this reference, e.g. to decide whether the rays of a beamlet hit the same shape. The `object` of the `Intersection` is set by the solver, the `intersect3d` method of a shape leaves it `nothing`.

The closest hit of a ray is found on two levels: the solver searches the system for the objects that the ray can reach, and each object it reaches computes its own hit.

The search, shown below, goes through `children(system)`, the top-level objects and groups of the system, one at a time and keeps the closest hit found so far. Once there is a hit, every part whose bounding sphere the ray enters only beyond that hit (plus `get_coincidence_tolerance()`, so that touching surfaces are still found) is skipped, see [Bounding sphere pretest](@ref). A group is entered only if the ray reaches its sphere; its objects are then searched in the same way. `trace_all` starts without a hit. `trace_one` starts with the hit of the hinted object, so it can skip more parts from the start, see [Hints](@ref). A [`StaticSystem`](@ref) has no groups while tracing, since it flattens them when it is built.

```mermaid
flowchart TB
    EA(["<b>trace_all</b><br/>best hit so far: none"])
    EO(["<b>trace_one</b><br/>best hit so far: hinted object"])
    N{"another object or group<br/>of children(system)?"}
    H{"is it the<br/>hinted object?"}
    BH{"best hit<br/>so far?"}
    BY{"<b>isbeyond</b><br>does its sphere start<br/>behind the best hit?<br/><i></i>"}
    K{"object<br/>or group?"}
    GM{"does the ray miss<br/>the group's sphere?"}
    O["<b>intersect3d(object, ray)</b><br/>the object computes its hit,<br/>see the next diagram"]
    M["<b>merge_hits</b><br/>the closest hit becomes <br>the best hit"]
    R(["best hit with its coincident<br/>hits, or nothing"])
    B1(["tested already,<br/>next part"])
    B2(["cannot be closer,<br/>next part"])
    B3(["skip the whole group,<br/>next part"])
    B4(["search its objects the<br/>same way, then next part"])
    B5(["next part"])

    EA --> N
    EO --> N
    N -- yes --> H
    N -- no --> R
    H -- yes --> B1
    H -- no --> BH
    BH -- yes --> BY
    BH -- no --> K
    BY -- yes --> B2
    BY -- no --> K
    K -- "object" --> O --> M --> B5
    K -- "group (System only)" --> GM
    GM -- yes --> B3
    GM -- no --> B4

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class N,H,BH,K,M purple
    class BY,GM blue
    class O green
    class EA,EO,R,B1,B2,B3,B4,B5 terminal
```

Where surfaces touch, e.g. the faces of a cemented doublet, the ray hits several of them at the same point. `merge_hits` keeps the hits within `get_coincidence_tolerance()` of the closest one as a chain, see [Coincident hits](@ref).

Each object that the search reaches computes its hit with `intersect3d(object, ray)`, shown below. The default method tests the bounding sphere of the object first and then its shape or parts. A component can define its own method instead, as the plate beamsplitters do, e.g. the [`RectangularPlateBeamsplitter`](@ref). Once there is a hit, the search skips the object by this sphere, so the sphere must enclose every hit of such a method as well, see [Bounding sphere pretest](@ref).

```mermaid
flowchart TB
    IO(["<b>intersect3d(object, ray)</b>"])
    D{"does the object have<br/>its own intersect3d?"}
    OWN["<b>own method</b><br/>e.g. plate beamsplitters,<br/>hits must lie in its sphere"]
    SPH{"does the ray miss the<br/>object's bounding sphere?"}
    T{"one shape or<br/>several parts?<br/><i>shape_trait_of</i>"}
    SH["<b>intersect3d(shape, ray)</b><br/>SDF or mesh test,<br/>see below"]
    P["<b>intersect3d(parts, ray)</b><br/>each part as above,<br/>merged by <b>merge_hits</b>"]
    OBJ["<b>object!(hit, object)</b><br/>labels the hits<br/>with this object"]
    NO(["nothing"])
    HIT(["hit of the object with its<br/>coincident hits, or nothing"])

    IO --> D
    D -- yes --> OWN --> HIT
    D -- no --> SPH
    SPH -- yes --> NO
    SPH -- no --> T
    T -- "one shape" --> SH --> OBJ
    T -- "several parts" --> P --> OBJ
    OBJ --> HIT

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class D,T,P,OBJ purple
    class SPH,SH blue
    class OWN green
    class IO,NO,HIT terminal
```

The shape backend is selected by dispatch on the type of the shape, see [Geometry representation](@ref).

```mermaid
flowchart TB
    SH{"<b>intersect3d(shape, ray)</b>"}
    SDF["<b>AbstractSDF</b><br/>Ray marching "]
    MESH["<b>AbstractMesh</b><br/>Möller-Trumbore algorithm"]
    OWN["<b>own AbstractShape</b><br/>own intersect3d method"]
    IX(["<b>return ::Intersection or ::Nothing</b>"])

    SH -- isa AbstractSDF --> SDF
    SH -- isa AbstractMesh --> MESH
    SH -- "other shape" --> OWN
    SDF --> IX
    MESH --> IX
    OWN --> IX

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class SDF,MESH blue
    class OWN green
    class SH purple
    class IX terminal
```

### Coincident hits

Touching elements, e.g. the two lenses of a cemented doublet or a coating on its substrate, share a surface. A ray that reaches it hits the surfaces of both within the accuracy of the intersection algorithms. Instead of picking one of them, the solver keeps all hits within `get_coincidence_tolerance()` (1 nm by default) of the closest hit as a chain, linked by the `next` field of the [`BeamletOptics.Intersection`](@ref):

```@docs; canonical=false
BeamletOptics.coincident_hits
BeamletOptics.merge_hits
```

Hence the hit of the ray does not depend on the order of the objects in the system or on a [`BeamletOptics.Hint`](@ref). Each hit holds its own object, shape and normal. One of them interacts with the ray: a surface that is no refracting body comes first, e.g. a coating, a mirror or a detector, then the surface of the body that the ray enters, then that of the body it leaves:

```@docs; canonical=false
BeamletOptics.first_to_interact
BeamletOptics.interaction_rank
```

The `object` of every hit is the object of the system (or of an object group) that was hit. For an object of several parts, e.g. a [`DoubletLens`](@ref), the part follows from the shape:

```@docs; canonical=false
BeamletOptics.part_of(::BeamletOptics.AbstractObject, ::Any)
```

### Bounding sphere pretest

Since the tracing is non-sequential, every ray is tested against every shape of the system, and most of these tests are misses. A shape can therefore provide a sphere that encloses it. A ray that does not pass through this sphere is a miss without a call of `intersect3d`. This is optional: a shape without a sphere is tested as described above.

```@docs; canonical=false
BeamletOptics.bounding_sphere_of(::BeamletOptics.AbstractShape)
```

All [Signed Distance Functions (SDFs)](@ref) and [Meshes](@ref) have a bounding sphere. An object of several shapes and an object group have the sphere around the spheres of their parts, such that a ray that misses it is tested against none of them. The table holds spheres of floating-point shapes only: for other number types, e.g. dual numbers under `ForwardDiff`, the pretest is skipped and the shapes are tested directly.

The spheres are not stored in the shapes or objects. [`solve_system!`](@ref) computes them once per call, for the poses that the objects have at that time, via `bounding_sphere_of(x)` into a table, which the intersection code reads per ray via `bounding_sphere_of(table, x)`. The solve of a beam group computes the table once for all its beams, and the solves of its beams keep it instead of computing it again. Outside of a solve there is no table, i.e. a direct call of `intersect3d(object, ray)` tests the object without a sphere.

The table lives in the task-local storage of the solving task, not in a `ScopedValue`, because the intersection code reads it for every object and ray and this read is several times cheaper. A task that is spawned within a solve does not see the table of its parent. It passes the table on by calling [`BeamletOptics.with_bounding_spheres`](@ref)`(f, table)`, as the solve of a beam group does. The system must not change while a solve runs, since the spheres are those of the poses at the start. The types in the sphere code are asserted, and the methods for the parts of an object are kept separate from those for a shape, because inference would otherwise widen their argument types and the shapes would be tested by code compiled for abstract types.

```@docs; canonical=false
BeamletOptics.bounding_sphere_of(::BeamletOptics.BoundingSphereTable, ::Any)
BeamletOptics.SingleBoundingSphere
BeamletOptics.MultiBoundingSphere
BeamletOptics.NoBoundingSphere
```

The spheres of a system can be shown via `render!(ax, BeamletOptics.bounding_sphere_of(x))` for one sphere or the `show_bounding_sphere` keyword of [`render!`](@ref).

## Interactions

Optical interactions are performed after the point of intersection has been determined. The `interact3d` interface allows users to implement algorithms that calculate or try to mimic optical effects. The fidelity of the algorithm is effectively only limited by the amount of information that can be passed into the `interact3d` interface. The method is defined as follows:

```@docs; canonical=false
BeamletOptics.interact3d(::BeamletOptics.AbstractSystem, ::BeamletOptics.AbstractObject, ::BeamletOptics.AbstractBeam, ::Union{BeamletOptics.AbstractRay, Int})
```

As with the `intersect3d` method, a predefined return type must be provided in order to make the [`solve_system!`](@ref) interface work. The `AbstractInteraction` is used in order to create "building blocks" from which the output beam is constructed.

```@docs; canonical=false
BeamletOptics.AbstractInteraction
```

An object whose `interact3d` method stores data in the object, e.g. a detector its hits, also implements [`initialize!`](@ref), which discards this data before a new solve.

### Transitions

At a boundary, the ray passes from one medium into another. The parts that the ray hits there play two roles: a part with a medium is a **body**, which the ray leaves or enters, e.g. a lens or a prism; a part without one is a **thin layer** between the media, e.g. a thin beamsplitter or a polarizing film. The [`BeamletOptics.transition`](@ref) at the hit of a ray holds both media and the layer. Hence a component needs no knowledge of the elements that touch it: the surface of a lens cemented to another lens refracts from one glass into the other, a beamsplitter coating on a substrate transmits into the substrate, and a film between two plates transmits from one plate into the other.

The physics at a boundary has one place: the [`BeamletOptics.surface_model`](@ref) of the transition. It is the model that the thin layer at the boundary carries, else an uncoated surface. Bodies and thin layers both hand their rays to it, so a coating acts the same whether it stands free or lies on a substrate. The diagram shows this and how the solver assigns the enclosure of the new rays.

```mermaid
flowchart TB
    IA(["<b>interact3d</b>(system, object, beam, ray)"])
    MS{"<b>shape_trait_of</b>(object)"}
    PO["<b>part_of</b>(hit)<br/>the part that was hit,<br/>e.g. the back lens of a doublet"]
    K{"own interact3d<br/>of the part?"}
    HK{"fallback: <b>medium</b> or<br/><b>surface_model</b> declared?"}
    W(["warning, nothing:<br/>the object blocks the beam"])
    TR["<b>transition</b>(system, ray)<br/><b>boundary_hits</b>: left body, entered body,<br/>thin layer; media from the enclosure"]
    BOUND["<b>interact3d</b>(system, t, beam, at)<br/>beamlet: to the model of a thin layer,<br/>else per component ray"]
    SM["<b>surface_model</b>(t)<br/>model of the thin layer, else Uncoated:<br/>its channels, see Surface channels"]
    OWN["<b>own interact3d</b><br/>e.g. mirror, detector (a detector<br/>that does not stop: as a layer<br/>without a model)"]
    ENC["<b>enclose_new!</b><br/>enclosure of each new ray at the hit:<br/>back into the medium before the boundary<br/>or into the medium behind it"]
    E([next tracing step])

    IA --> MS
    MS -- MultiShape --> PO --> K
    MS -- SingleShape --> K
    K -- "no" --> HK
    HK -- "neither" --> W
    HK -- "body or<br/>thin layer" --> TR --> BOUND --> SM --> ENC
    K -- "own" --> OWN --> ENC
    ENC --> E

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class TR,BOUND,SM,ENC blue
    class OWN green
    class MS,K,HK purple
    class IA,E,W terminal
```

```@docs; canonical=false
BeamletOptics.AbstractMedium
BeamletOptics.IsotropicMedium
BeamletOptics.medium
BeamletOptics.transition
BeamletOptics.Transition
```

Every ray carries the bodies that enclose it, which the transitions update. Hence the medium behind a boundary is the body around the one that the ray leaves, e.g. the water of a cuvette around a lens:

```@docs; canonical=false
BeamletOptics.Enclosure
```

A component creates its new rays without an enclosure. After each interaction, the solver assigns it to every new ray that starts at the hit, the continued ray and the first rays of child beams, by the side of the boundary into which the ray propagates: back into the medium before it, e.g. after a reflection, or into the medium behind it.

```@docs; canonical=false
BeamletOptics.enclose!
```

The physics at the boundary belongs to the boundary, not to one of the bodies. A body, e.g. any [`BeamletOptics.AbstractRefractiveOptic`](@ref), and a thin layer, e.g. a [`ThinBeamsplitter`](@ref), hand the ray to the [`BeamletOptics.surface_model`](@ref) of the transition:

```@docs; canonical=false
BeamletOptics.interact3d(::BeamletOptics.AbstractSystem, ::BeamletOptics.Transition, ::BeamletOptics.AbstractBeam, ::Union{BeamletOptics.AbstractRay, Int})
BeamletOptics.surface_model(::BeamletOptics.Transition)
BeamletOptics.surface_model(::BeamletOptics.AbstractObject)
BeamletOptics.AbstractSurfaceModel
```

#### Surface channels

A surface model does not build rays. It describes the **channels** of the boundary: the directions into which the ray leaves it, reflected back into the medium before it or transmitted into the medium behind it, and how each transfers the field. The solver builds the new rays from them, once for every beam type: the continued ray, child beams, or child beamlets whose component rays leave their own hits on the same side as the chief ray. A model implements only `channels`, and optionally `splitting_of`; it needs no code per beam type. The channels of the bare boundary come first, so a coating changes the Fresnel channels instead of computing refraction again (green: implemented by a model author, purple: dispatch).

```mermaid
flowchart TB
    T(["<b>interact3d</b>(system, model, t, beam, ray)"])
    B["<b>channels</b>(t, ray)<br/>bare boundary: (transmitted, reflected),<br/>refraction, TIR, Fresnel"]
    M["<b>channels</b>(model, t, ray, cs)<br/>the model changes the channels,<br/>e.g. r and t of a coating"]
    S{"<b>splitting_of</b>(model)"}
    D["<b>Dominant</b><br/>the channel with the<br/>largest amplitude factor"]
    A["<b>AllChannels</b><br/>every channel that<br/>carries light"]
    R["<b>realize</b> / <b>beamlet_interaction</b><br/>one method per beam type"]
    BI(["BeamInteraction:<br/>one channel continues the beam"])
    CH(["children!, nothing:<br/>a child beam per channel"])
    NO(["nothing:<br/>no channel carries light"])

    T --> B --> M --> S
    S -- Dominant --> D --> R
    S -- AllChannels --> A --> R
    R --> BI
    R --> CH
    R --> NO

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class B,D,A,R blue
    class M green
    class S purple
    class T,BI,CH,NO terminal
```

```@docs; canonical=false
BeamletOptics.interact3d(::BeamletOptics.AbstractSystem, ::BeamletOptics.AbstractSurfaceModel, ::BeamletOptics.Transition, ::BeamletOptics.Beam, ::BeamletOptics.AbstractRay)
BeamletOptics.BoundaryChannel
BeamletOptics.AbstractSide
BeamletOptics.Transmitted
BeamletOptics.Reflected
BeamletOptics.channels
BeamletOptics.transfer
BeamletOptics.AbstractSplitting
BeamletOptics.splitting_of
BeamletOptics.Dominant
BeamletOptics.AllChannels
BeamletOptics.select_channels
BeamletOptics.realize
BeamletOptics.beamlet_interaction
BeamletOptics.power_fraction
```

The models of BMO:

```@docs; canonical=false
BeamletOptics.Uncoated
BeamletOptics.cross_boundary
BeamletOptics.BeamsplitterCoating
BeamletOptics.JonesFilm
```

A new refracting component declares its body via `medium`, a new thin layer its model via `surface_model(layer::MyLayer)`, like the [`ThinBeamsplitter`](@ref). Neither needs an `interact3d` method: without one, a body or a thin layer interacts at the boundary, and a beamlet at a thin layer is handed to its model, which may split it. An object that is neither, e.g. a mirror, implements `interact3d`. New physics at a boundary is a new `AbstractSurfaceModel` with a `channels` method, see [Surface channels](@ref).

The `interact3d` return type holds one new `beam` segment per interaction, which the solver appends to the current beam. Further beams, e.g. the reflected part at a beamsplitter, are attached to the beam as children by the component, and the solver queues them after the current beam, see [Tracing systems](@ref). After e.g. a lens surface air-to-glass interaction, the most likely next step is the refraction at the glass-to-air surface of the same lens. An interaction can tell the solver so with a `Hint`, which speeds up the next intersection.

### Re-emitting components

Some components do not continue a beam where it ended, but hand it to another solver and emit the result elsewhere, e.g. a fiber that takes the light at its input facet and releases it at its output facet. Their `interact3d` method samples the incoming beamlet segment as a field, lets the other solver compute the field at the exit and starts new beamlets there. The exchange format is the `PlaneField` of [OpticsBase.jl](https://github.com/StackEnjoyer/OpticsBase.jl) (see [OpticsBase interoperability](@ref)):

```julia
function BeamletOptics.interact3d(::AbstractSystem, fiber::MyFiber, agb::AstigmaticGaussianBeamlet, id::Int)
    f = PlaneField(agb, id; size, spacing, origin, axes)      # the segment that hit the fiber, on its input plane
    g = other_solver(f)                                       # PlaneField at the exit, with the full phase
    BeamletOptics.relaunch!(agb, [GaussianModeDecomposition(g)])
    return nothing
end
```

The new beamlets are attached with [`BeamletOptics.relaunch!`](@ref) instead of `children!`. Children added by `children!` continue the parent: its length, optical path and reference plane are counted on. A re-emitted beamlet already carries the full phase in its amplitude and starts at another place, so it must count its optical path from its own start.

```@docs; canonical=false
BeamletOptics.relaunch!
```

`PlaneField(beamlet, id; ...)` is for components without a detector, e.g. in the `interact3d` of a component that hands the incoming beamlet to another solver: the segment that just hit the component is `id`. It builds a single hit of the segment and passes it to `PlaneField(hits; ...)`, which also takes the hits of several segments.

```@docs
OpticsBase.PlaneField(::AstigmaticGaussianBeamlet, ::Integer)
OpticsBase.PlaneField(::AbstractVector{<:Union{BeamletOptics.GaussianBeamletHit, BeamletOptics.AstigmaticGaussianBeamletHit}})
```

## Hints

The [`BeamletOptics.Hint`](@ref) interface allows developers to tell the non-sequential solver which object will probably be hit during the next cycle of the [Intersect-Interact-Repeat-Loop](@ref). The solver tests this object first. Its hit lets the solver skip all objects whose bounding sphere lies beyond it, which saves most of the intersection tests in a large system.

```@docs; canonical=false
BeamletOptics.Hint
```

A `Hint` never changes the result of the trace: the hinted object is not preferred over a closer one, and an object that is moved away from its hinted neighbour, e.g. the back lens of a doublet, leaves an air gap that the trace sees. Touching elements need no `Hint`, see [Coincident hits](@ref) and [Transitions](@ref). The diagram shows how the next tracing step uses the hint of an interaction.

```mermaid
flowchart TB
    S(["<b>push!(beam, interaction)</b><br/>next ray, hint(interaction)"])
    HQ{"did the interaction<br/>give a hint?<br/><i>tracing_step!</i>"}
    HS["<b>trace_one(system, ray, hint)</b><br/>intersect3d(object(hint), ray) first"]
    OB["<b>closest_hits</b><br/>of all other objects, those<br/>beyond the hint's hit skipped"]
    TA["<b>trace_all(system, ray)</b><br/>closest_hits of all objects"]
    FI["<b>first_to_interact</b><br/>the hit that interacts first"]
    E(["on to <b>interact3d</b>"])

    S --> HQ
    HQ -- yes --> HS --> OB --> FI
    HQ -- no --> TA --> FI
    FI --> E

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class HQ,HS,OB,TA,FI purple
    class S,E terminal
```

The solver picks the `interact3d` method by the *object* of the hint, not by its shape. For a [`BeamletOptics.MultiShape`](@ref) component the hint therefore names the composite as `object` and the part that is hit next as `shape`, e.g. `Hint(cbs, shape(cbs.coating))`; `Hint(object)` alone fails for such an object, since `shape(object)` is a tuple.

!!! tip
    Any component can use a `Hint`: a single lens hints itself when a ray enters it or is totally internally reflected, since the next intersection lies on the same lens.

## Tracing logic

### Tracing systems

In the initial state, it is assumed that the problem consists of `objects <: AbstractObject` (in a system) and a `beam <: AbstractBeam` with a defined starting position and direction. No additional information is provided, and the specific path of the beam is not known beforehand. Consequently, brute force tracing of the optical system is required, involving testing against each individual element to determine the trajectory of the beam.

A call of [`solve_system!`](@ref) traces the beam tree breadth-first: it traces one beam after another and queues the children that the components attach, up to the keyword `depth_max` of [`solve_system!`](@ref) branching levels. Children that carry less than the fraction `power_cutoff` of the power of the root beam are dropped before they are queued.

```mermaid
flowchart TB
    S(["<b>solve_system!</b>"])
    INI["<b>initialize!(system)</b><br/>initialize! each contained object"]
    BS["<b>with_bounding_spheres</b><br/>precalculate spheres of all objects,<br/>unless a beam group did"]
    Q["<b>popfirst!(queue)</b><br/>next (beam, depth),<br/>the root at depth 1"]
    EM["<b>empty!(current)</b><br/>back to its first ray"]
    TS{"<b>trace_system!</b><br/>(system, current)"}
    D{"depth #lt;= depth_max?"}
    PC["<b>drop_weak_children!</b><br/>children below power_cutoff<br/>are discarded"]
    ENQ["<b>push!(queue, ...)</b><br/>children(current) at depth + 1,<br/>attached via children! or relaunch!"]
    DROP["<b>_drop_beams!(current)</b><br/>children are discarded"]
    QE{"queue empty?"}
    BQ(["back to<br/><b>popfirst!(queue)</b>"])
    E([done])

    S -. "initialize = true" .-> INI
    S --> BS
    INI --> BS
    BS --> Q --> EM --> TS
    TS -- "method per beam type" --> D
    D -- yes --> PC --> ENQ --> QE
    D -- no --> DROP --> QE
    QE -- no --> BQ
    QE -- yes --> E

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class BS,Q,EM,ENQ,DROP blue
    class INI green
    class TS purple
    class S,E,BQ terminal
```

```@docs; canonical=false
BeamletOptics.trace_system!
```

For a [`Beam`](@ref), `trace_system!` runs the loop below until the beam has `r_max` rays, its last ray hits nothing or `interact3d` returns `nothing`; the methods for beamlets run one `tracing_step!` per component ray and one `interact3d` for the whole beamlet per step, and also stop if the rays hit different shapes or the optical invariant check fails.

```mermaid
flowchart TB
    S(["<b>trace_system!(system, beam)</b><br/>with a ::Beam"])
    RM{"fewer than<br/>r_max rays?"}
    HQ{"did the last interaction<br/>give a hint?<br/><i>tracing_step!</i>"}
    TO["<b>trace_one(system, ray, hint)</b><br/>closest hit, hinted object first"]
    TA["<b>trace_all(system, ray)</b><br/>closest hit, no object first"]
    HIT{"did the ray<br/>hit anything?"}
    FI["<b>first_to_interact</b><br/>of the hits at that point"]
    IA["<b>interact3d(system, object, beam, ray)</b><br/>the hit object computes the next ray"]
    NI{"did it return<br/>an interaction?"}
    EN["<b>enclose_new!</b><br/>enclosure of the new rays:<br/>next ray or children"]
    P["<b>push!(beam, interaction)</b><br/>append the next ray"]
    E(["done"])
    BR(["back to<br/><b>fewer than r_max rays?</b>"])

    S --> RM
    RM -- yes --> HQ
    RM -- no --> E
    HQ -- yes --> TO --> HIT
    HQ -- no --> TA --> HIT
    HIT -- yes --> FI --> IA --> NI
    HIT -- no --> E
    NI -- yes --> EN --> P --> BR
    NI -- "no, e.g. blocked" --> E

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class RM,HQ,TO,TA,HIT,FI,NI,EN,P purple
    class IA green
    class S,E,BR terminal
```

This non-sequential mode is comparatively safe in determining the "true" beam path, but will scale suboptimally in time-complexity with the amount of optical elements.

!!! info "Object order"
    Unlike with classic, surface-based ray tracers, the order in which objects are listed in the [`System`](@ref) object vector/tuple is not considered for the purpose of tracing.

## CPU and GPU support

Parallelizing the execution of a [`solve_system!`](@ref) call on the CPU is straight-forward, since the solution is independent for each beam or ray and the solver can run on multiple threads: the beams of a beam group are solved in parallel, see below. Objects that are mutated during a solve, e.g. detectors like the [`Detector`](@ref), need special consideration when implementing them, since multiple threads might be able to access the underlying memory, leading to race conditions. Specifically, this means ensuring e.g. atomic write and read access. The [`Detector`](@ref) takes a lock for its hits and, within a beam group solve, collects them in a buffer per chunk of beams.

With respect to GPU acceleration, this is not the case. Currently, all available implementations of [`solve_system!`](@ref) are highly branching algorithms which can not be implemented in a parallelized way easily. This will most likely require a specific new subtype of the [`BeamletOptics.AbstractSystem`](@ref) with determinable sequential properties. This is not a development goal as of the writing of this section.

The beams of a group are solved in chunks, and all beams of a group share one table of bounding spheres. It belongs to a task, so each task of the threads sets it for its beams. Each chunk collects its detector hits and stores them at its end. A free task takes the next chunk, because the beams of a group differ in cost (e.g. the inner rings of a source pass a lens, the outer ones miss it), and those of similar cost lie next to each other. Four chunks per thread balance most of the load. Each further chunk costs a channel access and a flush, which shows if all beams cost the same.



