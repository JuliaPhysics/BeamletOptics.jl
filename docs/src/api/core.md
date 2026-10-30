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
    S(["<b>solve_system!</b><br/>system + beam"]) --> I
    I["<b>1. Intersect</b><br/>closest hit of all objects,<br/>hinted object first;<br/>all surfaces at that point"] -- "coincident hits" --> F["<b>first_to_interact</b><br/>thin layer, then entered body,<br/>then left body"]
    F --> X["<b>2. Interact</b><br/>interact3d(system, object, beam, ray),<br/>dispatch on the object"]
    X -- "body or thin layer" --> T["<b>transition(system, ray)</b><br/>media on both sides, thin layer"]
    X -- "own method,<br/>e.g. mirror" --> R
    T --> R["<b>3. Repeat</b><br/>attach next ray or children,<br/>assign their enclosure"]
    R -- "optional Hint" --> I
    I -- "no hit or r_max" --> E([done])

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class I,F blue
    class X,T green
    class R purple
    class S,E terminal
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

Since an optical element can consist of multiple joint shapes, the return type must store which specific part of the object was hit.

The following diagram shows how `trace_all` reaches the shapes of the objects and groups of a system, including the [Bounding sphere pretest](@ref) (green: implemented by a component author, purple: dispatch). Where surfaces touch, e.g. the faces of a cemented doublet, the ray hits several of them at the same point. The hits within `get_coincidence_tolerance()` of the closest one are kept as a chain, see [Coincident hits](@ref).

```mermaid
flowchart TB
    TA(["<b>trace_all(system, ray)</b><br/><b>closest_hits</b> of the top-level<br/>objects and groups"])
    BSO["<b>bounding_sphere_of(shape)</b><br/><i>optional</i>"]
    TBL["<b>bounding_spheres(system)</b><br/>table, once per solve"]
    P{"<b>intersect3d(parts, ray)</b><br/>each part, <b>merge_hits</b>:<br/>closest hit and those<br/>at the same point"}
    BY{"<b>beyond</b>(sphere, ray, t)?<br/>sphere only behind<br/>the closest hit so far"}
    G["<b>intersect3d(group, ray)</b><br/>x = objects(group)"]
    OJ["<b>intersect3d(object, ray)</b><br/><i>own method possible</i>"]
    T{"<b>shape_trait_of(object)</b>"}
    SS["<b>SingleShape</b><br/>x = shape(object)"]
    MS["<b>MultiShape</b><br/>x = shape(object), a tuple"]
    PT{"<b>intersect3d(sphere, x, ray)</b><br/>sphere of x from the table,<br/>else NoBoundingSphere"}
    SH["<b>intersect3d(shape, ray)</b><br/>shape backend, see below"]
    OBJ["<b>object!(intersection, object)</b><br/>by the SingleShape and MultiShape<br/>methods, not by a group"]
    MISS(["nothing"])

    BSO --> TBL
    TBL -. "current_bounding_spheres()" .-> PT
    TA -- "each top-level part" --> BY
    BY -- "yes: skipped" --> MISS
    BY -- no --> P
    P -- "AbstractObjectGroup" --> G
    P -- "AbstractObject" --> OJ
    P -- "AbstractShape" --> SH
    OJ --> T
    T -- SingleShape --> SS --> PT
    T -- MultiShape --> MS --> PT
    G --> PT
    PT -- "ray misses the sphere" --> MISS
    PT -- "x is a shape" --> SH
    PT -- "x are parts" --> P
    SH --> OBJ

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class TBL,G,SS,MS,OBJ blue
    class BSO,OJ,SH green
    class P,T,PT,BY purple
    class TA,MISS terminal
```

The shape backend is selected by dispatch on the type of the shape, see [Geometry representation](@ref).

```mermaid
flowchart TB
    SH{"<b>intersect3d(shape, ray)</b>"}
    SDF["<b>AbstractSDF</b><br/>ray marching from outside<br/>or inside, up to the sphere exit<br/>if tested with its own sphere"]
    MESH["<b>AbstractMesh</b><br/>MoellerTrumboreAlgorithm<br/>for each face, closest"]
    OWN["<b>own AbstractShape</b><br/>own intersect3d method"]
    IX(["<b>Intersection(t, n, shape)</b><br/>or nothing"])

    SH -- AbstractSDF --> SDF
    SH -- AbstractMesh --> MESH
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

All [Signed Distance Functions (SDFs)](@ref) and [Meshes](@ref) have a bounding sphere. An object of several shapes and an object group have the sphere around the spheres of their parts, such that a ray that misses it is tested against none of them.

The spheres are not stored in the shapes or objects. [`solve_system!`](@ref) computes them once per call, for the poses that the objects have at that time, via `bounding_sphere_of(x)` into a table, which the intersection code reads per ray via `bounding_sphere_of(table, x)`. Outside of a solve there is no table, i.e. a direct call of `intersect3d(object, ray)` tests the object without a sphere.

The table lives in the task-local storage of the solving task, not in a `ScopedValue`, because the intersection code reads it for every object and ray and this read is several times cheaper. A task that is spawned within a solve does not see the table of its parent. It passes the table on by calling `with_bounding_spheres(g, table)`, as the solve of a beam group does. The system must not change while a solve runs, since the spheres are those of the poses at the start. The types in the sphere code are asserted, and the methods for the parts of an object are kept separate from those for a shape, because inference would otherwise widen their argument types and the shapes would be tested by code compiled for abstract types.

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

The diagram shows how a body and a thin layer interact at a boundary and how the solver assigns the enclosure of the new rays (green: implemented by a component author, purple: dispatch).

```mermaid
flowchart TB
    IA(["<b>interact3d</b>(system, object, beam, ray)"])
    MS{"<b>shape_trait_of</b>(object)"}
    PO["<b>part_of</b>(hit)<br/>the part that was hit,<br/>e.g. the back lens of a doublet"]
    K{"own method<br/>of the part?"}
    TR["<b>transition</b>(system, ray)<br/><b>boundary_hits</b>: left body, entered body,<br/>thin layer; media from the enclosure"]
    BODY["<b>AbstractRefractiveOptic</b><br/><b>interact3d</b>(system, surface_model(t), t, beam, ray)"]
    UNC["<b>Uncoated</b><br/><b>cross_boundary</b>: refraction, TIR,<br/>Fresnel for PolarizedRay"]
    LAYER["<b>thin layer</b><br/>interact3d(system, layer, t, beam, ray)<br/>e.g. ThinBeamsplitter: split,<br/>transmit! and reflect!"]
    OWN["<b>own interact3d</b><br/>e.g. mirror, detector"]
    ENC["<b>enclose_new!</b><br/>enclosure of each new ray at the hit:<br/>back into the medium before the boundary<br/>or into the medium behind it"]
    E([next tracing step])

    IA --> MS
    MS -- MultiShape --> PO --> K
    MS -- SingleShape --> K
    K -- "body" --> TR --> BODY --> UNC --> ENC
    K -- "thin layer" --> LAYER --> ENC
    TR -. "t" .-> LAYER
    K -- "own" --> OWN --> ENC
    ENC --> E

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class TR,UNC,ENC blue
    class BODY,LAYER,OWN green
    class MS,K purple
    class IA,E terminal
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

The physics at the boundary belongs to the boundary, not to one of the bodies. A body, e.g. any [`BeamletOptics.AbstractRefractiveOptic`](@ref), hands the ray to the [`BeamletOptics.surface_model`](@ref) of the transition, which is an uncoated surface so far. A thin layer implements the interaction between the media of the transition for itself:

```@docs; canonical=false
BeamletOptics.interact3d(::BeamletOptics.AbstractSystem, ::BeamletOptics.Uncoated, ::BeamletOptics.Transition, ::BeamletOptics.Beam, ::BeamletOptics.AbstractRay)
BeamletOptics.AbstractSurfaceModel
BeamletOptics.cross_boundary
```

A new refracting component declares its body via `medium` and needs no `interact3d` method of its own if it subtypes `AbstractRefractiveOptic`. A new thin layer implements `interact3d(system, layer::MyLayer, t::Transition, beam, ray)` and calls it from its `interact3d` method with `transition(system, ray)`, like the [`ThinBeamsplitter`](@ref).

The `interact3d` return type limits the interface to only accepting one new `beam` segment per interaction at the moment. After e.g. a lens surface air-to-glass interaction, the most likely next step is the refraction at the glass-to-air surface of the same lens. An interaction can tell the solver so with a `Hint`, which speeds up the next intersection.

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

```@docs
OpticsBase.PlaneField(::Union{GaussianBeamlet, AstigmaticGaussianBeamlet}, ::Integer)
```

## Hints

The [`BeamletOptics.Hint`](@ref) interface allows developers to tell the non-sequential solver which object will probably be hit during the next cycle of the [Intersect-Interact-Repeat-Loop](@ref). The solver tests this object first. Its hit lets the solver skip all objects whose bounding sphere lies beyond it, which saves most of the intersection tests in a large system.

```@docs; canonical=false
BeamletOptics.Hint
```

A `Hint` never changes the result of the trace: the hinted object is not preferred over a closer one, and an object that is moved away from its hinted neighbour, e.g. the back lens of a doublet, leaves an air gap that the trace sees. Touching elements need no `Hint`, see [Coincident hits](@ref) and [Transitions](@ref). The diagram shows how the next tracing step uses the hint of an interaction (purple: dispatch).

```mermaid
flowchart TB
    P["<b>push!(beam, interaction)</b><br/>next ray, hint(interaction)"]
    TST{"<b>tracing_step!</b><br/>(system, ray, hint)"}
    HS["<b>trace_one</b><br/>intersect3d(object(hint), ray) first"]
    OB["<b>closest_hits</b><br/>of all other objects, those<br/>beyond the hint's hit skipped"]
    TA["<b>trace_all</b><br/>closest_hits of all objects"]
    FI["<b>first_to_interact</b><br/>the hit that interacts first"]

    P --> TST
    TST -- "::Hint" --> HS --> OB --> FI
    TST -- "::Nothing" --> TA --> FI

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class P,HS,OB,TA,FI blue
    class TST purple
```

## Tracing logic

### Tracing systems

In the initial state, is is assumed that the problem consists of `objects <: AbstractObject` (in a system) and a `beam <: AbstractBeam` with a defined starting position and direction. No additional information is provided, and the specific path of the beam is not known beforehand. Consequently, brute force tracing of the optical system is required, involving testing against each individual element to determine the trajectory of the beam.

A call of [`solve_system!`](@ref) traces the beam tree breadth-first: it traces one beam after another and queues the children that the components attach, up to `depth_max` branching levels (green: implemented by a component author, purple: dispatch).

```mermaid
flowchart TB
    S(["<b>solve_system!</b><br/>system + beam"])
    INI["<b>initialize!(system)</b><br/>initialize!(object) of each object,<br/>if initialize = true"]
    BS["<b>with_bounding_spheres</b><br/>spheres of all objects,<br/>once per solve"]
    Q["<b>popfirst!(queue)</b><br/>next (beam, depth),<br/>the root at depth 1"]
    EM["<b>empty!(current)</b><br/>back to its first ray"]
    TS{"<b>trace_system!</b><br/>(system, current)"}
    D{"depth #lt;= depth_max?"}
    ENQ["<b>push!(queue, ...)</b><br/>children(current) at depth + 1,<br/>attached via children! or relaunch!"]
    DROP["<b>_drop_beams!(current)</b><br/>children are discarded"]
    QE{"queue empty?"}
    E([done])

    S -. "initialize = true" .-> INI
    S --> BS
    INI --> BS
    BS --> Q --> EM --> TS
    TS -- "method per beam type" --> D
    D -- yes --> ENQ --> QE
    D -- no --> DROP --> QE
    QE -- no --> Q
    QE -- yes --> E

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class BS,Q,EM,ENQ,DROP blue
    class INI green
    class TS purple
    class S,E terminal
```

```@docs; canonical=false
BeamletOptics.trace_system!
```

For a [`Beam`](@ref), `trace_system!` runs the loop below until the beam has `r_max` rays, its last ray hits nothing or `interact3d` returns `nothing`; the methods for beamlets run one `tracing_step!` per component ray and one `interact3d` for the whole beamlet per step, and also stop if the rays hit different shapes (an [`AstigmaticGaussianBeamlet`](@ref) also stops if the optical invariant check fails).

```mermaid
flowchart TB
    S(["<b>trace_system!</b><br/>system + beam::Beam"])
    RM{"length(rays(beam))<br/>#lt; r_max?"}
    TST{"<b>tracing_step!</b><br/>(system, last ray, hint)"}
    TO["<b>trace_one</b><br/>hinted object first,<br/>then all others"]
    TA["<b>trace_all</b><br/>intersect3d of all top-level<br/>objects and groups"]
    HIT{"intersection(ray)<br/>is nothing?<br/>else first_to_interact"}
    IA["<b>interact3d</b><br/>(system, object, beam, ray)"]
    NI{"interaction<br/>is nothing?"}
    EN["<b>enclose_new!</b><br/>enclosure of the new rays:<br/>next ray or children"]
    P["<b>push!(beam, interaction)</b><br/>next ray, hint(interaction)"]
    E([return])

    S --> RM
    RM -- yes --> TST
    TST -- "::Hint" --> TO
    TST -- "::Nothing" --> TA
    TO --> HIT
    TA --> HIT
    HIT -- no --> IA --> EN --> NI
    NI -- no --> P --> RM
    RM -- no --> E
    HIT -- yes --> E
    NI -- yes --> E

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class TO,TA,P,EN blue
    class IA green
    class TST purple
    class S,E terminal
```

This non-sequential mode is comparatively safe in determining the "true" beam path, but will scale suboptimally in time-complexity with the amount of optical elements.

!!! info "Object order"
    Unlike with classic, surface-based ray tracers, the order in which objects are listed in the [`System`](@ref) object vector/tuple is not considered for the purpose of tracing.

## CPU and GPU support

Parallizing the execution of a [`solve_system!`](@ref) call on the CPU is straight-forward for systems that do not feature objects which can be mutated during runtime, e.g. detectors like the [`Detector`](@ref). For each beam or ray the solution is independent and the solver can run on multiple threads. Special consideration needs to be taken when implementing mutable elements as mentioned above, since multiple threads might be able to access the underlying memory, leading to race conditions. Specifically, this means ensuring e.g. atomic write and read access.

With respect to GPU acceleration, this is not the case. Currently, all available implementations of [`solve_system!`](@ref) are highly branching algorithms which can not be implemented in a parallized way easily. This will most likley require a specific new subtype of the [`BeamletOptics.AbstractSystem`](@ref) with determinable sequential properties. This is not a development goal as of the writing of this section.

The beams of a group are solved in chunks, and all beams of a group share one table of bounding spheres. It belongs to a task, so each task of the threads sets it for its beams. Each chunk collects its detector hits and stores them at its end. A free task takes the next chunk, because the beams of a group differ in cost (e.g. the inner rings of a source pass a lens, the outer ones miss it), and those of similar cost lie next to each other. Four chunks per thread balance most of the load. Each further chunk costs a channel access and a flush, which shows if all beams cost the same.



