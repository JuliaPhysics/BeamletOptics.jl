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
    I["<b>1. Intersect</b><br/>hinted shape first,<br/>else all objects"] -- hit --> X["<b>2. Interact</b><br/>compute<br/>AbstractInteraction"]
    X --> R["<b>3. Repeat</b><br/>attach next<br/>ray or beam segment"]
    R -- "optional Hint" --> I
    I -- "no hit or r_max" --> E([done])

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class I blue
    class X green
    class R purple
    class S,E terminal
```

This scheme is loosely referred to as the **Intersect-Interact-Repeat-Loop** and consists of the following steps:

1. Calculate the closest **Intersection** between a ray/beam and the objects within the system
2. Calculate the optical **Interaction** that occurs at the surface or within the volume of the element
3. Attach the next part of the ray chain or beam tree
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

The following diagram shows how `trace_all` reaches the shapes of the objects and groups of a system, including the [Bounding sphere pretest](@ref) (green: implemented by a component author, purple: dispatch).

```mermaid
flowchart TB
    TA(["<b>trace_all(system, ray)</b>"])
    BSO["<b>bounding_sphere_of(shape)</b><br/><i>optional</i>"]
    TBL["<b>bounding_spheres(system)</b><br/>table, once per solve"]
    P{"<b>intersect3d(parts, ray)</b><br/>each part, closest hit wins"}
    G["<b>intersect3d(group, ray)</b><br/>x = objects(group)"]
    OJ["<b>intersect3d(object, ray)</b><br/><i>own method possible,<br/>e.g. PlateBeamsplitter</i>"]
    T{"<b>shape_trait_of(object)</b>"}
    SS["<b>SingleShape</b><br/>x = shape(object)"]
    MS["<b>MultiShape</b><br/>x = shape(object), a tuple"]
    PT{"<b>intersect3d(sphere, x, ray)</b><br/>sphere of x from the table,<br/>else NoBoundingSphere"}
    SH["<b>intersect3d(shape, ray)</b><br/>shape backend, see below"]
    OBJ["<b>object!(intersection, object)</b><br/>by the SingleShape and MultiShape<br/>methods, not by a group"]
    MISS(["nothing"])

    BSO --> TBL
    TBL -. "current_bounding_spheres()" .-> PT
    TA --> P
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
    class P,T,PT purple
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
BeamletOptics.interact3d(::BeamletOptics.AbstractSystem, ::BeamletOptics.AbstractObject, ::BeamletOptics.AbstractBeam, ::BeamletOptics.AbstractRay)
```

As with the `intersect3d` method, a predefined return type must be provided in order to make the [`solve_system!`](@ref) interface work. The `AbstractInteraction` is used in order to create "building blocks" from which the output beam is constructed.

```@docs; canonical=false
BeamletOptics.AbstractInteraction
```

An object whose `interact3d` method stores data in the object, e.g. a detector its hits, also implements [`initialize!`](@ref), which discards this data before a new solve.

The `interact3d` return type limits the interface to only accepting one new `beam` segment per interaction at the moment. The developer needs to take into account that after e.g. a lens surface air-to-glass interaction, the solver "forgets" that the next logical step is to immediately test against the lens again, since the most likely step will be the refraction at the glass-to-air surface. In order to alleviate this issue, the `Hint` type can be used.

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

As mentioned in the previous section, the [`BeamletOptics.Hint`](@ref) interface allows developers to manipulate the non-sequential solver algorithm into testing against a specific component and shape during the next cycle of the [Intersect-Interact-Repeat-Loop](@ref). This interface has very high priority during intersection testing.

```@docs; canonical=false
BeamletOptics.Hint
```

The main reason for this is the intersection ambiguity encountered at interfaces between air-tight component interfaces, e.g. [Plate beamsplitters](@ref) or cemented [Doublet lenses](@ref). This is caused by the fact that for a ray with a starting point at this interface, technically both shapes are being "touched" at the same time. Additional program logic considering the ray direction of propagation can not always resolve this ambiguity. Therefore the task of providing additional information to the solver via the `Hint` interface is placed as a **burden on the developer**.

The diagram below shows this for the [`CubeBeamsplitter`](@ref), whose `interact3d` method finds the hit part by identity, hands the beam to it and hints the coating; the `DoubletLens` hints the other lens in the same way (green: implemented by a component author, purple: dispatch).

```mermaid
flowchart TB
    S(["<b>interact3d(system, cbs, beam, ray)</b><br/>object of the intersection: cbs,<br/>shape: the part that was hit"])
    F{"shape(intersection(ray))<br/>=== shape of which part?"}
    PR["<b>interact3d(system, cbs.front, beam, ray)</b><br/>the prism refracts<br/>(cbs.back alike)"]
    H["<b>hint!(interaction,</b><br/><b>Hint(cbs, shape(cbs.coating)))</b>"]
    CO["<b>interact3d(system, cbs.coating, beam, ray)</b><br/>children!: transmitted and reflected,<br/>refractive_index! of the cube"]
    N(["return nothing<br/>the solver queues the children"])
    P["<b>push!(beam, interaction)</b>"]
    TST{"<b>tracing_step!</b><br/>(system, ray, hint)"}
    HS{"<b>trace_one</b><br/>hinted shape hit?"}
    OB["<b>object!(intersection, object(hint))</b><br/>no other object is tested"]
    TA["<b>trace_all(system, ray)</b>"]

    S --> F
    F -- "cbs.front or cbs.back" --> PR --> H --> P --> TST
    F -- "cbs.coating" --> CO --> N
    TST -- "::Hint" --> HS
    HS -- yes --> OB
    HS -- no --> TA

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class P,OB,TA,HS blue
    class F,H,PR,CO green
    class TST purple
    class S,N terminal
```

!!! tip
    Use of the `Hint` interface is primarily intended for multishape objects with joint surfaces.

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
    TO["<b>trace_one</b><br/>hinted shape first"]
    TA["<b>trace_all</b><br/>intersect3d of all top-level<br/>objects and groups"]
    HIT{"intersection(ray)<br/>is nothing?"}
    IA["<b>interact3d</b><br/>(system, object, beam, ray)"]
    NI{"interaction<br/>is nothing?"}
    P["<b>push!(beam, interaction)</b><br/>next ray, hint(interaction)"]
    E([return])

    S --> RM
    RM -- yes --> TST
    TST -- "::Hint" --> TO
    TST -- "::Nothing" --> TA
    TO -- "hinted shape missed" --> TA
    TO -- hit --> HIT
    TA --> HIT
    HIT -- no --> IA --> NI
    NI -- no --> P --> RM
    RM -- no --> E
    HIT -- yes --> E
    NI -- yes --> E

    classDef blue fill:#4063D826,stroke:#4063D8,stroke-width:2px
    classDef green fill:#38982626,stroke:#389826,stroke-width:2px
    classDef purple fill:#9558B226,stroke:#9558B2,stroke-width:2px
    classDef terminal fill:transparent,stroke:#CB3C33,stroke-width:1.5px,stroke-dasharray:4 3
    class TO,TA,P blue
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



