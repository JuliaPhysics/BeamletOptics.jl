# Systems & object groups

## System

```julia
System([obj1, obj2, group1, ...])   # flexible, Vector-backed
StaticSystem([obj1, obj2, ...])     # Tuple-backed: faster repeated solves, longer compile, fixed content
```

- Tracing is **non-sequential**: object order is irrelevant, and every beam segment is tested against all objects.
- `ObjectGroup`s are flattened (`BeamletOptics.objects(system)`).
- Prefer `System` while building and iterating and for large scenes. Use `StaticSystem` for small,
  fixed systems in tight loops.
- A `System` can be changed after construction (a `StaticSystem` can not):

```julia
system = System()                 # empty
push!(system, lens, mirror)       # objects or groups; an object already in the system throws
delete!(system, mirror)           # by identity, no-op if absent; an object inside a group throws
pop!(system)                      # removes and returns the last top-level object
popat!(system, i)                 # removes and returns the i-th top-level object (a group is one)
```

  Solve beams and sources again after a change: `solve_system!(system, beam)`.
  `System(v)` copies the vector `v`.

```julia
solve_system!(system, beam_or_source; r_max = 100, depth_max = 100,
              check_invariant = true, threshold = get_invariant_threshold())
```

Sources and other beam groups also take `progress = true`: a progress bar in the terminal once
tracing has run for `get_progress_threshold()` s (default 5 s).

## ObjectGroup

```julia
group = ObjectGroup([lens1, lens2, mount])   # members may be objects or other groups
```

- An `ObjectGroup` is a kinematic container: each member is traced and interacts on its own. Read
  the members with `BeamletOptics.objects(group)`. A component whose parts interact as one (e.g. a
  cemented doublet) is an `AbstractObject` subtype with the `MultiShape` trait.
- `position(group)` is the group `center`. It starts at the **origin**, not at the centroid;
  `orientation(group)` starts as identity.
- `translate3d!(group, Δ)` moves all members and the center.
- `translate_to3d!(group, p)` moves the members so that the center lands on `p`.
- `rotate3d!` and `x|y|zrotate3d!` rotate all members about the center.
- `set_pivot3d!(group, p)` moves the center (the pivot) without moving the members.
- Groups can be nested, e.g. for zoom lenses or interferometer arms with mounts.
- `push!(group, objs...)` and `delete!(group, obj)` change the top level of the group like those of
  a `System`: by identity, poses and pivot unchanged; an object already in the group (at any depth)
  or a mix of static and movable objects throws, so does `delete!` of an object inside a subgroup.
  `pop!(group)` and `popat!(group, i)` remove and return the last or the `i`-th top-level part.
  The group does not know its systems: `delete!(system, obj)` before `push!(group, obj)` if the
  system holds `obj` at its top level. Only `ObjectGroup` can be changed, not own group types.
  A `StaticSystem` keeps the objects its groups had when it was built: trace a changed group in a
  `System` or build a new `StaticSystem`. Added objects join the current frame of the group:
  `reset_rotation3d!`/`reset_translation3d!` of the group later move them by the group's whole
  rotation/translation so far. Do not change a group while a system that holds it is solved.

Typical pattern (see `templates/05_periscope_objectgroup.jl`):

```julia
periscope = ObjectGroup([m1, m2])
set_pivot3d!(periscope, [0, 0, 50mm])
zrotate3d!(periscope, deg2rad(10))
```

## Assemblies: group or one multi-shape object

Decide by whether the parts touch:

| Parts of the assembly | Use | `interact3d` |
|---|---|---|
| separate, with air between them (telescope, periscope, lens group, arm with mounts) | `ObjectGroup`, or an own subtype of `BeamletOptics.AbstractObjectGroup` | none: each part interacts on its own |
| touching and interacting as one (cemented, a `Hint` across a shared face, a shared medium) | an `AbstractObject` subtype with `shape_trait_of(::T) = MultiShape()` | required, for every beam model it must support |

- An own group type gives named parts and an own constructor. It needs the fields `objects`,
  `center::Point3{T}` and `dir::SMatrix{3,3,T,9}`, is a `mutable struct`, and its constructor
  checks `allequal(BeamletOptics.is_static, objects)`. Moving, `set_pivot3d!`, systems and rendering
  work as for an `ObjectGroup`. See `templates/08_custom_group.jl`.
- Never wrap separate parts in a `MultiShape` object without `interact3d`: the solver hands every
  hit to the wrapper, finds no method, warns `No interact3d method defined` and stops the ray.
- A `MultiShape` object finds the part that was hit via `shape(intersection(ray))` and forwards to
  it (see `CubeBeamsplitter`, `DoubletLens`). Beamlets reach this ray method through their generic
  `interact3d(system, obj, beamlet, id::Int)`, once per ray of the beamlet. If a part has its own
  beamlet method (detectors, beamsplitters), the object needs a beamlet method that forwards to it.
