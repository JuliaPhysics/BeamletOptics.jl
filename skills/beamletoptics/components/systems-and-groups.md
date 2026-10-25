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
  the members with `BeamletOptics.objects(group)`. Members that touch (e.g. two plates against each
  other) are traced correctly: the solver resolves the shared face. A component with named parts and
  a fixed pivot (e.g. a cemented doublet) is an `AbstractObject` subtype with the `MultiShape` trait.
- `position(group)` is the group `center`. It starts at the **origin**, not at the centroid;
  `orientation(group)` starts as identity.
- `translate3d!(group, Δ)` moves all members and the center.
- `translate_to3d!(group, p)` moves the members so that the center lands on `p`.
- `rotate3d!` and `x|y|zrotate3d!` rotate all members about the center.
- `set_pivot3d!(group, p)` moves the center (the pivot) without moving the members.
- Groups can be nested, e.g. for zoom lenses or interferometer arms with mounts.

Typical pattern (see `templates/05_periscope_objectgroup.jl`):

```julia
periscope = ObjectGroup([m1, m2])
set_pivot3d!(periscope, [0, 0, 50mm])
zrotate3d!(periscope, deg2rad(10))
```

## Assemblies: group or one multi-shape object

Touching parts need no special handling: at a shared face the solver collects the surfaces of all
touching parts and derives the media on both sides (`BeamletOptics.transition`). Decide by how the
assembly is used:

| Parts of the assembly | Use | `interact3d` |
|---|---|---|
| a kinematic unit of objects (telescope, periscope, lens group, arm with mounts, also touching plates) | `ObjectGroup`, or an own subtype of `BeamletOptics.AbstractObjectGroup` | none: each part interacts on its own |
| one component with named parts (a cemented doublet, a cube beamsplitter) | an `AbstractObject` subtype with `shape_trait_of(::T) = MultiShape()` and objects as parts | none: the solver hands each hit to the part that was hit (`BeamletOptics.part_of`) |
| one component built from plain shapes | as above, with shapes as parts | required, for every beam model it must support |

- An own group type gives named parts and an own constructor. It needs the fields `objects`,
  `center::Point3{T}` and `dir::SMatrix{3,3,T,9}`, is a `mutable struct`, and its constructor
  checks `allequal(BeamletOptics.is_static, objects)`. Moving, `set_pivot3d!`, systems and rendering
  work as for an `ObjectGroup`. See `templates/08_custom_group.jl`.
- A `MultiShape` object of plain shapes without `interact3d` stops every ray that hits it with the
  warning `No interact3d method defined`. Its method finds the part that was hit via
  `shape(intersection(ray))`.
- A `MultiShape` object of objects forwards rays and beamlets to the part that was hit, also the
  beamlet methods of the part (detectors, beamsplitters).
- Media and boundaries (not exported, may still change): a body declares its medium with
  `BeamletOptics.medium(x::MyBody) = BeamletOptics.IsotropicMedium(n)` (every `AbstractRefractiveOptic`
  has one). In an own `interact3d`, `t = BeamletOptics.transition(system, ray)` holds the media on both
  sides of the hit (`BeamletOptics.refractive_index_in(t)`, `BeamletOptics.refractive_index_out(t)` at
  the ray's wavelength) and the thin layer at it; `BeamletOptics.interact3d(system,
  BeamletOptics.surface_model(t), t, beam, ray)` refracts or totally reflects at the boundary. A
  `Hint` only speeds up the trace.
- Rays carry the bodies that enclose them (`BeamletOptics.enclosure(ray)`), so bodies within bodies
  work (a lens in a water cuvette). The solver assigns the enclosure to the new rays of an own
  `interact3d` (continued ray and children starting at the hit); nothing to pass on. A source that
  starts inside a body is taken to be in the medium of the system.
