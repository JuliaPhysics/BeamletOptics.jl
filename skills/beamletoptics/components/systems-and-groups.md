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
delete!(system, mirror)           # by identity; an object inside a group throws, remove the group
pop!(system)                      # removes and returns the last top-level object
popat!(system, i)                 # removes and returns the i-th top-level object (a group is one)
```

  Beams solved before a change keep their old path: `empty!(beam)`, then
  `solve_system!(system, beam; retrace = false)`.

```julia
solve_system!(system, beam_or_source; r_max = 100, retrace = true, depth_max = 100,
              check_invariant = true, threshold = get_invariant_threshold())
```

Sources and other beam groups also take `progress = true`: a progress bar in the terminal once
tracing has run for `get_progress_threshold()` s (default 5 s).

## ObjectGroup

```julia
group = ObjectGroup([lens1, lens2, mount])   # members may be objects or other groups
```

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
