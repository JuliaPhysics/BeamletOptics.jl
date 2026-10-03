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

  Beams and sources solved before a change keep their old path: `empty!(beam)`, then
  `solve_system!(system, beam)`. `System(v)` copies the vector `v`.

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

## Clear aperture (vignetting)

```julia
D = clear_aperture(system, pos, dir; λ = 1e-6)   # largest collimated bundle diameter [m] without vignetting
vignetted(source, axial_beam)                    # indices of the vignetted beams of a traced source
```

A ray is vignetted if it hits another sequence of objects than the axial ray (found by tracing, no
object order assumed) or leaves the system early. Collimated input only; sampling-based
(`rings`, `azimuths`), assumes a disc-shaped clear region (no central obstruction). Throws if the
axial ray hits nothing; returns `d_max` (default 1 m) if nothing vignettes. Use `0.9D` for a
vignetting-free source.
