# Developer API

Packages built on BeamletOptics, e.g. the interactive
[BeamletOpticsGUI](https://github.com/StackEnjoyer/BeamletOpticsGUI.jl), need a few functions that are
not part of the user-facing API. They are declared `public` (see `Base.ispublic`) but not exported:
access them as `BeamletOptics.name`. They are stable within a minor version like the exported API.
The names of the render handle protocol are documented in
[Render handle protocol](@ref "Render handle protocol").

## Kinematic queries

Whether an object can be moved is decided by its kinematic trait; `is_static` is the query used by
tools that offer to move objects.

```@docs; canonical=false
BeamletOptics.is_static
```

## Detectors and sources

Used by tools that display detector counts and that regenerate the rays of a source at a lower or
higher count.

```@docs; canonical=false
BeamletOptics.hit_count
BeamletOptics.source_wavelength
BeamletOptics.min_num_rays
```

## Progress and cancellation

Used by tools that show the progress of a long solve and let the user cancel it.

```@docs; canonical=false
BeamletOptics.ProgressSink
BeamletOptics.PROGRESS_SINK
BeamletOptics.progress_state
BeamletOptics.is_cancelled
```
