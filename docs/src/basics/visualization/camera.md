# Scene and camera helpers

Alongside `render!`, a small set of `LScene`-specific helpers is provided for framing and
annotating a 3D scene once a backend is loaded: [`get_view`](@ref), [`set_view`](@ref),
[`set_orthographic`](@ref), [`hide_axis`](@ref), [`look_at!`](@ref), [`arrow!`](@ref) and
[`render_lcs!`](@ref). Like `render!`, each throws a
[`BeamletOptics.MissingBackendError`](@ref) if called before a suitable backend has been
loaded. Refer to the **Reference** page for their full docstrings.

To freeze a camera position found interactively, use `get_view` and `set_view`; to derive the
viewpoint from the scene's own geometry, use `set_view` with eye, look-at and up, or `look_at!`.

!!! note
    These helpers need an `LScene`. They do not work with an `Axis3`.

```@docs; canonical=false
get_view
set_view
set_orthographic
hide_axis
look_at!
arrow!
render_lcs!
```
