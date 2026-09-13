# Visualization

As mentioned in other sections of this documentation, the [Makie](https://docs.makie.org) backend can be used in order to generate 2D/3D renderings of optical systems and results generated with this package. Refer to the extensive `Makie` documentation and the **Examples** and the **Tutorials** sections of this package for a variety of showcases on how to visualize your simulation.

## Rendering elements

The main function provided for visualization purposes is the [`render!`](@ref) function. 

```@docs; canonical=false
render!(::Any, ::BeamletOptics._RenderTypes)
```

If a suitable backend is loaded, additional dispatched `render!` functions will become available. For instance, this allows the plotting of a [`GaussianBeamlet`](@ref).

## Loading the extension

Refer to the following snippet for an example on how the extension loading behaves. When only BMO is loaded, the `render!` function becomes available but will throw an [`BeamletOptics.MissingBackendError`](@ref) when trying to plot something.

```julia
julia> using BeamletOptics

julia> methods(render!)
# 1 method for generic function "render!" from BeamletOptics:
 [1] render!(::Any, ::Union{BeamletOptics.AbstractSystem, BeamletOptics.AbstractBeam, BeamletOptics.AbstractObject, BeamletOptics.AbstractObjectGroup, BeamletOptics.AbstractRay, BeamletOptics.AbstractShape}, kwargs...)
     @ C:\Users\anon\.julia\dev\BeamletOptics\src\Render.jl:56

julia> axis = nothing;

julia> mirror = RoundPlanoMirror(25e-3, 5e-3);

julia> render!(axis, mirror)
ERROR: It appears no suitable Makie backend is loaded in this session.
Stacktrace:
 [1] render!(::Nothing, ::RoundPlanoMirror{Float64})
   @ BeamletOptics c:\Users\anon\.julia\dev\BeamletOptics\src\Render.jl:46
 [2] top-level scope
   @ REPL[5]:1
```

Once a backend has been loaded, additional dispatched versions of `render!` become available.

```julia
julia> using GLMakie

julia> methods(render!)
# 21 methods for generic function "render!" from BeamletOptics:
  [1] render!(ax::Union{Axis3, LScene}, s::BeamletOptics.UnionSDF; kwargs...)
     @ BeamletOpticsMakieExt C:\Users\anon\.julia\dev\BeamletOptics\ext\RenderSDF.jl:32
  [2] render!(axis::Union{Axis3, LScene}, css::BeamletOptics.ConcaveSphericalSurfaceSDF; color, kwargs...)
     @ BeamletOpticsMakieExt C:\Users\anon\.julia\dev\BeamletOptics\ext\RenderLenses.jl:1
  [3] render!(axis::Union{Axis3, LScene}, css::BeamletOptics.ConvexSphericalSurfaceSDF; color, kwargs...)
     @ BeamletOpticsMakieExt C:\Users\anon\.julia\dev\BeamletOptics\ext\RenderLenses.jl:31
  [4] render!(axis::Union{Axis3, LScene}, acyl::BeamletOptics.AbstractAcylindricalSurfaceSDF; color, kwargs...)
     @ BeamletOpticsMakieExt C:\Users\anon\.julia\dev\BeamletOptics\ext\RenderCylinderLenses.jl:1
  [5] etc...
```

## Camera and scene helpers

Alongside `render!`, a small set of `LScene`-specific helpers is provided for framing and
annotating a 3D scene once a backend is loaded: [`get_view`](@ref), [`set_view`](@ref),
[`set_orthographic`](@ref), [`hide_axis`](@ref), [`look_at!`](@ref), [`arrow!`](@ref) and
[`render_lcs!`](@ref). Like `render!`, each throws a
[`BeamletOptics.MissingBackendError`](@ref) if called before a suitable backend has been
loaded. Refer to the **Reference** page for their full docstrings.

The `get_view`/`set_view(ls, matrix)` pair is meant for interactive use: rotate the scene
by hand, call `get_view(ax)`, and paste the printed matrix back into the script as a
literal passed to `set_view`. This is the pattern used throughout this package's own
tutorials to freeze a camera position found interactively. The `set_view(ls, eye, lookat,
up)` and `look_at!` forms are the reproducible alternative, useful when the viewpoint
should be derived from the scene's own geometry instead of copy-pasted.