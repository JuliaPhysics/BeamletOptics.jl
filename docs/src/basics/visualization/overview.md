# Visualization

As mentioned in other sections of this documentation, the [Makie](https://docs.makie.org) backend can be used in order to generate 2D/3D renderings of optical systems and results generated with this package. Refer to the extensive `Makie` documentation and the **Examples** and the **Tutorials** sections of this package for a variety of showcases on how to visualize your simulation.

This section documents the [`render!`](@ref) function sorted by the type of the object being rendered:

- [Rays and beams](beams.md): single rays, beams of rays and beam groups
- [Gaussian beamlets](gaussian.md): stigmatic and astigmatic Gaussian beamlets
- [Components and systems](components.md): optical components, their default colours, polarizers, systems and the underlying shapes
- [Scene and camera](camera.md): helpers for framing and annotating a 3D scene

## Backends

The rendering functionality is provided by a package extension that is loaded together with a `Makie` backend:

- **GLMakie** is recommended. It opens interactive windows in which a 3D scene can be rotated and zoomed, supports both `LScene` and `Axis3`, and handles regular 2D plots (e.g. `lines`, `scatter`, `heatmap`) as well.
- **CairoMakie** is intended for static, high-quality PNG/PDF output. Only `Axis3` is supported for 3D renderings (see the [`render!`](@ref) docstring below).

## A first plot

A minimal plot needs a `Figure`, a 3D axis (`LScene` or `Axis3`) and one `render!` call per object that should be shown. Systems and beams are rendered the same way:

```julia
using GLMakie, BeamletOptics

mirror = RoundPlanoMirror(25e-3, 5e-3)
system = System([mirror])
beam = Beam([0, -0.1, 0], [0, 1.0, 0], 1e-6)
solve_system!(system, beam)

fig = Figure()
ax = LScene(fig[1, 1]) # or Axis3(fig[1, 1], aspect=:data)
render!(ax, system)
render!(ax, beam)
display(fig)
```

## How `render!` works

`render!` is a single generic function with many methods. The method is selected by the type of the second argument, i.e. the object that is rendered, so a [`Beam`](@ref) is drawn differently from a [`GaussianBeamlet`](@ref) or a lens. Keyword arguments and their defaults therefore differ between methods and are documented on the following pages next to the type they belong to.

- Typing `?render!` in the REPL prints every docstring of `render!` at once.
- Typing `?render!(ax, beam)`, with `ax` and `beam` defined, prints only the docstrings of the methods that match those argument types.
- `methods(render!)` lists all methods that are currently available.

Passing a renderable type for which no method exists throws a `RenderNotImplementedError`. Calling `render!` before a backend has been loaded throws a [`BeamletOptics.MissingBackendError`](@ref).

```@docs; canonical=false
render!(::Any, ::BeamletOptics._RenderTypes)
```

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
 [1] render!(::Nothing, ::Mirror{Float64, BeamletOptics.PlanoSurfaceSDF{Float64}})
   @ BeamletOptics c:\Users\anon\.julia\dev\BeamletOptics\src\Render.jl:46
 [2] top-level scope
   @ REPL[5]:1
```

Once a backend has been loaded, additional dispatched versions of `render!` become available.

```julia
julia> using GLMakie

julia> methods(render!)
# ... methods for generic function "render!" from BeamletOptics:
  [1] render!(axis::Union{Axis3, LScene}, gauss::GaussianBeamlet{T}; show_beams, show_pos, r_res, z_res, flen, color, transparency, kwargs...) where T
     @ BeamletOpticsMakieExt C:\Users\anon\.julia\dev\BeamletOptics\ext\RenderGaussian.jl:26
  [2] render!(axis::Union{Axis3, LScene}, beam_group::BeamletOptics.AbstractBeamGroup; render_every, kwargs...)
     @ BeamletOpticsMakieExt C:\Users\anon\.julia\dev\BeamletOptics\ext\RenderBeam.jl:152
  [3] render!(axis::Union{Axis3, LScene}, beam::Beam; flen, show_polarization, pol_λ, pol_amplitude, pol_ppl, pol_color, pol_linewidth, kwargs...)
     @ BeamletOpticsMakieExt C:\Users\anon\.julia\dev\BeamletOptics\ext\RenderBeam.jl:111
  ⋮
```
