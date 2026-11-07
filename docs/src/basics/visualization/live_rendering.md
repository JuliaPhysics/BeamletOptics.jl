# Live rendering

`render!` generates new plots on every call. For animations and interactive applications, where
components are moved and the system is re-solved many times per second, use
[`live_render!`](@ref) instead. It returns a handle that re-synchronizes the existing plots with
the current simulation state via [`update_render!`](@ref):

```julia
using GLMakie, BeamletOptics

fig = Figure()
ax = LScene(fig[1, 1])

hsys = live_render!(ax, system)   # geometry is generated once
hbeam = live_render!(ax, beam)    # all ray segments bundled into a single plot

zrotate3d!(mirror, 1e-3)
solve_system!(system, beam)
update_render!(hsys)              # moved components: only their model transformation changes
update_render!(hbeam)             # beam path: point buffer is replaced in place
```

```@docs; canonical=false
BeamletOptics.live_render!(::Any, ::BeamletOptics._RenderTypes)
BeamletOptics.update_render!(::Any)
BeamletOptics.remove_render!(::Any)
```

`live_render!` draws the same plots as `render!`, see [Rays and beams](@ref "Rendering rays and beams") and [Gaussian beamlets](@ref "Rendering Gaussian beamlets"); it only keeps the data they are drawn from. The envelope meshes of Gaussian beamlets use coarser default resolutions `r_res` and `z_res` than `render!` (24/40 instead of 50/100 for a `GaussianBeamlet`, 10/8 instead of 64/100 per beamlet of an `AstigmaticBeamGroup`). Change them via the keywords or with [`BeamletOptics.render_settings!`](@ref).

Which `live_render!` method draws a thing, which handle it returns and how `update_render!` refreshes it are described in [Live rendering dispatch](@ref) and [Updating a handle](@ref).

## Changing a live scene

An object that is added to a system at runtime, see [Changing a system](@ref), is rendered into the
handle of the system, and removed from it again, by:

```@docs; canonical=false
BeamletOptics.live_render!(::BeamletOptics.AbstractSystemRenderHandle, ::BeamletOptics.ObjectOrGroup)
BeamletOptics.remove_render!(::BeamletOptics.AbstractSystemRenderHandle, ::BeamletOptics.ObjectOrGroup)
```

```julia
push!(system, lens)
live_render!(hsys, lens)     # draws the lens, hsys now updates and picks it
delete!(system, lens)
remove_render!(hsys, lens)   # deletes its plots and its handle
```

## Changing the shape of an object

`update_render!` moves the plots of an object but keeps its geometry. After the shape of a rendered
object changed, e.g. the leg length of a prism of an own component type, draw it again with
`rerender!`, called on the handle of the object or of the whole system:

```@docs; canonical=false
BeamletOptics.rerender!(::Any)
```

```julia
hprism = live_render!(ax, prism)
translate3d!(prism, [0, 0.01, 0])
update_render!(hprism)   # moved: the plots follow the new pose
# ... the shape of prism changes
rerender!(hprism)        # new plots of the new shape, at the current pose
```

Beam handles report how they were drawn, and can change it at runtime, e.g. the drawn length of the
final rays or the color (a single color or the color of the wavelength of each ray), without
creating new plots:

```@docs; canonical=false
BeamletOptics.render_settings
BeamletOptics.render_settings!
```

Overlays that are not part of a component, e.g. markers, are drawn with the function form of
`live_render!`, which returns a handle for the plots `draw` creates and moves them together with `x`:

```@docs; canonical=false
BeamletOptics.live_render!(::Function, ::Any, ::Any)
```

The object that owns a picked plot is found by `pick_object`:

```@docs; canonical=false
BeamletOptics.pick_object(::Any, ::Any)
```

Packages built on BeamletOptics access the handles through the [Render handle protocol](@ref).
