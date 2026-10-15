#=
The render handle protocol of the live rendering, see RenderLive.jl. Packages built on
BeamletOptics, e.g. a GUI, use only the protocol below; the concrete handles of the `Makie`
extension are internal.
=#

"""
    AbstractRenderHandle

Supertype of all handles returned by [`live_render!`](@ref). A handle references the rendered
object, system or beam and its plots, which can be re-synchronized via [`update_render!`](@ref)
and deleted via [`remove_render!`](@ref). Concrete handles are implemented by the `Makie`
extension.

# Render handle protocol

Every handle implements [`rendered`](@ref) and [`render_plots`](@ref). The subtypes add:

- [`AbstractObjectRenderHandle`](@ref): one movable thing, e.g. an object or a marker drawn via
  `live_render!(draw, ax, x)`
- [`AbstractSystemRenderHandle`](@ref): the object handles of a system, with its hierarchy of
  groups
- [`AbstractBeamRenderHandle`](@ref): a ray, beam or beam group, with [`render_settings`](@ref) and
  [`render_settings!`](@ref)

Code that uses handles, e.g. a GUI, relies on this protocol only, and may add own subtypes, e.g.
a system handle that combines the handles of several systems.
"""
abstract type AbstractRenderHandle end

"""
    AbstractObjectRenderHandle <: AbstractRenderHandle

Handle of one movable thing, returned by `live_render!(ax, obj)` for an `AbstractObject` or an
object group and by `live_render!(draw, ax, x)` for anything with a pose. The plots are drawn once;
[`update_render!`](@ref) applies the current pose of [`rendered`](@ref)`(h)` as their model
matrix (and draws them again if the parts of a `MultiShape` object moved relative to each other).
[`pick_object`](@ref) returns `rendered(h)` for the plots of [`pickable_plots`](@ref).
"""
abstract type AbstractObjectRenderHandle <: AbstractRenderHandle end

"""
    AbstractSystemRenderHandle <: AbstractRenderHandle

Handle of a system, returned by `live_render!(ax, sys)`: one [`AbstractObjectRenderHandle`](@ref)
per object, object groups rendered per object, such that each object of a group can be moved on
its own. A subtype implements

- [`rendered`](@ref)`(h)`: the system
- [`render_children`](@ref)`(h)`: the object handles
- [`render_parent`](@ref)`(h, obj)`: the group that holds `obj`, or `nothing` at the top level
- `push!(h, oh::AbstractObjectRenderHandle)`: adds the object handle `oh` at the top level, e.g. of
  an object added to the scene, and returns `h`
- `delete!(h, oh::AbstractObjectRenderHandle)`: removes `oh` from `h` and returns `h`; its plots
  stay, see [`remove_render!`](@ref)

[`render_plots`](@ref), [`update_render!`](@ref), [`remove_render!`](@ref) and
[`pick_object`](@ref) of a system handle are defined in terms of these.

A handle that knows its axis, like the one returned by `live_render!(ax, sys)`, also implements
`live_render!(h, obj; kwargs...)`, which renders an object (or object group) into `h` at runtime,
and extends `remove_render!(h, obj)` by the hierarchy of the groups.
"""
abstract type AbstractSystemRenderHandle <: AbstractRenderHandle end

"""
    AbstractBeamRenderHandle <: AbstractRenderHandle

Handle of a ray, beam or beam group, returned by `live_render!(ax, beam)`. Implements
[`rendered`](@ref), [`render_plots`](@ref), [`render_settings`](@ref) and
[`render_settings!`](@ref).
[`update_render!`](@ref) draws the current rays of the beam, e.g. after
[`solve_system!`](@ref). Keep the handle on a root beam or a beam group: child beams and all rays
after the first are new objects after every solve, so a handle on one of them does not follow the
solved beam.
"""
abstract type AbstractBeamRenderHandle <: AbstractRenderHandle end

"""
    rendered(h::AbstractRenderHandle)

The object, system or beam that the render handle `h` shows, see [`AbstractRenderHandle`](@ref).
"""
function rendered end

"""
    render_plots(h::AbstractRenderHandle) -> AbstractVector

The plots of the render handle `h`; of a system handle, the plots of all its object handles. The
plots of an object handle may be replaced by [`update_render!`](@ref), hence do not keep them.
"""
function render_plots end

render_plots(h::AbstractSystemRenderHandle) =
    reduce(vcat, (render_plots(c) for c in render_children(h)); init = Any[])

"""
    render_children(h::AbstractSystemRenderHandle) -> AbstractVector{<:AbstractObjectRenderHandle}

The object handles of the system handle `h`, one per rendered object (the objects of groups, not
the groups). Change them via `push!` and `delete!` of `h`, not via the returned vector.
"""
function render_children end

"""
    render_parent(h::AbstractSystemRenderHandle, obj)

The object group of the system of `h` that holds `obj` (an object or a group), or `nothing` if
`obj` is at the top level or not in `h`.
"""
function render_parent end

"""
    push!(h::AbstractSystemRenderHandle, oh::AbstractObjectRenderHandle) -> h

Adds the object handle `oh` at the top level of the system handle `h`, e.g. of an object added to
the scene after [`live_render!`](@ref), such that [`update_render!`](@ref) and
[`pick_object`](@ref) of `h` include it.
"""
Base.push!(::AbstractSystemRenderHandle, ::AbstractObjectRenderHandle)

"""
    delete!(h::AbstractSystemRenderHandle, oh::AbstractObjectRenderHandle) -> h

Removes the object handle `oh` from the system handle `h`. Its plots stay in the axis, delete
them via [`remove_render!`](@ref)`(oh)`.
"""
Base.delete!(::AbstractSystemRenderHandle, ::AbstractObjectRenderHandle)

"""
    render_settings(h::AbstractBeamRenderHandle) -> NamedTuple

The settings with which the beam handle `h` draws its beam, at least `flen` (length of a final
ray without intersection [m]) and `render_every` (every how many beams of a beam group are drawn,
`1` for other beams), e.g. to find the drawn segments of a beam, and `color` (the color of the rays
or of the envelope as it was given, e.g. `:blue` or `:wavelength`). The handles of Gaussian beamlets
and of beam groups of beamlets also report the resolution of the envelope, `r_res` and `z_res`.
Change them via [`render_settings!`](@ref).
"""
function render_settings end

"""
    render_settings!(h::AbstractBeamRenderHandle; kwargs...) -> h

Changes the settings of the beam handle `h` and draws the beam again with them, without creating
new plots: the plots of [`render_plots`](@ref)`(h)` stay in the axis, with their other attributes
(visibility, opacity, line width, clip planes, ...) unchanged. The keywords are among the keys of
[`render_settings`](@ref)`(h)`:

- `flen`: length of a final ray without intersection [m], positive and finite
- `render_every`: every how many beams of a beam group are drawn, a positive integer (no effect on
  other beams)
- `r_res`, `z_res`: radial and longitudinal resolution of an envelope mesh, integers of at least 2
- `color`: the color of the rays or of the envelope, a single color that `Makie` knows (e.g. `:red`,
  `(:red, 0.3)`), or `:wavelength` or `(:wavelength, alpha)` for the color of the wavelength of each
  ray, see [`wavelength_color`](@ref)

Any other keyword, e.g. `r_res` for a ray, throws an `ArgumentError`, and then `h` is not changed.
The beam is drawn as it is, i.e. nothing is solved, like [`update_render!`](@ref). The overlays of
the handle, e.g. `show_beams` or `show_pos`, follow the settings; the generating rays of
`show_beams` and the polarization curve keep their colors.

The plots of a handle that draw the `color` hold one color per vertex, which the handle sets
together with the positions. Change the color via `render_settings!`, not via the `color` attribute
of the plots, which the next update overwrites. A `color` that is neither a single color nor the
wavelength, e.g. a vector of colors, is passed on to `Makie` by `live_render!` and can not be
changed here.

```julia
h = live_render!(ax, beam; flen = 0.1)
render_settings!(h; flen = 0.5)             # the final ray is now 0.5 m long
render_settings!(h; color = :wavelength)    # each ray in the color of its wavelength
render_settings!(h; color = :orange)        # and all in one color again
```

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
render_settings!(::Any; kwargs...) = throw(MissingBackendError())

"""
    pickable_plots(x, plots) -> AbstractVector

The plots among the `plots` of `x` that select `x` in [`pick_object`](@ref), by default all.
Add a method for the type of `x` to exclude plots, e.g. the outline of a marker, whose clicks
should reach what lies behind it.
"""
pickable_plots(x, plots) = plots
