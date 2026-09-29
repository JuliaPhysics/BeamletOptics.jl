abstract type RenderException <: Exception end

message(e::RenderException) = e.msg
showerror(io::IO, e::RenderException) = print(io, message(e))

"""
    MissingBackendError

Custom `Exception` type that indicates that the `Makie` extension of `BeamletOptics` has not been loaded correctly.
"""
mutable struct MissingBackendError <: RenderException
    msg::String
    function MissingBackendError()
        msg = "It appears no suitable Makie backend is loaded in this session."
        return new(msg)
    end
end

"""A collection of all types from BMO which might be renderable in principle."""
const _RenderTypes = Union{
    AbstractRay,
    AbstractBeam,
    AbstractShape,
    AbstractObject,
    AbstractObjectGroup,
    AbstractSystem,
}

"""
    render!(axis, thing; kwargs...)

The `render!` function allows for the visualization of optical system and beams under the condition that a suitable backend is loaded.
This means that either one of the following packages must be loaded in combination with `BeamletOptics` via `using`:

1. GLMakie 
    - preferred for 3D viewing
    - use `LScene` or `Axis3` environments
2. CairoMakie 
    - preferred for the generation of high-quality .pngs
    - only `Axis3` is supported

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.

# Implementations reqs.

All concrete implementations of `render!` must adhere to the following minimal interface:

`render!(axis, thing; kwargs...)`

- `axis`: an axis type of the union of `LScene` or `Axis3`
- `thing`: an abstract or concrete object or beam type
- `kwargs`: custom or `Makie` keyword arguments that are passed to the underlying backend

Refer to the `BeamletOptics` extension docs for `Makie` for more information.
"""
render!(::Any, ::_RenderTypes, kwargs...) = throw(MissingBackendError())

"""
    get_view(ls)

Returns the current camera view matrix of an `LScene` `ls`, if a suitable backend is loaded.

Handy at the REPL to freeze a view found interactively: rotate the scene by hand, call
`get_view(ax)`, and paste the printed matrix into the script as a literal passed to
[`set_view`](@ref).

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
get_view(::Any) = throw(MissingBackendError())

"""
    set_view(ls, view::AbstractMatrix)
    set_view(ls, eye, lookat, up)

Sets the camera of an `LScene` `ls`, if a suitable backend is loaded, either directly from
a view matrix (e.g. one obtained via [`get_view`](@ref)) or from an eye/lookat/up triple.
See also [`look_at!`](@ref) to aim the camera at a known point instead of specifying the
triple directly.

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
set_view(::Any, ::AbstractMatrix) = throw(MissingBackendError())
set_view(::Any, ::Any, ::Any, ::Any) = throw(MissingBackendError())

"""
    set_orthographic(ls)

Switches the camera of an `LScene` `ls` to an orthographic projection, if a suitable
backend is loaded. If not, a [`MissingBackendError`](@ref) will be thrown.
"""
set_orthographic(::Any) = throw(MissingBackendError())

"""
    hide_axis(ls, hide::Bool=true)

Hides the axis markers of an `LScene` `ls`, if a suitable backend is loaded. Can be
toggled via `hide`. If no suitable backend is loaded, a [`MissingBackendError`](@ref)
will be thrown.
"""
hide_axis(::Any, ::Bool) = throw(MissingBackendError())

"""
    arrow!(ax, pos, dir; scale=1, kwargs...)

Draws a single 3D arrow from `pos` pointing along `dir` into `ax`, if a suitable backend
is loaded, scaled to a fixed on-screen length (independent of `dir`'s own norm) so it
stays legible next to CAD geometry. If no suitable backend is loaded, a
[`MissingBackendError`](@ref) will be thrown.
"""
arrow!(::Any, ::AbstractVector, ::AbstractVector; kwargs...) = throw(MissingBackendError())

"""
    look_at!(ax, target, offset; up = [0, 0, 1])

Aims the camera of `ax` at `target` from `target + offset`, if a suitable backend is
loaded. A deterministic replacement for manually orbiting the scene to find a viewpoint,
handy for reproducible close-up figures. If no suitable backend is loaded, a
[`MissingBackendError`](@ref) will be thrown.
"""
look_at!(::Any, ::AbstractVector, ::AbstractVector; kwargs...) = throw(MissingBackendError())

"""
    render_lcs!(ax, pos, lcs; scale = 10, show_labels = false)
    render_lcs!(ax, object; scale = 10, show_labels = false)

Draws the local coordinate system of an object (or of an explicit `pos`/orientation pair)
into `ax` as a red/green/yellow arrow triad, if a suitable backend is loaded. Useful to
make the reference frame of an imported CAD mesh visible in the scene. If no suitable
backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
render_lcs!(::Any, ::AbstractArray = zeros(3), ::AbstractMatrix = Matrix{Float64}(I, 3, 3); kwargs...) =
    throw(MissingBackendError())
render_lcs!(::Any, ::AbstractObject; kwargs...) = throw(MissingBackendError())


#=
Live rendering and the render handle protocol. Packages built on BeamletOptics, e.g. a GUI, use
only the protocol below; the concrete handles of the `Makie` extension are internal.
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
- [`AbstractBeamRenderHandle`](@ref): a ray, beam or beam group, with [`render_settings`](@ref)

Code that uses handles, e.g. a GUI, relies on this protocol only, and may add own subtypes, e.g.
a system handle that combines the handles of several systems.
"""
abstract type AbstractRenderHandle end

"""
    AbstractObjectRenderHandle <: AbstractRenderHandle

Handle of one movable thing, returned by `live_render!(ax, obj)` for an `AbstractObject` and by
`live_render!(draw, ax, x)` for anything with a pose. The plots are drawn once;
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
"""
abstract type AbstractSystemRenderHandle <: AbstractRenderHandle end

"""
    AbstractBeamRenderHandle <: AbstractRenderHandle

Handle of a ray, beam or beam group, returned by `live_render!(ax, beam)`. Implements
[`rendered`](@ref), [`render_plots`](@ref) and [`render_settings`](@ref).
[`update_render!`](@ref) draws the current rays of the beam, e.g. after
[`solve_system!`](@ref).
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
`1` for other beams), e.g. to find the drawn segments of a beam.
"""
function render_settings end

"""
    pickable_plots(x, plots) -> AbstractVector

The plots among the `plots` of `x` that select `x` in [`pick_object`](@ref), by default all.
Add a method for the type of `x` to exclude plots, e.g. the outline of a marker, whose clicks
should reach what lies behind it.
"""
pickable_plots(x, plots) = plots

"""
    live_render!(axis, thing; kwargs...)

Renders `thing` into the `axis` like [`render!`](@ref), but returns an
[`AbstractRenderHandle`](@ref) that can be updated in place via [`update_render!`](@ref). Intended
for animations and interactive applications, e.g. moving components that require the system to be
solved repeatedly.

- objects and systems: the geometry is generated once, kinematic changes are applied as a model
  transformation of the existing plots
- rays, beams and beam groups: all segments are bundled into a single plot
- `GaussianBeamlet`: the envelope of all segments is bundled into a single mesh

Keyword arguments are passed on as for [`render!`](@ref).

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
live_render!(::Any, ::_RenderTypes; kwargs...) = throw(MissingBackendError())

"""
    live_render!(draw, axis, x) -> AbstractObjectRenderHandle

Live-renders `x` via the function `draw`: `draw()` plots `x` in its current pose into the `axis`,
and the returned handle moves these plots with `x` like those of an object, see
[`live_render!`](@ref), e.g. a marker of a thing that `render!` does not draw. `x` has a `position`
and an `orientation` (or a `direction`), e.g. a source or an own movable type (see
[`kinematic_trait_of`](@ref)):

```julia
h = live_render!(ax, src) do
    scatter!(ax, [Point3f(position(src))]; color = :orange)
end
```

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
live_render!(::Function, ::Any, ::Any) = throw(MissingBackendError())

"""
    update_render!(handle)

Re-synchronizes the plots of the `handle` with the current state of the rendered object or beam,
e.g. after moving components or calling [`solve_system!`](@ref).

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
update_render!(::Any; kwargs...) = throw(MissingBackendError())

function update_render!(h::AbstractSystemRenderHandle)
    foreach(update_render!, render_children(h))
    return h
end

"""
    remove_render!(handle)

Deletes all plots of the `handle` from its axis.

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
remove_render!(::Any) = throw(MissingBackendError())

function remove_render!(h::AbstractSystemRenderHandle)
    foreach(remove_render!, render_children(h))
    return nothing
end

"""
    pick_object(handle, plot)

Returns the object of a live-rendered object or system `handle` that is visualized by the `plot`
(one of its [`pickable_plots`](@ref), or a child plot of one), or `nothing` if the `plot` does not
belong to the `handle`. For a system handle, the top-level object, i.e. the outermost group of
the picked object, see [`render_parent`](@ref).

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
pick_object(::Any, ::Any) = throw(MissingBackendError())

function pick_object(h::AbstractSystemRenderHandle, plot)
    for child in render_children(h)
        obj = pick_object(child, plot)
        isnothing(obj) && continue
        parent = render_parent(h, obj)
        while !isnothing(parent)
            obj = parent
            parent = render_parent(h, obj)
        end
        return obj
    end
    return nothing
end

"""
    look_colors() -> Dict{Symbol, RGBf}

The colors of the material classes (`:refractive`, `:reflective`, `:coating`, `:polarizer`,
`:detector`, `:mechanics`, `:interface`) of the active look, see [`set_render_look`](@ref), e.g. to
recolor the rendered objects of a class for a dark background.

Needs the `Makie` extension, i.e. a loaded Makie backend.
"""
function look_colors end

"""
    studio_lighting!(ax::Union{LScene, Axis3}; preset = :studio)

Sets up a CAD-like lighting rig in the 3D view `ax`, if a suitable backend is loaded: an ambient
light, a key light from the upper right front, a fill light from the left and a rim light from
behind, all relative to the camera. Backends with a single directional light (e.g. CairoMakie)
get the ambient and the key light only. `preset = :none` leaves the lights unchanged.
The GUI of the package BeamletOpticsGUI applies the rig by default, scenes created via [`render!`](@ref) call it
explicitly.

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
studio_lighting!(::Any; kwargs...) = throw(MissingBackendError())

"""
    set_render_look(look::Symbol)

Sets the look of all subsequently rendered objects, `:modern` (default) or `:cad`. The `:modern`
look renders clear glass with faint silhouettes, metallic mirrors and neutral mechanics without
edge lines, the `:cad` look saturated materials with feature edge lines. Explicit kwargs of [`render!`](@ref), e.g.
`color`, `material` or `edges`, override the look.

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
set_render_look(::Any) = throw(MissingBackendError())
