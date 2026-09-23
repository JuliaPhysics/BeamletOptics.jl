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

"""
    AbstractRenderHandle

Supertype of all handles returned by [`live_render!`](@ref). A handle references the rendered
object or beam and its plots, which can be re-synchronized via [`update_render!`](@ref).
Concrete handles are implemented by the `Makie` extension.
"""
abstract type AbstractRenderHandle end

"""
    live_render!(axis, thing; kwargs...)

Renders `thing` into the `axis` like [`render!`](@ref), but returns an [`AbstractRenderHandle`](@ref)
that can be updated in place via [`update_render!`](@ref). Intended for animations and interactive
applications, e.g. moving components that require the system to be solved repeatedly.

- objects and systems: the geometry is generated once, kinematic changes are applied as a model
  transformation of the existing plots
- rays, beams and beam groups: all segments are bundled into a single plot
- `GaussianBeamlet`: the envelope of all segments is bundled into a single mesh

Keyword arguments are passed on as for [`render!`](@ref).

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
live_render!(::Any, ::_RenderTypes; kwargs...) = throw(MissingBackendError())

"""
    update_render!(handle)

Re-synchronizes the plots of the `handle` with the current state of the rendered object or beam,
e.g. after moving components or calling [`solve_system!`](@ref).

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
update_render!(::Any; kwargs...) = throw(MissingBackendError())

"""
    remove_render!(handle)

Deletes all plots of the `handle` from its axis.

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
remove_render!(::Any) = throw(MissingBackendError())

"""
    pick_object(handle, plot)

Returns the object of a live-rendered object or system `handle` that is visualized by the `plot`,
or `nothing` if the `plot` does not belong to the `handle`.

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
pick_object(::Any, ::Any) = throw(MissingBackendError())

"""
    kinematic_controls!(axis, handle; kwargs...)

Enables mouse and keyboard controls for moving and rotating the objects of a live-rendered system
`handle` within the `axis`. Returns a controller that can be removed via `close`.

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
kinematic_controls!(::Any, ::Any; kwargs...) = throw(MissingBackendError())
