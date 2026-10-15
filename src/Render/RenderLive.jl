# Live rendering via the render handles of RenderHandles.jl, implemented in ext/RenderLive.jl and
# the beam files of the extension

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
    live_render!(h::AbstractSystemRenderHandle, obj; kwargs...) -> Vector{<:AbstractObjectRenderHandle}

Live-renders `obj`, an object or object group, into the axis of the system handle `h` and adds it
to `h`, e.g. an object that was added to the system via `push!(system, obj)` after `h` was created.
An object group is rendered
per object, with its hierarchy known to [`render_parent`](@ref), like the groups of
`live_render!(ax, system)`. Returns the new object handles, one per rendered object. The system of
`h` is not changed. An `obj` that `h` already renders throws an `ArgumentError`.

Keyword arguments are passed on as for [`render!`](@ref). Implemented by the handle returned by
`live_render!(ax, system)`. A system handle that does not implement it throws an `ArgumentError`:
render `obj` via `live_render!(ax, obj)` and add its handle via `push!(h, oh)` instead.
"""
function live_render!(h::AbstractSystemRenderHandle, ::ObjectOrGroup; kwargs...)
    throw(ArgumentError(
        "a $(nameof(typeof(h))) can not render objects, use live_render!(ax, obj) and push!(h, oh) instead"))
end

"""
    update_render!(handle)

Re-synchronizes the plots of the `handle` with the current state of the rendered object or beam,
e.g. after moving components or calling [`solve_system!`](@ref). Keep render handles on root beams and
groups: child beams and all rays after the first are new objects after every solve.

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
    remove_render!(h::AbstractSystemRenderHandle, obj)

Removes `obj`, an object or object group, from the system handle `h`: deletes the plots of `obj`
(of all objects of an object group) from the axis and removes their object handles from `h`, e.g. after `delete!(system, obj)`.
The system of `h` is not changed. Nothing happens for an `obj` that `h` does not render. An object
within a group can not be removed on its own and throws an `ArgumentError`: remove the group
instead.
"""
function remove_render!(h::AbstractSystemRenderHandle, obj::ObjectOrGroup)
    parent = render_parent(h, obj)
    isnothing(parent) || throw(ArgumentError(
        "the $(nameof(typeof(obj))) is part of a $(nameof(typeof(parent))) of the system handle, remove the group instead"))
    leaves = collect(Leaves(obj))
    for oh in filter(oh -> _is_leaf_of(rendered(oh), leaves), render_children(h))
        remove_render!(oh)
        delete!(h, oh)
    end
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
