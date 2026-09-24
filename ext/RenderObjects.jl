"""
    render!(ax, object; kwargs...)

Renders the `object` into the specified `ax`is. Additional `kwargs` can be piped through to the backend.

# Examples

It is recommended to use the following snippet in order to generate plots:

```julia
using GLMakie, BeamletOptics

fig = Figure()
ax = LScene(fig[1,1]) # or Axis3
render!(ax, my_BMO_obj; color=:white)
```

Additional keyword arguments can be passed. Refer to the `Makie` and `BeamletOptics` documentation for supported options
for each `object`.
"""
render!(ax::_RenderEnv, object::BMO.AbstractObject; kwargs...) = _render!(ax, object; kwargs...)

# Dispatch helper fct. for RenderPresets.jl, do not remove
_render!(ax::_RenderEnv, obj::BMO.AbstractObject; kwargs...) = render!(ax, BMO.shape_trait_of(obj), obj; kwargs...)

"""
    render!(ax, ::SingleShape, obj; kwargs...)

Renders the shape of `obj` as one mesh if it has an analytic tessellation (see `_has_mesh`),
otherwise via the `render!` method of the shape. The `show_normals` kwargs of
`render!(ax, ::AbstractMesh)` also select the latter.
"""
function render!(ax::_RenderEnv, ::BMO.SingleShape, obj; kwargs...)
    s = BMO.shape(obj)
    if _has_mesh(s) && !haskey(kwargs, :show_normals) && !haskey(kwargs, :show_normals_length)
        _render_mesh!(ax, s; kwargs...)
    else
        render!(ax, s; kwargs...)
    end
    return nothing
end

"""
    render!(ax, ::MultiShape, obj; kwargs...)

Renders all parts of `obj`. The analytic meshes of the parts are merged into one mesh plot per set
of plot attributes (e.g. color), see `_plot_collected!`. Nested `MultiShape` parts are merged
into the meshes of the outermost object.
"""
function render!(ax::_RenderEnv, ::BMO.MultiShape, obj; kwargs...)
    if !isnothing(_MESH_COLLECTOR[])
        for _obj in BMO.shape(obj)
            render!(ax, _obj; kwargs...)
        end
        return nothing
    end
    parts = Any[]
    with(_MESH_COLLECTOR => parts) do
        for _obj in BMO.shape(obj)
            render!(ax, _obj; kwargs...)
        end
    end
    _plot_collected!(ax, parts)
    return nothing
end

"""
    render!(ax::_RenderEnv, sys::AbstractSystem)

Render all objects contained in the `sys`tem.
"""
function render!(ax::_RenderEnv, sys::BMO.AbstractSystem; kwargs...)
    # Avoid use of objects(sys)
    for _obj in sys.objects
        render!(ax, _obj; kwargs...)
    end
    return nothing
end
