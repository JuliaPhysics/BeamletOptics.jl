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
render!(ax::_RenderEnv, obj::BMO.AbstractObject; kwargs...) = _render!(ax, obj; kwargs...)

# Dispatch helper fct. for RenderPresets.jl, do not remove
_render!(ax::_RenderEnv, obj::BMO.AbstractObject; kwargs...) = render!(ax, BMO.shape_trait_of(obj), obj; kwargs...)

function render!(ax::_RenderEnv, ::BMO.SingleShape, obj; kwargs...)
    render!(ax, BMO.shape(obj); kwargs...)
    return nothing
end

function render!(ax::_RenderEnv, ::BMO.MultiShape, obj; kwargs...)
    for _obj in BMO.shape(obj)
        render!(ax, _obj; kwargs...)
    end
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
