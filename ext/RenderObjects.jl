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

function render!(ax::_RenderEnv, ::BMO.SingleShape, obj; kwargs...)
    render!(ax, BMO.shape(obj); kwargs...)
    return nothing
end

function render!(ax::_RenderEnv, ::BMO.SingleShape, lens::BMO.Lens; color = RGBf(0.678, 0.847, 0.902), kwargs...)
    s = BMO.shape(lens)
    if s isa BMO.UnionSDF && !isempty(BMO.coatings(lens))
        for sub_sdf in s.sdfs
            if sub_sdf isa BMO.ConcaveSphericalSurfaceSDF || sub_sdf isa BMO.ConvexSphericalSurfaceSDF || sub_sdf isa BMO.AbstractAsphericalSurfaceSDF
                sub_opt_axis = orientation(sub_sdf)[:, 2]
                lens_opt_axis = orientation(lens)[:, 2]
                local_n = if dot(sub_opt_axis, lens_opt_axis) < 0
                    SVector{3, Float64}(0.0, -1.0, 0.0) # front surface
                else
                    SVector{3, Float64}(0.0, 1.0, 0.0)  # back surface
                end
                c_model = BMO.get_matching_coating(BMO.coatings(lens), s, [0.0, 0.0, 0.0], local_n)
                sub_color = coating_color(c_model, color)
                render!(ax, sub_sdf; color = sub_color, kwargs...)
            else
                render!(ax, sub_sdf; color = color, kwargs...)
            end
        end
    else
        render!(ax, s; color = color, kwargs...)
    end
    return nothing
end

function render!(ax::_RenderEnv, ::BMO.SingleShape, mirror::BMO.Mirror; color = :silver, kwargs...)
    if !isempty(BMO.coatings(mirror))
        c_model = BMO.get_matching_coating(BMO.coatings(mirror), BMO.shape(mirror), [0.0, 0.0, 0.0], [0.0, -1.0, 0.0])
        m_color = c_model isa BMO.Uncoated ? color : coating_color(c_model, color)
        render!(ax, BMO.shape(mirror); color = m_color, kwargs...)
    else
        render!(ax, BMO.shape(mirror); color = color, kwargs...)
    end
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
