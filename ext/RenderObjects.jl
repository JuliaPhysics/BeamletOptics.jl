function render!(ax::_RenderEnv, obj::BMO.AbstractObject; kwargs...)
    render!(ax, BMO.shape_trait_of(obj), obj; kwargs...)
    return nothing
end

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
