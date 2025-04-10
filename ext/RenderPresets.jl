# function render_dummy_mesh!(axis::_RenderEnv, d::NonInteractableObject; transparency = false, kwargs...)
#     mesh = d.shape
#     mesh!(axis, vertices(mesh), faces(mesh); transparency, color = :grey, kwargs...)
#     return nothing
# end

render!(ax::_RenderEnv, dummy::BMO.NonInteractableObject; transparency=false, color=:grey, kwargs...) = render!(
    ax,
    BMO.shape(dummy);
    transparency,
    color,
    kwargs...
)