"""
    render_ray!(axis, ray::AbstractRay; color=:blue, flen=1.0)

Renders a `ray` as a 3D line. If the ray has no intersection, the substitute length `flen` is used.
"""
function render!(ax::_RenderEnv, ray::BMO.AbstractRay; color = :blue, flen = 1.0, show_pos=false)
    if isnothing(BMO.intersection(ray))
        len = flen
    else
        len = length(BMO.intersection(ray))
    end
    temp = BMO.position(ray) + len * BMO.direction(ray)

    lines!(ax,
        [BMO.position(ray)[1], temp[1]],
        [BMO.position(ray)[2], temp[2]],
        [BMO.position(ray)[3], temp[3]],
        color=color,
        linewidth=1.0,
        transparency=true)
    if show_pos
        scatter!(ax, ray.pos; color)
    end

    return nothing
end

"""
    render_beam!(axis, beam::Beam; color=:blue, flen=1.0, show_pos=false)

Render the entire `beam` into the specified 3D-`axis`. A `color` can be specified.
"""
function render!(ax::_RenderEnv, beam::Beam; color = :blue, flen = 1.0, show_pos=false)
    for child in PreOrderDFS(beam)
        for ray in BMO.rays(child)
            render!(ax, ray; color, flen, show_pos)
        end
    end
    return nothing
end