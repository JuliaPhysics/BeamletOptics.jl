"""
    render!(axis, ray; kwargs...)

Renders a `ray` as a 3D line into the specified `axis`.

# Keyword args

- `flen`: plotted length of the infinite ray in case of no intersection in [m], default is 1.0
- `show_pos`: marks the starting position of the `ray` with a sphere, default is false

# Makie kwargs

- `color`: default is :blue
- `linewidth`: default is 1.0
- `transparency`: default is true

Additional kwargs can be passed into the line plot.
"""
function render!(
        axis::_RenderEnv,
        ray::BMO.AbstractRay;
        # kwargs
        flen = 1.0,
        show_pos = false,
        # Makie kwargs
        color = :blue,
        linewidth = 1.0,
        transparency = true,
        kwargs...
    )
    if isnothing(BMO.intersection(ray))
        len = flen
    else
        len = length(BMO.intersection(ray))
    end
    temp = BMO.position(ray) + len * BMO.direction(ray)

    lines!(axis,
        [BMO.position(ray)[1], temp[1]],
        [BMO.position(ray)[2], temp[2]],
        [BMO.position(ray)[3], temp[3]];
        color,
        linewidth,
        transparency,
        kwargs...
    )
    if show_pos
        scatter!(axis, ray.pos; color)
    end

    return nothing
end

"""
    render!(axis, beam; kwargs...)

Render the entire `beam` of rays into the specified 3D-`axis`.

# Keyword args

Refer to the plotting method of the `AbstractRay` for a list of keyword arguments.
"""
function render!(
        axis::_RenderEnv,
        beam::Beam;
        kwargs...
    )
    for child in PreOrderDFS(beam)
        for ray in BMO.rays(child)
            render!(axis, ray; kwargs...)
        end
    end
    return nothing
end