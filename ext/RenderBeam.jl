"""
    render!(axis, ray; kwargs...)

Renders a `ray` as a 3D line into the specified `axis`.

# Keyword args

- `flen = 1.0`: plotted length of the infinite ray in case of no intersection in [m]
- `show_pos = false`: marks the starting position of the `ray` with a sphere

# Makie kwargs

- `color = :blue`: ray color
- `linewidth = 1.0`: ray line width
- `transparency = true`: ray transparency

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
    temp = position(ray) + len * BMO.direction(ray)

    lines!(axis,
        [position(ray)[1], temp[1]],
        [position(ray)[2], temp[2]],
        [position(ray)[3], temp[3]];
        color,
        linewidth,
        transparency,
        kwargs...
    )
    
    if show_pos
        # start point
        scatter!(axis, ray.pos; color)

        # end point
        if !isnothing(BMO.intersection(ray))
            scatter!(axis, temp; color)
        end
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

"""
    render!(axis, beam_group; kwargs...)

Renders the [`BeamletOptics.AbstractBeamGroup`](@ref) into the specified `axis`.

# Keywords arguments

- `render_every = 5`: renders only every e.g. fifth individual beam in the group

Refer to the plotting method of the `AbstractRay` for further keyword arguments.
"""
function render!(
        axis::_RenderEnv,
        beam_group::BMO.AbstractBeamGroup;
        # kwargs
        render_every::Int=5,
        # Makie kwargs
        kwargs...
    )
    numEl = length(BMO.beams(beam_group))
    for i = 1:render_every:numEl
        render!(axis, BMO.beams(beam_group)[i]; kwargs...)
    end
    return nothing
end