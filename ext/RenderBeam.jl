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
    temp = position(ray) + flen * BMO.direction(ray)

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
        scatter!(axis, position(ray); color)
    end

    return nothing
end

"""
    render!(axis, segment::RaySegment; kwargs...)

Renders a `RaySegment` as a 3D line from the ray position to its intersection (or `flen` if no intersection).
"""
function render!(
        axis::_RenderEnv,
        segment::BMO.RaySegment;
        # kwargs
        flen = 1.0,
        show_pos = false,
        # Makie kwargs
        color = :blue,
        linewidth = 1.0,
        transparency = true,
        kwargs...
    )
    hit = BMO.intersection(segment)
    if isnothing(hit)
        len = flen
    else
        len = length(hit)
    end
    ray = segment.ray
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
        scatter!(axis, position(ray); color)
        if !isnothing(hit)
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
        if isempty(child.segments)
            render!(axis, child.head_ray; kwargs...)
        else
            for seg in child.segments
                render!(axis, seg; kwargs...)
            end
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