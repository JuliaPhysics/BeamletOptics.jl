function _collect_ray_segments!(
        pts::Vector{Point3f},
        scatter_pts::Vector{Point3f},
        ray::BMO.AbstractRay,
        flen::Real,
        show_pos::Bool
    )
    p1 = Point3f(position(ray)...)
    p2 = Point3f((position(ray) + flen * BMO.direction(ray))...)
    push!(pts, p1, p2)
    if show_pos
        push!(scatter_pts, p1)
    end
    return nothing
end

function _collect_ray_segments!(
        pts::Vector{Point3f},
        scatter_pts::Vector{Point3f},
        segment::BMO.RaySegment,
        flen::Real,
        show_pos::Bool
    )
    hit = BMO.intersection(segment)
    len = isnothing(hit) ? flen : length(hit)
    ray = segment.ray
    p1 = Point3f(position(ray)...)
    p2 = Point3f((position(ray) + len * BMO.direction(ray))...)
    push!(pts, p1, p2)
    if show_pos
        push!(scatter_pts, p1)
        if !isnothing(hit)
            push!(scatter_pts, p2)
        end
    end
    return nothing
end

function _collect_ray_segments!(
        pts::Vector{Point3f},
        scatter_pts::Vector{Point3f},
        beam::Beam,
        flen::Real,
        show_pos::Bool
    )
    for child in PreOrderDFS(beam)
        if isempty(child.segments)
            _collect_ray_segments!(pts, scatter_pts, child.head_ray, flen, show_pos)
        else
            for seg in child.segments
                _collect_ray_segments!(pts, scatter_pts, seg, flen, show_pos)
            end
        end
    end
    return nothing
end

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
    pts = Point3f[]
    scatter_pts = Point3f[]
    _collect_ray_segments!(pts, scatter_pts, ray, flen, show_pos)

    linesegments!(axis, pts;
        color,
        linewidth,
        transparency,
        kwargs...
    )

    if show_pos && !isempty(scatter_pts)
        scatter!(axis, scatter_pts; color)
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
    pts = Point3f[]
    scatter_pts = Point3f[]
    _collect_ray_segments!(pts, scatter_pts, segment, flen, show_pos)

    linesegments!(axis, pts;
        color,
        linewidth,
        transparency,
        kwargs...
    )

    if show_pos && !isempty(scatter_pts)
        scatter!(axis, scatter_pts; color)
    end

    return nothing
end

"""
    render!(axis, beam; kwargs...)

Render the entire `beam` of rays into the specified 3D-`axis`.

# Keyword args

- `flen = 1.0`: plotted length of the infinite ray in case of no intersection in [m]
- `show_pos = false`: marks the starting position (and intersections) with spheres

# Makie kwargs

- `color = :blue`: ray color
- `linewidth = 1.0`: ray line width
- `transparency = true`: ray transparency

Additional kwargs can be passed into the line plot.
"""
function render!(
        axis::_RenderEnv,
        beam::Beam;
        # kwargs
        flen = 1.0,
        show_pos = false,
        # Makie kwargs
        color = :blue,
        linewidth = 1.0,
        transparency = true,
        kwargs...
    )
    pts = Point3f[]
    scatter_pts = Point3f[]
    _collect_ray_segments!(pts, scatter_pts, beam, flen, show_pos)

    if !isempty(pts)
        linesegments!(axis, pts;
            color,
            linewidth,
            transparency,
            kwargs...
        )
    end

    if show_pos && !isempty(scatter_pts)
        scatter!(axis, scatter_pts; color)
    end

    return nothing
end

"""
    render!(axis, beam_group; kwargs...)

Renders the [`BeamletOptics.AbstractBeamGroup`](@ref) into the specified `axis`.

# Keywords arguments

- `render_every = 5`: renders only every e.g. fifth individual beam in the group
- `flen = 1.0`: plotted length of the infinite ray in case of no intersection in [m]
- `show_pos = false`: marks the starting position (and intersections) with spheres

# Makie kwargs

- `color = :blue`: ray color
- `linewidth = 1.0`: ray line width
- `transparency = true`: ray transparency

Additional kwargs can be passed into the line plot.
"""
function render!(
        axis::_RenderEnv,
        beam_group::BMO.AbstractBeamGroup;
        # kwargs
        render_every::Int = 5,
        flen = 1.0,
        show_pos = false,
        # Makie kwargs
        color = :blue,
        linewidth = 1.0,
        transparency = true,
        kwargs...
    )
    pts = Point3f[]
    scatter_pts = Point3f[]
    bms = BMO.beams(beam_group)
    numEl = length(bms)
    for i = 1:render_every:numEl
        _collect_ray_segments!(pts, scatter_pts, bms[i], flen, show_pos)
    end

    if !isempty(pts)
        linesegments!(axis, pts;
            color,
            linewidth,
            transparency,
            kwargs...
        )
    end

    if show_pos && !isempty(scatter_pts)
        scatter!(axis, scatter_pts; color)
    end

    return nothing
end