"""
    _render_transmission_axis!(ax, pf::BMO.PolarizationFilter; color, linewidth, rim=nothing, rim_radius=nothing)

Draws the [`transmission_axis`](@ref) of `pf` as a film line across the filter, offset slightly along
the filter normal on both sides to remain visible regardless of viewing direction. If `rim = (a, b)` is
given (distances along `-normal`/`+normal`), also draws the two rim marks used by the [`LinearPolarizer`](@ref)
preset, from the front outer face to the back outer face at each end of the axis, at `rim_radius`
(defaults to the film extent along the axis).
"""
function _render_transmission_axis!(ax::_RenderEnv, pf::BMO.PolarizationFilter; color, linewidth, rim=nothing, rim_radius=nothing)
    t = BMO.transmission_axis(pf)
    n = BMO.orientation(pf)[:, 2]
    c = BMO.position(pf)
    verts = BMO.vertices(BMO.shape(pf))
    r = maximum(abs.((verts .- c') * t))
    δ = 1e-3 * r
    for s in (δ, -δ)
        p1 = c - r * t + s * n
        p2 = c + r * t + s * n
        lines!(ax, [p1[1], p2[1]], [p1[2], p2[2]], [p1[3], p2[3]]; color, linewidth, transparency=true)
    end
    if !isnothing(rim)
        a, b = rim
        ε = 5e-3
        for sgn in (1, -1)
            p = c + sgn * something(rim_radius, r) * (1 + ε) * t
            q1 = p - a * n
            q2 = p + b * n
            lines!(ax, [q1[1], q2[1]], [q1[2], q2[2]], [q1[3], q2[3]]; color, linewidth, transparency=true)
        end
    end
    return nothing
end

"""
    render!(ax, pf::PolarizationFilter; kwargs...)

Renders the [`PolarizationFilter`](@ref) `pf` into the specified `ax`, optionally with its transmission axis.

# Keyword args

- `show_transmission_axis = true`: draws the [`transmission_axis`](@ref) as a line on the filter
- `axis_color = :black`: color of the transmission axis line
- `axis_linewidth = 2`: line width of the transmission axis line

Remaining kwargs are passed on to the Makie mesh plot.
"""
function render!(ax::_RenderEnv, pf::PolarizationFilter; show_transmission_axis::Bool=true,
        axis_color=:black, axis_linewidth=2, kwargs...)
    _render!(ax, pf; kwargs...)
    show_transmission_axis && _render_transmission_axis!(ax, pf; color=axis_color, linewidth=axis_linewidth)
    return nothing
end

"""
    render!(ax, lipo::LinearPolarizer; kwargs...)

Renders the [`LinearPolarizer`](@ref) `lipo` into the specified `ax`, optionally with its transmission axis.
The two substrate halves render like refractive optics, and the filter film is drawn in translucent green.

# Keyword args

- `show_transmission_axis = true`: draws the [`transmission_axis`](@ref) as a line on the filter
- `axis_color = :black`: color of the transmission axis line
- `axis_linewidth = 2`: line width of the transmission axis line

Remaining kwargs are passed on to the Makie mesh plots of the substrate halves.
"""
function render!(ax::_RenderEnv, lipo::LinearPolarizer; show_transmission_axis::Bool=true,
        axis_color=:black, axis_linewidth=2, kwargs...)
    render!(ax, lipo.front; kwargs...)
    render!(ax, lipo.back; kwargs...)
    _render!(ax, lipo.filter; transparency=true, alpha=0.75, color=:green)
    if show_transmission_axis
        t_f = thickness(lipo.front)
        t_b = thickness(lipo.back)
        _render_transmission_axis!(ax, lipo.filter; color=axis_color, linewidth=axis_linewidth, rim=(t_f, t_b), rim_radius=BMO.diameter(BMO.shape(lipo.front)) / 2)
    end
    return nothing
end