"""
    render!(axis, ray; kwargs...)

Renders a `ray` as a 3D line into the specified `axis`.

# Keyword args

- `flen = 1.0`: plotted length of the infinite ray in case of no intersection in [m]
- `show_pos = false`: marks the starting position of the `ray` with a sphere

# Polarization kwargs

- `show_polarization = false`: overlay the E-field curve for a `PolarizedRay` (throws an
  `ArgumentError` for a non-polarized ray)
- `pol_λ = nothing`: visualization wavelength [m], default = total plotted length / 20.
  This is a plotting parameter, not the physical ray wavelength; the curve shows the
  `t = 0` snapshot `Re{E⊥·exp(i·k·s)}` along the accumulated optical path `s`.
- `pol_amplitude = nothing`: curve amplitude [m] at the maximum |E⊥|, default `pol_λ/4`
- `pol_ppl = 32`: sample points per `pol_λ` along the curve
- `pol_color = :crimson`: field curve color
- `pol_linewidth = 2.0`: field curve line width

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
        # Polarization kwargs
        show_polarization = false,
        pol_λ = nothing,
        pol_amplitude = nothing,
        pol_ppl = 32,
        pol_color = :crimson,
        pol_linewidth = 2.0,
        # Makie kwargs
        color = :blue,
        linewidth = 1.0,
        transparency = true,
        kwargs...
    )
    show_polarization && _check_polarized(ray)

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

    if show_polarization
        pts = _polarization_points(ray; flen, λ_vis = pol_λ, amplitude = pol_amplitude, ppl = pol_ppl)
        _render_field_curve!(axis, pts; color = pol_color, linewidth = pol_linewidth)
    end

    return nothing
end

"""
    render!(axis, beam; kwargs...)

Render the entire `beam` of rays into the specified 3D-`axis`.

# Keyword args

Refer to the plotting method of the `AbstractRay` for a list of keyword arguments.

# Polarization kwargs

- `show_polarization = false`: overlay the E-field curve for a `Beam{T, <:PolarizedRay}`
  (throws an `ArgumentError` for a beam of non-polarized rays). The curve is drawn once
  for the whole beam tree, not per ray, to preserve phase continuity.
- `pol_λ = nothing`: visualization wavelength [m], default = total plotted length / 20.
  This is a plotting parameter, not the physical ray wavelength; the curve shows the
  `t = 0` snapshot `Re{E⊥·exp(i·k·s)}` along the accumulated optical path `s`.
- `pol_amplitude = nothing`: curve amplitude [m] at the maximum |E⊥|, default `pol_λ/4`
- `pol_ppl = 32`: sample points per `pol_λ` along the curve
- `pol_color = :crimson`: field curve color
- `pol_linewidth = 2.0`: field curve line width
"""
function render!(
        axis::_RenderEnv,
        beam::Beam;
        # kwargs
        flen = 1.0,
        # Polarization kwargs
        show_polarization = false,
        pol_λ = nothing,
        pol_amplitude = nothing,
        pol_ppl = 32,
        pol_color = :crimson,
        pol_linewidth = 2.0,
        kwargs...
    )
    show_polarization && _check_polarized(beam)

    for child in PreOrderDFS(beam)
        for ray in BMO.rays(child)
            render!(axis, ray; flen, kwargs...)
        end
    end

    if show_polarization
        pts = _polarization_points(beam; flen, λ_vis = pol_λ, amplitude = pol_amplitude, ppl = pol_ppl)
        _render_field_curve!(axis, pts; color = pol_color, linewidth = pol_linewidth)
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