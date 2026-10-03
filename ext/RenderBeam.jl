using Makie: Observable, AbstractPlot, lift

#=
One drawing path per beam type: a pure geometry function per beam type (segments of rays, mesh of a
Gaussian envelope, polarization curve) and one plot function per beam type, which draws from
`Observable`s of this geometry. `render!` calls the plot function and discards the observables,
`live_render!` keeps them in a `BeamRenderHandle`, whose `update_render!` refills them.
=#

"""
    _BeamPlots

Plots of a beam, as returned by the plot function of a beam type, e.g. [`_plot_rays!`](@ref). Each
entry of `fills` pairs an `Observable` that a plot draws from with the function that returns its
current geometry, see [`_observe!`](@ref).
"""
struct _BeamPlots
    plots::Vector{AbstractPlot}
    fills::Vector{Pair{Observable, Function}}
end

_BeamPlots() = _BeamPlots(AbstractPlot[], Pair{Observable, Function}[])

"""
    _observe!(bp::_BeamPlots, geometry, obs = Observable(geometry()))

Registers the observable `obs` of the geometry function `geometry()` in `bp`, such that
`update_render!` refills it, and returns `obs`.
"""
function _observe!(bp::_BeamPlots, geometry::Function, obs::Observable = Observable(geometry()))
    push!(bp.fills, obs => geometry)
    return obs
end

"""Adds the plot `p` to `bp` and returns it."""
_add!(bp::_BeamPlots, p::AbstractPlot) = (push!(bp.plots, p); p)

"""
    BeamRenderHandle{B} <: AbstractBeamRenderHandle

Live rendering handle of a ray, [`Beam`](@ref), beam group, [`GaussianBeamlet`](@ref) or
[`BMO.AstigmaticGaussianBeamlet`](@ref), see [`live_render!`](@ref). The plots are the ones that
`render!` draws, e.g. all ray segments in a single `linesegments` plot or the envelope in a single
`mesh`; [`update_render!`](@ref) refills the observables they draw from, the number of segments may
change between updates.
"""
mutable struct BeamRenderHandle{B} <: AbstractBeamRenderHandle
    thing::B
    axis::_RenderEnv
    plots::Vector{AbstractPlot}
    fills::Vector{Pair{Observable, Function}}
    settings::NamedTuple
end

BeamRenderHandle(thing, axis, bp::_BeamPlots, settings::NamedTuple) =
    BeamRenderHandle(thing, axis, bp.plots, bp.fills, settings)

function Base.show(io::IO, h::BeamRenderHandle{B}) where {B}
    print(io, "BeamRenderHandle{", B, "}(", length(h.plots), " plots)")
end

function update_render!(h::BeamRenderHandle)
    for (obs, geometry) in h.fills
        obs[] = geometry()
    end
    return h
end

function remove_render!(h::BeamRenderHandle)
    foreach(p -> delete!(h.axis, p), h.plots)
    empty!(h.plots)
    return nothing
end

rendered(h::BeamRenderHandle) = h.thing
render_plots(h::BeamRenderHandle) = h.plots
render_settings(h::BeamRenderHandle) = h.settings

"""
    _render_every(thing, render_every)

Every how many beams of `thing` are drawn: `render_every` for a beam group, `1` otherwise.
"""
_render_every(::BMO.AbstractBeamGroup, render_every) = render_every
_render_every(_, render_every) = 1

#=
Geometry of rays, beams and beam groups of `Beam`s
=#

"""
    _foreach_ray(f, thing, render_every = 1)

Calls `f(ray)` for all rays of a ray, [`Beam`](@ref) (the whole tree) or every `render_every`-th
beam of a beam group of `Beam`s.
"""
_foreach_ray(f, ray::BMO.AbstractRay, render_every = 1) = (f(ray); nothing)

function _foreach_ray(f, beam::Beam, render_every = 1)
    for child in PreOrderDFS(beam)
        for ray in BMO.rays(child)
            f(ray)
        end
    end
    return nothing
end

function _foreach_ray(f, bg::BMO.AbstractBeamGroup, render_every = 1)
    bms = BMO.beams(bg)
    if !isempty(bms) && !(first(bms) isa Beam)
        throw(ArgumentError("rendering via linesegments is only implemented for beam groups of " *
            "`Beam`s; groups of `AstigmaticGaussianBeamlet`s are rendered via the dedicated " *
            "`render!(axis, ::AstigmaticBeamGroup; ...)` mesh method"))
    end
    for i in 1:render_every:length(bms)
        _foreach_ray(f, bms[i])
    end
    return nothing
end

"""Pushes the start and the end point of the `ray` onto `pts`, `flen` long without intersection."""
function _push_ray_segment!(pts::Vector{Point3f}, ray::BMO.AbstractRay; flen)
    isect = BMO.intersection(ray)
    len = isnothing(isect) ? flen : length(isect)
    push!(pts, Point3f(position(ray)))
    push!(pts, Point3f(position(ray) + len * BMO.direction(ray)))
    return nothing
end

"""Pushes the start point of the `ray` onto `pts`, and its end point if it has an intersection."""
function _push_ray_ends!(pts::Vector{Point3f}, ray::BMO.AbstractRay)
    push!(pts, Point3f(position(ray)))
    isect = BMO.intersection(ray)
    isnothing(isect) || push!(pts, Point3f(position(ray) + length(isect) * BMO.direction(ray)))
    return nothing
end

"""
    _ray_segments(thing; flen, render_every = 1) -> Vector{Point3f}

Start and end points of all ray segments of `thing` (see [`_foreach_ray`](@ref)), as drawn by a
`linesegments` plot. A final ray without intersection is `flen` long.
"""
function _ray_segments(thing; flen, render_every::Int = 1)
    pts = Point3f[]
    _foreach_ray(ray -> _push_ray_segment!(pts, ray; flen), thing, render_every)
    return pts
end

"""
    _ray_ends(thing; render_every = 1) -> Vector{Point3f}

Start points of all rays of `thing` and the end points of those with an intersection, see
`show_pos`.
"""
function _ray_ends(thing; render_every::Int = 1)
    pts = Point3f[]
    _foreach_ray(ray -> _push_ray_ends!(pts, ray), thing, render_every)
    return pts
end

"""
    _ray_colors(thing, color::_ByWavelength; render_every = 1, ends = false) -> Vector{RGBAf}

The color of each vertex of [`_ray_segments`](@ref) (two per ray), or of [`_ray_ends`](@ref) (one or
two per ray) with `ends = true`: the color of the wavelength of its ray.
"""
function _ray_colors(thing, color::_ByWavelength; render_every::Int = 1, ends::Bool = false)
    cols = RGBAf[]
    _foreach_ray(thing, render_every) do ray
        c = _wavelength_rgba(BMO.wavelength(ray), color)
        push!(cols, c)
        (!ends || !isnothing(BMO.intersection(ray))) && push!(cols, c)
    end
    return cols
end

"""
    _segments!(bp, color, thing; flen, render_every) -> (points, color)
    _ends!(bp, color, thing; render_every) -> (points, color)

Registers the observable of the segments (or of the end points of `show_pos`) of `thing` in `bp`
and returns it with the `color` attribute of its plot: `color` itself, or for a
[`_ByWavelength`](@ref) an observable of one color per vertex, which is refilled together with the
points, such that both always have the same length.
"""
_segments!(bp::_BeamPlots, color, thing; flen, render_every) =
    (_observe!(bp, () -> _ray_segments(thing; flen, render_every)), color)

_ends!(bp::_BeamPlots, color, thing; render_every) = (_observe!(bp, () -> _ray_ends(thing; render_every)), color)

function _segments!(bp::_BeamPlots, color::_ByWavelength, thing; flen, render_every)
    both = _observe!(bp, () -> (_ray_segments(thing; flen, render_every), _ray_colors(thing, color; render_every)))
    return lift(first, both), lift(last, both)
end

function _ends!(bp::_BeamPlots, color::_ByWavelength, thing; render_every)
    both = _observe!(bp,
        () -> (_ray_ends(thing; render_every), _ray_colors(thing, color; render_every, ends = true)))
    return lift(first, both), lift(last, both)
end

"""
    _plot_rays!(axis, thing; kwargs...) -> _BeamPlots

The plot function of rays, [`Beam`](@ref)s and beam groups of `Beam`s: all segments in a single
`linesegments` plot, the positions of `show_pos` in a single `scatter` plot and the polarization
curve of `show_polarization` in a single `lines` plot. See `render!(axis, ray)` for the kwargs.
"""
function _plot_rays!(
        axis::_RenderEnv,
        thing;
        flen,
        render_every::Int,
        show_pos,
        show_polarization,
        pol_λ,
        pol_amplitude,
        pol_ppl,
        pol_color,
        pol_linewidth,
        color,
        linewidth,
        transparency,
        kwargs...
    )
    show_polarization && _check_polarized(thing)
    bp = _BeamPlots()
    coloring = _coloring(color)
    segments, segment_color = _segments!(bp, coloring, thing; flen, render_every)
    _add!(bp, linesegments!(axis, segments; color = segment_color, linewidth, transparency, kwargs...))
    if show_pos
        ends, end_color = _ends!(bp, coloring, thing; render_every)
        _add!(bp, scatter!(axis, ends; color = end_color))
    end
    if show_polarization
        curve = _observe!(bp,
            () -> _field_curve(thing; flen, render_every, λ_vis = pol_λ, amplitude = pol_amplitude, ppl = pol_ppl))
        _add!(bp, lines!(axis, curve; color = pol_color, linewidth = pol_linewidth))
    end
    return bp
end

"""The kwargs of `render!` and `live_render!` of rays, beams and beam groups, with their defaults."""
_ray_defaults(kwargs) = (;
    flen = 1.0,
    render_every = 5,
    show_pos = false,
    show_polarization = false,
    pol_λ = nothing,
    pol_amplitude = nothing,
    pol_ppl = 32,
    pol_color = :crimson,
    pol_linewidth = 2.0,
    color = :blue,
    linewidth = 1.0,
    transparency = true,
    kwargs...)

"""
    render!(axis, ray; kwargs...)

Renders a `ray` as a 3D line into the specified `axis`.

# Keyword args

- `flen = 1.0`: plotted length of the infinite ray in case of no intersection in [m]
- `show_pos = false`: marks the starting position of the `ray` (and its end point, if it has an
  intersection) with a point

# Polarization kwargs

- `show_polarization = false`: overlay the E-field curve for a `PolarizedRay` (throws an
  `ArgumentError` for a non-polarized ray)
- `pol_λ = nothing`: visualization wavelength [m], default = total plotted length / 20.
  This is a plotting parameter, not the physical ray wavelength; the curve shows the
  `t = 0` snapshot `Re{E⊥·exp(i·k·s)}` along the accumulated optical path `s`. Values
  finer than `plotted length / 2000` are clamped with a warning.
- `pol_amplitude = nothing`: curve amplitude [m] at the maximum |E⊥|, default `pol_λ/4`
- `pol_ppl = 32`: sample points per `pol_λ` along the curve
- `pol_color = :crimson`: field curve color
- `pol_linewidth = 2.0`: field curve line width

# Makie kwargs

- `color = :blue`: ray color. `color = :wavelength` draws each ray in the display color of its
  wavelength, see [`wavelength_color`](@ref), e.g. to show the rays of white light after a prism;
  `color = (:wavelength, 0.3)` sets the opacity
- `linewidth = 1.0`: ray line width
- `transparency = true`: ray transparency

Additional kwargs are passed into the `linesegments` plot of the ray.
"""
function render!(axis::_RenderEnv, ray::BMO.AbstractRay; kwargs...)
    _plot_rays!(axis, ray; _ray_defaults(kwargs)...)
    return nothing
end

"""
    render!(axis, beam; kwargs...)

Render the entire `beam` of rays into the specified 3D-`axis`. All ray segments are drawn by a
single `linesegments` plot.

# Keyword args

Refer to the plotting method of the `AbstractRay` for a list of keyword arguments.

# Polarization kwargs

- `show_polarization = false`: overlay the E-field curve for a `Beam{T, <:PolarizedRay}`
  (throws an `ArgumentError` for a beam of non-polarized rays). The curve is drawn once
  for the whole beam tree, not per ray, to preserve phase continuity.
- `pol_λ = nothing`: visualization wavelength [m], default = total plotted length / 20.
  This is a plotting parameter, not the physical ray wavelength; the curve shows the
  `t = 0` snapshot `Re{E⊥·exp(i·k·s)}` along the accumulated optical path `s`. Values
  finer than `plotted length / 2000` are clamped with a warning.
- `pol_amplitude = nothing`: curve amplitude [m] at the maximum |E⊥|, default `pol_λ/4`
- `pol_ppl = 32`: sample points per `pol_λ` along the curve
- `pol_color = :crimson`: field curve color
- `pol_linewidth = 2.0`: field curve line width
"""
function render!(axis::_RenderEnv, beam::Beam; kwargs...)
    _plot_rays!(axis, beam; _ray_defaults(kwargs)...)
    return nothing
end

"""
    render!(axis, beam_group; kwargs...)

Renders the [`BeamletOptics.AbstractBeamGroup`](@ref) into the specified `axis`. The ray segments of
all rendered beams are drawn by a single `linesegments` plot.

# Keywords arguments

- `render_every = 5`: renders only every e.g. fifth individual beam in the group

Refer to the plotting method of the `AbstractRay` for further keyword arguments. With
`show_polarization = true`, each rendered beam has its own field curve.
"""
function render!(axis::_RenderEnv, beam_group::BMO.AbstractBeamGroup; kwargs...)
    _plot_rays!(axis, beam_group; _ray_defaults(kwargs)...)
    return nothing
end

"""
    live_render!(axis, beam; kwargs...)

Live-renders a ray, [`Beam`](@ref) or beam group like `render!`, i.e. all ray segments in a single
`linesegments` plot, and returns a `BeamRenderHandle`, see [`live_render!`](@ref). The number of
segments may change between updates.

# Keyword args

- `flen = 1.0`: length of the final segment in case of no intersection [m]
- `render_every = 5`: renders only every e.g. fifth beam of a beam group

All other kwargs of `render!(axis, ray)` (`show_pos`, polarization and Makie kwargs) apply with
the same defaults.
"""
function live_render!(axis::_RenderEnv, thing::Union{BMO.AbstractRay, Beam, BMO.AbstractBeamGroup}; kwargs...)
    kw = _ray_defaults(kwargs)
    bp = _plot_rays!(axis, thing; kw...)
    settings = (; flen = Float64(kw.flen), render_every = _render_every(thing, kw.render_every))
    return BeamRenderHandle(thing, axis, bp, settings)
end
