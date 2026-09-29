#=
Detector panels of the live view: specs, construction, update and the panel options on the
cards of the detectors
=#

# Number of updates shown in the history of a detector panel
const _HISTORY_LENGTH = 300

# Floor of the logarithmic intensity scale, relative to the maximum
const _LOG_FLOOR = 1e-4

"""
    DetectorPanel

Detector panel of a `LiveView`, which shows the spot diagram or the intensity of a
`Detector` in a 2D `Axis`, see `live_view`. The metrics of the last update are shown in the
subtitle of the axis and stored in `metrics`, the centroid is marked by a cross. Optionally, a
`history` axis shows the power (or number of hits) and the centroid over the last updates, and a
`profiles` axis the intensity along x and z through the centroid.
"""
mutable struct DetectorPanel
    pd::BMO.Detector
    name::String
    # :auto, :spot or :intensity
    mode::Symbol
    # kwargs of `intensity`, i.e. without the panel options below
    kwargs::NamedTuple
    # :linear or :log, fixed color range or `nothing`
    colorscale::Symbol
    colorrange::Any
    ax::Axis
    xy::Observable{Vector{Point2f}}
    heat_x::Observable{Vector{Float32}}
    heat_y::Observable{Vector{Float32}}
    heat_I::Observable{Matrix{Float32}}
    scatter_plot::AbstractPlot
    heat_plot::AbstractPlot
    # metrics of the last update and centroid cross [mm]
    metrics::Any
    centroid::Observable{Vector{Point2f}}
    centroid_plot::AbstractPlot
    # history of the power (or number of hits) and of the centroid [mm] over the update index, in
    # two axes with a common x-axis, empty without history
    history_axes::Vector{Axis}
    history_count::Int
    history_value::Observable{Vector{Point2f}}
    history_cx::Observable{Vector{Point2f}}
    history_cz::Observable{Vector{Point2f}}
    # intensity profiles along x and z through the centroid [mm]
    profiles_ax::Union{Nothing, Axis}
    profile_x::Observable{Vector{Point2f}}
    profile_z::Observable{Vector{Point2f}}
    # last error of the update, logged only once
    last_error::Union{Nothing, String}
end

Base.show(io::IO, p::DetectorPanel) = print(io, "DetectorPanel(", p.name, ", mode = ", p.mode, ")")

"""Returns all `Detector`s of the `systems`, deduplicated by identity, in the order of discovery."""
function _find_detectors(systems)
    seen = Base.IdSet{Any}()
    dets = BMO.Detector[]
    for sys in systems, obj in BMO.objects(sys)
        obj isa BMO.Detector || continue
        obj in seen && continue
        push!(seen, obj)
        push!(dets, obj)
    end
    return dets
end

const _PANEL_MODES = (:auto, :spot, :intensity)

_panel_spec(pd::BMO.Detector) = (pd, :auto, (;))
_panel_spec(p::Pair{<:BMO.Detector, Symbol}) = (p.first, p.second, (;))
_panel_spec(p::Pair{<:BMO.Detector, <:Tuple{Symbol, NamedTuple}}) = (p.first, p.second...)
function _panel_spec(x)
    throw(ArgumentError("invalid detector panel $x, use `pd`, `pd => mode` or `pd => (mode, kwargs)`"))
end

"""Returns the `(detector, mode, kwargs)` of all panels of the `detectors` kwarg of `live_view`."""
function _panel_specs(detectors, systems)
    if detectors isa Symbol
        detectors == :auto || throw(ArgumentError("detectors must be :auto or a vector, got :$detectors"))
        return [(pd, :auto, (;)) for pd in _find_detectors(systems)]
    end
    specs = [_panel_spec(d) for d in detectors]
    for (_, mode, _) in specs
        mode in _PANEL_MODES || throw(ArgumentError("detector panel mode must be one of $_PANEL_MODES, got :$mode"))
    end
    return specs
end

# Options of a detector panel, which are not passed to `intensity`
const _PANEL_OPTIONS = (:colorscale, :colorrange, :history, :profiles)

"""
    _panel_options(kwargs)

Splits the `kwargs` of a detector panel into the panel options `(; colorscale, colorrange,
history, profiles)` and the remaining kwargs, which are passed to `intensity`.
"""
function _panel_options(kwargs::NamedTuple)
    colorscale = get(kwargs, :colorscale, :linear)
    colorscale in (:linear, :log) ||
        throw(ArgumentError("colorscale of a detector panel must be :linear or :log, got $(repr(colorscale))"))
    colorrange = get(kwargs, :colorrange, nothing)
    opts = (; colorscale, colorrange, history = Bool(get(kwargs, :history, false)),
        profiles = Bool(get(kwargs, :profiles, false)))
    rest = NamedTuple(k => v for (k, v) in pairs(kwargs) if !(k in _PANEL_OPTIONS))
    return opts, rest
end

function DetectorPanel(parent, pd::BMO.Detector, name::String, mode::Symbol, kwargs::NamedTuple)
    opts, kwargs = _panel_options(kwargs)
    grid = GridLayout(parent)
    # The metrics are shown in the subtitle, right above the axis, since the axis does not fill its
    # cell with `DataAspect`
    ax = Axis(grid[1, 1]; title = "$name: no hits", xlabel = "x [mm]", ylabel = "z [mm]",
        aspect = DataAspect(), subtitlesize = 11, subtitlecolor = :gray25)
    n = get(kwargs, :n, 100)
    xy = Observable(Point2f[])
    heat_x = Observable(zeros(Float32, n))
    heat_y = Observable(zeros(Float32, n))
    heat_I = Observable(zeros(Float32, n, n))
    heat_plot = heatmap!(ax, heat_x, heat_y, heat_I; colorrange = (0.0f0, 1.0f0), visible = false)
    scatter_plot = scatter!(ax, xy; markersize = 3, color = :black, visible = false)
    centroid = Observable(Point2f[])
    centroid_plot = scatter!(ax, centroid; marker = :cross, markersize = 12, color = :red)
    history_value, history_cx, history_cz = Observable(Point2f[]), Observable(Point2f[]), Observable(Point2f[])
    history_axes = Axis[]
    row = 2
    if opts.history
        value_ax = Axis(grid[row, 1]; height = 80, xlabel = "update", ylabel = "",
            xlabelsize = 11, ylabelsize = 11, xticklabelsize = 10, yticklabelsize = 10)
        # Centroid on a second y-axis on the right, i.e. a twin axis
        centroid_ax = Axis(grid[row, 1]; height = 80, yaxisposition = :right, ylabel = "c [mm]",
            ylabelsize = 11, yticklabelsize = 10, backgroundcolor = :transparent)
        hidexdecorations!(centroid_ax)
        hidespines!(centroid_ax)
        linkxaxes!(value_ax, centroid_ax)
        lines!(value_ax, history_value; color = :black)
        lines!(centroid_ax, history_cx; color = :red, linewidth = 1)
        lines!(centroid_ax, history_cz; color = :blue, linewidth = 1)
        push!(history_axes, value_ax, centroid_ax)
        row += 1
    end
    profile_x, profile_z = Observable(Point2f[]), Observable(Point2f[])
    profiles_ax = nothing
    if opts.profiles
        # x profile red, z profile blue, like the centroid of the history
        profiles_ax = Axis(grid[row, 1]; height = 80, xlabel = "x (red), z (blue) [mm]",
            ylabel = "I [W/m²]", xlabelsize = 11, ylabelsize = 11, xticklabelsize = 10,
            yticklabelsize = 10)
        lines!(profiles_ax, profile_x; color = :red)
        lines!(profiles_ax, profile_z; color = :blue)
    end
    rowgap!(grid, 4)
    return DetectorPanel(pd, name, mode, kwargs, opts.colorscale, opts.colorrange, ax, xy, heat_x,
        heat_y, heat_I, scatter_plot, heat_plot, nothing, centroid, centroid_plot,
        history_axes, 0, history_value, history_cx, history_cz, profiles_ax, profile_x, profile_z,
        nothing)
end

"""Formats `x` with 3 decimal places."""
function _fmt3(x::Real)
    s = string(round(x, digits = 3))
    occursin(r"[eE]", s) && return s
    i = findfirst('.', s)
    isnothing(i) && return s * ".000"
    return s * "0"^max(0, 3 - (length(s) - i))
end

_is_beamlet_hits(h) = h isa AbstractVector{<:BMO.AbstractBeamletHit}

function _resolve_mode(mode::Symbol, h)
    mode == :auto || return mode
    return _is_beamlet_hits(h) ? :intensity : :spot
end

function _clear_panel!(p::DetectorPanel)
    p.scatter_plot.visible[] && (p.scatter_plot.visible[] = false)
    p.heat_plot.visible[] && (p.heat_plot.visible[] = false)
    isempty(p.xy[]) || (p.xy[] = Point2f[])
    isempty(p.centroid[]) || (p.centroid[] = Point2f[])
    _clear_profiles!(p)
    p.metrics = nothing
    p.ax.subtitle[] = ""
    return nothing
end

function _clear_profiles!(p::DetectorPanel)
    isempty(p.profile_x[]) || (p.profile_x[] = Point2f[])
    isempty(p.profile_z[]) || (p.profile_z[] = Point2f[])
    return nothing
end

function _hits_title(p::DetectorPanel, h)
    n = length(h)
    return _is_beamlet_hits(h) ? "$(p.name): $n beamlets" : "$(p.name): $n rays"
end

"""
    _spot_metrics(pts)

Returns the metrics `(; n, cx, cz, rms, rmax)` of the spot diagram `pts` (local x/z [m]): the
number of hits, the centroid, the RMS radius `sqrt(mean(|p - c|²))` and the geometric radius, i.e.
the largest distance from the centroid [m].
"""
function _spot_metrics(pts)
    n = length(pts)
    cx = sum(q -> Float64(q[1]), pts) / n
    cz = sum(q -> Float64(q[2]), pts) / n
    r2 = [(Float64(q[1]) - cx)^2 + (Float64(q[2]) - cz)^2 for q in pts]
    return (; n, cx, cz, rms = sqrt(sum(r2) / n), rmax = sqrt(maximum(r2)))
end

"""
    _intensity_metrics(x, z, I)

Returns the metrics `(; P, cx, cz, wx, wz, peak)` of the intensity `I[i, j]` at `x[i]`, `z[j]`
[m]: the power (integrated like `optical_power`), the centroid, the 1/e² radii `2σ` from the second
moments along x and z and the peak intensity. The centroid and the radii are `NaN` without
intensity.
"""
function _intensity_metrics(x, z, I)
    P = BMO.trapz((x, z), I)
    peak = maximum(I)
    S = sum(I)
    S > 0 || return (; P, cx = NaN, cz = NaN, wx = NaN, wz = NaN, peak)
    # Marginal distributions along x and z, the uniform cell area cancels
    Ix, Iz = vec(sum(I; dims = 2)), vec(sum(I; dims = 1))
    cx, cz = dot(Ix, x) / S, dot(Iz, z) / S
    σx2 = sum(Ix .* (x .- cx) .^ 2) / S
    σz2 = sum(Iz .* (z .- cz) .^ 2) / S
    return (; P, cx, cz, wx = 2 * sqrt(σx2), wz = 2 * sqrt(σz2), peak)
end

"""Formats the length `x` like `_length_string` with its sign, below 1 pm as `0 nm`."""
function _signed_length_string(x)
    abs(x) < 1e-12 && return "0 nm"
    return x < 0 ? "-" * _length_string(-x) : _length_string(x)
end

"""Formats the metrics `m` of a panel, see `_spot_metrics` and `_intensity_metrics`, in two lines."""
function _metrics_string(m)
    if haskey(m, :rms)
        return "N = $(m.n), c = ($(_signed_length_string(m.cx)), $(_signed_length_string(m.cz)))\n" *
               "rms $(_length_string(m.rms)), max $(_length_string(m.rmax))"
    end
    s = "P = $(_fmt3(1e3 * m.P)) mW, peak $(_fmt_sigdigits(m.peak)) W/m²"
    isfinite(m.cx) || return s * "\nno intensity"
    return s * "\nc = ($(_signed_length_string(m.cx)), $(_signed_length_string(m.cz))), " *
           "w = ($(_length_string(m.wx)), $(_length_string(m.wz)))"
end

"""
    _set_metrics!(p, m; record = true)

Shows the metrics `m` of the panel `p` in the subtitle of its axis and marks the centroid. With
`record`, the power (or the number of hits) and the centroid are added to the history of the
panel, if any.
"""
function _set_metrics!(p::DetectorPanel, m; record::Bool = true)
    p.metrics = m
    p.ax.subtitle[] = _metrics_string(m)
    c = isfinite(m.cx) ? [Point2f(1e3 * m.cx, 1e3 * m.cz)] : Point2f[]
    p.centroid[] == c || (p.centroid[] = c)
    record && _record_history!(p, m)
    return nothing
end

"""
    _record_history!(p, m; draw = true)

Adds the power (or the number of hits) and the centroid of the metrics `m` to the history of the
panel `p`, if it has one. Without `draw`, only the data is extended and the plots are not updated,
e.g. for a panel that is not shown, see `_draw_history!`.
"""
function _record_history!(p::DetectorPanel, m; draw::Bool = true)
    isempty(p.history_axes) && return nothing
    p.history_count += 1
    k = p.history_count
    value = haskey(m, :P) ? 1e3 * m.P : m.n
    for (obs, v) in ((p.history_value, value), (p.history_cx, 1e3 * m.cx), (p.history_cz, 1e3 * m.cz))
        push!(obs[], Point2f(k, v))
        length(obs[]) > _HISTORY_LENGTH && popfirst!(obs[])
    end
    draw && _draw_history!(p, m)
    return nothing
end

"""
Updates the plots of the history of the panel `p` from its data (see `_record_history!`), the
label of the power or number of hits from the metrics `m` of the last record.
"""
function _draw_history!(p::DetectorPanel, m = p.metrics)
    isempty(p.history_axes) && return nothing
    foreach(notify, (p.history_value, p.history_cx, p.history_cz))
    label = (isnothing(m) || haskey(m, :P)) ? "P [mW]" : "N"
    p.history_axes[1].ylabel[] == label || (p.history_axes[1].ylabel[] = label)
    foreach(autolimits!, p.history_axes)
    return nothing
end

function _update_spot!(p::DetectorPanel, h; preview = false, record = true)
    pts = BMO.spot_diagram(p.pd)
    p.xy[] = [Point2f(1e3 * q[1], 1e3 * q[2]) for q in pts]
    p.heat_plot.visible[] && (p.heat_plot.visible[] = false)
    p.scatter_plot.visible[] || (p.scatter_plot.visible[] = true)
    p.ax.title[] = _hits_title(p, h) * (preview ? " (preview)" : "")
    _set_metrics!(p, _spot_metrics(pts); record)
    _clear_profiles!(p)
    autolimits!(p.ax)
    _pad_degenerate_limits!(p.ax, p.xy[])
    return nothing
end

"""
    _pad_degenerate_limits!(ax, xy; ε = 1e-6)

Sets the limits of the spot diagram `ax` if the x or z extent of the spots `xy` [mm] is below `ε`,
e.g. for a single ray, since `DataAspect` gives degenerate limits for a zero-width data range. A
degenerate axis gets half the extent of the other axis around its center, or 1 µm if both are
degenerate. Otherwise the limits of `autolimits!` are kept.
"""
function _pad_degenerate_limits!(ax::Axis, xy; ε = 1e-6)
    isempty(xy) && return nothing
    xmin, xmax = extrema(q -> Float64(q[1]), xy)
    zmin, zmax = extrema(q -> Float64(q[2]), xy)
    Δx, Δz = xmax - xmin, zmax - zmin
    (Δx < ε || Δz < ε) || return nothing
    h = max(Δx, Δz) / 2
    h = h < ε ? 1e-3 : h
    hx = Δx < ε ? h : Δx / 2 * 1.05
    hz = Δz < ε ? h : Δz / 2 * 1.05
    cx, cz = (xmin + xmax) / 2, (zmin + zmax) / 2
    limits!(ax, cx - hx, cx + hx, cz - hz, cz + hz)
    return nothing
end

"""Grid size of the intensity of the panel `p`, reduced to a preview while moving objects."""
function _panel_n(p::DetectorPanel, coarse::Bool)
    n = get(p.kwargs, :n, 100)
    return coarse ? max(16, n ÷ 4) : n
end

"""
    _display_intensity(I, colorscale) -> (values, colorrange)

Returns the values of the heatmap of the intensity `I` and their color range: `I` for the linear
scale, `log10(I)` with a floor of `_LOG_FLOOR` times the maximum for the logarithmic scale.
"""
function _display_intensity(I, colorscale::Symbol)
    Imax = Float64(maximum(I))
    colorscale == :linear && return I, (0.0, Imax > 0 ? Imax : 1.0)
    floor = Imax > 0 ? _LOG_FLOOR * Imax : _LOG_FLOOR
    return log10.(max.(I, floor)), (log10(floor), log10(floor / _LOG_FLOOR))
end

"""Shows the intensity `I` along x and z through the centroid of the metrics `m` of the panel `p`."""
function _update_profiles!(p::DetectorPanel, x, z, I, m)
    isnothing(p.profiles_ax) && return nothing
    if !isfinite(m.cx)
        _clear_profiles!(p)
        return nothing
    end
    i, j = argmin(abs.(x .- m.cx)), argmin(abs.(z .- m.cz))
    p.profile_x[] = [Point2f(1e3 * x[k], I[k, j]) for k in eachindex(x)]
    p.profile_z[] = [Point2f(1e3 * z[k], I[i, k]) for k in eachindex(z)]
    autolimits!(p.profiles_ax)
    return nothing
end

"""Marks a panel whose hits are shown as a spot diagram, see `_panel_field`."""
struct _SpotField end

"""
    _panel_field(p::DetectorPanel, coarse::Bool)

Computes what the panel `p` shows after a solve, without changing any plot, such that it can run
in a background task (see `_solve!`): `nothing` without hits, `_SpotField()` for a spot diagram,
the intensity `(x, z, I)` on the grid of the panel (coarse: see `_panel_n`), or the exception of a
failed computation, which `_update_panel!` reports. A cancelled solve is rethrown, see
`BMO._is_cancelled`.
"""
function _panel_field(p::DetectorPanel, coarse::Bool)
    try
        h = BMO.hits(p.pd)
        isnothing(h) && return nothing
        _resolve_mode(p.mode, h) == :intensity || return _SpotField()
        return _intensity_field(p, h, coarse)
    catch e
        BMO._is_cancelled(e) && rethrow()
        return e
    end
end

_intensity_field(p::DetectorPanel, _, coarse::Bool) =
    BMO.intensity(p.pd; merge(p.kwargs, (; n = _panel_n(p, coarse)))...)

# In the auto mode, the spot diagram is shown if the field of the hits can not be evaluated
function _intensity_field(p::DetectorPanel, h::AbstractVector{<:BMO.AstigmaticGaussianBeamletHit},
        coarse::Bool)
    try
        return invoke(_intensity_field, Tuple{DetectorPanel, Any, Bool}, p, h, coarse)
    catch e
        (BMO._is_cancelled(e) || p.mode != :auto) && rethrow()
        return _SpotField()
    end
end

function _update_intensity!(p::DetectorPanel, h, (x, z, I); coarse = false, preview = false,
        record = true)
    # The plot is updated lazily, hence the grid size may change
    p.heat_x[] = Float32.(1e3 .* x)
    p.heat_y[] = Float32.(1e3 .* z)
    values, range = _display_intensity(I, p.colorscale)
    p.heat_I[] = Float32.(values)
    p.heat_plot.colorrange[] = Float32.(something(p.colorrange, range))
    p.scatter_plot.visible[] && (p.scatter_plot.visible[] = false)
    p.heat_plot.visible[] || (p.heat_plot.visible[] = true)
    # Optical power from the computed intensity like optical_power, avoids a second field evaluation
    m = _intensity_metrics(x, z, I)
    if h isa AbstractVector{<:BMO.GaussianBeamletHit}
        p.ax.title[] = "$(p.name): P = $(_fmt3(1e3 * m.P)) mW"
    else
        p.ax.title[] = _hits_title(p, h)
    end
    (coarse || preview) && (p.ax.title[] *= " (preview)")
    _set_metrics!(p, m; record)
    _update_profiles!(p, x, z, I, m)
    autolimits!(p.ax)
    return nothing
end

"""
    _update_panel!(p[, field]; coarse = false, preview = false, record = true)

Updates the plots, the title and the metrics of the panel `p` after the systems have been solved,
with the `field` computed by `_panel_field` (computed here if not given). `coarse` computes the
intensity on a coarse grid, `preview` marks the title after a preview solve, see `_resolve!`.
`record` adds the metrics to the history of the panel.
"""
function _update_panel!(p::DetectorPanel, field; coarse = false, preview = false, record = true)
    try
        _show_field!(p, field, BMO.hits(p.pd); coarse, preview, record)
        p.last_error = nothing
    catch e
        p.last_error = _log_once(e, p.last_error, "update of the panel \"$(p.name)\"")
        p.ax.title[] = "$(p.name): error"
    end
    return nothing
end

_update_panel!(p::DetectorPanel; coarse = false, kwargs...) =
    _update_panel!(p, _panel_field(p, coarse); coarse, kwargs...)

function _show_field!(p::DetectorPanel, ::Nothing, _; preview, _...)
    _clear_panel!(p)
    p.ax.title[] = "$(p.name): no hits" * (preview ? " (preview)" : "")
    return nothing
end

_show_field!(p::DetectorPanel, ::_SpotField, h; preview, record, _...) =
    _update_spot!(p, h; preview, record)
_show_field!(p::DetectorPanel, field::Tuple, h; kwargs...) = _update_intensity!(p, h, field; kwargs...)
_show_field!(::DetectorPanel, e::Exception, _; _...) = throw(e)

#=
Panels shown by the layout: the compact layout shows, and thus computes, all panels after each
solve; the app layout only the panel of its active tab, see `LiveDock.jl`
=#

"""
    _computed_panels(gui, preview::Bool) -> Vector

The detector panels of the `gui` whose fields are computed after a solve (a preview solve if
`preview`), see `_compute`: all panels, unless the layout shows only some of them. Their results
are shown by `_apply_panel!`.
"""
_computed_panels(gui::LiveView, ::Bool) = gui.panels

"""
    _shown_panels(gui) -> Vector

The detector panels of the `gui` that are shown, i.e. refined after a coarse preview, see
`_on_idle!`: all panels, unless the layout shows only some of them.
"""
_shown_panels(gui::LiveView) = gui.panels

"""
    _apply_panel!(gui, p, field; coarse, preview)

Shows the `field` of the panel `p` (see `_panel_field`) after a solve of the `gui`, the metrics of a
full solve are recorded in its history, see `_update_panel!`.
"""
_apply_panel!(::LiveView, p::DetectorPanel, field; coarse, preview) =
    _update_panel!(p, field; coarse, preview, record = !preview)

# Hooks of the layout around a solve: before a solve empties the detectors, and after its result
# was shown, see `_start_job` and `_apply!`
_on_solve_started!(::LiveView) = nothing
_on_applied!(::LiveView) = nothing

#=
Detectors: the options of their panel, see `_set_panel_options!`, on their cards (`_panel_row`,
see `card_rows` in LiveCardRows.jl)
=#

"""
The row of the options of the panel of a detector on its card: a button that cycles the mode (see
`_PANEL_MODES`) and shows the current one, and a toggle of the logarithmic color scale. The row is
the same for all detectors; for a detector without a panel (see the `detectors` kwarg of
`live_view`), the button reads "no panel" and the inputs only show a message.
"""
_panel_row() = CardRow("panel",
    CardWidget(Button; name = :panel_mode, label = "auto",
        value = (gui, pd) -> _panel_mode_label(_panel_of(gui, pd)),
        on = (gui, pd, _) -> _cycle_panel_mode!(gui, _panel_of(gui, pd))),
    CardWidget(Toggle; name = :panel_log, value = (gui, pd) -> _is_log(_panel_of(gui, pd)),
        on = (gui, pd, v) -> _set_panel_log!(gui, _panel_of(gui, pd), v)),
    "log")

# The modes of a detector panel, in the order of the button of `_panel_row`
const _PANEL_MODES = (:auto, :spot, :intensity)

_panel_mode_label(::Nothing) = "no panel"
_panel_mode_label(p::DetectorPanel) = string(p.mode)
_is_log(::Nothing) = false
_is_log(p::DetectorPanel) = p.colorscale == :log

function _cycle_panel_mode!(gui::LiveView, p::DetectorPanel)
    i = something(findfirst(==(p.mode), _PANEL_MODES), 0)
    _set_panel_options!(gui, p; mode = _PANEL_MODES[mod1(i + 1, length(_PANEL_MODES))])
    gui.status.text[] = "panel $(p.name): $(p.mode)"
    return nothing
end
_cycle_panel_mode!(gui::LiveView, ::Nothing) = _no_panel(gui)

_set_panel_log!(gui::LiveView, p::DetectorPanel, log::Bool) =
    _set_panel_options!(gui, p; colorscale = log ? :log : :linear)
_set_panel_log!(gui::LiveView, ::Nothing, _) = _no_panel(gui)

_no_panel(gui::LiveView) = (gui.status.text[] = "this detector has no panel, see the detectors kwarg"; nothing)

"""Returns the first detector panel of `pd` in the `gui`, or `nothing`."""
_panel_of(gui::LiveView, pd) = (i = findfirst(p -> p.pd === pd, gui.panels);
    isnothing(i) ? nothing : gui.panels[i])

"""
    _set_panel_options!(gui, p; mode = p.mode, colorscale = p.colorscale)

Sets the mode (`:auto`, `:spot` or `:intensity`) and the color scale (`:linear` or `:log`) of the
detector panel `p` of the `gui` and shows its current hits again, see `_refresh_panel!`.
"""
function _set_panel_options!(gui::LiveView, p::DetectorPanel; mode::Symbol = p.mode,
        colorscale::Symbol = p.colorscale)
    (p.mode == mode && p.colorscale == colorscale) && return nothing
    p.mode = mode
    p.colorscale = colorscale
    _refresh_panel!(gui, p)
    return nothing
end

"""
Shows the current hits of the detector panel `p` again, e.g. after its options changed, without
solving and without recording its history.
"""
_refresh_panel!(::LiveView, p::DetectorPanel) = _update_panel!(p; coarse = false, record = false)

