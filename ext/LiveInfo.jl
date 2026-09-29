#=
Info cards of the live view: the results of the beam inspection and of the measurement on a card
at their points in the 3D view, built like the cards of the components (see `_ComponentCard`)
from the declarations `card_rows` of the types below
=#

"""
    _InfoItem

A result of the live view that is shown on a card at its points in the 3D view instead of a
component: an inspected point of a beam (`_BeamPoint`) or a measurement (`_Measurement`). It has
the `points` [m] that the card is attached to, a `title`, the `kind` of its icon (see `_icon`) and
declares the rows of its card by [`card_rows`](@ref), without actions. The card is transient (see
`_show_info!`): the next inspection or measurement replaces it, a click elsewhere or `Esc` removes
it, unless it is pinned, which keeps it like the card of a component, docked in the app layout.
"""
abstract type _InfoItem end

"""The inspected point `info` of a beam, see `_inspect_beam`, and its marker once it is kept."""
struct _BeamPoint <: _InfoItem
    info::NamedTuple
    plots::Vector{AbstractPlot}
end
_BeamPoint(info::NamedTuple) = _BeamPoint(info, AbstractPlot[])

"""
The measurement between the measured `points` (`(; point, obj)`, one or two) with their `names`
and the `result` of `_measure` (`nothing` for one point), and its line and markers once it is kept.
"""
struct _Measurement <: _InfoItem
    points::Vector{Any}
    result::Any
    names::Vector{String}
    plots::Vector{AbstractPlot}
end
_Measurement(points, result, names) = _Measurement(points, result, names, AbstractPlot[])

_points(x::_BeamPoint) = [Point3f(x.info.point)]
_points(x::_Measurement) = [Point3f(m.point) for m in x.points]
_title(::_BeamPoint) = "Beam"
_title(::_Measurement) = "Measurement"

# The item takes the place of an object on a card: its label, icon, pose and bounding box
_label(::LiveView, x::_InfoItem) = _title(x)
_label(::AppView, x::_InfoItem) = _title(x)
_tree_kind(::_BeamPoint) = :trace
_tree_kind(::_Measurement) = :measure
_pose(x::_InfoItem) = (Point3{Float64}(sum(_points(x)) / length(_points(x))), Matrix{Float64}(I, 3, 3))
BMO.position(x::_InfoItem) = Vector{Float64}(_pose(x)[1])
_card_corners(::LiveView, ::_ComponentCard, x::_InfoItem) = _points(x)

card_actions(::_InfoItem) = ()

# An inspected point: position, direction, path and optical path length, and for Gaussian beamlets
# the radius and the radius of curvature, see `_inspection_string`
function card_rows(x::_BeamPoint)
    i = x.info
    rows = [_text_row("at", :at, (_, _) -> _point_string(i.point)),
        _text_row("dir", :dir, (_, _) -> _direction_string(i.direction)),
        _text_row("path", :path, (_, _) -> "$(_fmt3(1e3 * i.length)) mm"),
        _text_row("OPL", :opl, (_, _) -> "$(_fmt3(1e3 * i.opl)) mm")]
    if !isnothing(i.w)
        r = iszero(i.R) ? "∞" : _signed_length_string(1 / i.R)
        push!(rows, _text_row("w", :w, (_, _) -> _length_string(i.w)), _text_row("R", :R, (_, _) -> r))
    end
    return Tuple(rows)
end

# A measurement: the first point until the second one is clicked, then the distance, its components
# and the angle between the optical axes of two components
function card_rows(x::_Measurement)
    length(x.points) == 1 && return (
        _text_row("at", :at, (_, _) -> "$(x.names[1]), $(_point_string(x.points[1].point))"),
        _text_row("next", :next, (_, _) -> "click the second point"))
    a, b = x.points
    d = b.point .- a.point
    rows = [_text_row("from", :from, (_, _) -> x.names[1]), _text_row("to", :to, (_, _) -> x.names[2]),
        _text_row("d", :distance, (_, _) -> "$(round(1e3 * x.result.distance, digits = 6)) mm"),
        _text_row("Δ", :delta, (_, _) -> _point_string(d))]
    isnothing(x.result.angle) ||
        push!(rows, _text_row("angle", :angle, (_, _) -> _angle_string(x.result.angle)))
    return Tuple(rows)
end

"""Returns the transient card of the `gui` (see `_show_info!`), or `nothing`."""
function _info_card(gui::LiveView)
    i = findfirst(c -> c.transient, gui.cards.all)
    return isnothing(i) ? nothing : gui.cards.all[i]
end

"""
    _show_info!(gui, x::_InfoItem)

Shows `x` on the transient card of the `gui`, a floating card at the points of `x` that replaces the
card of the previous item: shown like a pinned card, but with its pin off. Its pin keeps the card,
see `_keep_info!`; `_release_info!` removes it.
"""
function _show_info!(gui::LiveView, x::_InfoItem)
    c = something(_info_card(gui), _spare_card!(gui))
    c.pinned, c.transient, c.obj, c.key, c.pose, c.collapsed = true, true, x, nothing, nothing, false
    _show_head!(c)
    _update_cards!(gui)
    return nothing
end

"""Removes the transient card of the `gui` if it shows an item of the type `T`, see `_show_info!`."""
function _release_info!(gui::LiveView, ::Type{T}) where {T}
    c = _info_card(gui)
    (isnothing(c) || !(c.obj isa T)) && return nothing
    c.pinned, c.transient, c.obj = false, false, nothing
    _hide_card!(c)
    _show_head!(c)
    return nothing
end

"""
    _keep_info!(gui, c)

Keeps the transient card `c` of the `gui` when its pin is clicked: it becomes a pinned card, which
the next inspection or measurement does not replace. In the app layout, whose pinned cards are
docked (see `_pin!`), the item is pinned there and the floating card is removed.
"""
function _keep_info!(gui::LiveView, c::_ComponentCard)
    _keep_plots!(gui, c.obj)
    c.transient = false
    _show_head!(c)
    _update_cards!(gui)
    _on_pinned!(gui)
    return nothing
end
function _keep_info!(gui::AppView, c::_ComponentCard)
    x = c.obj
    _keep_plots!(gui, x)
    _release_info!(gui, _InfoItem)
    _pin!(gui, x)
    _on_pinned!(gui)
    return nothing
end

"""
    _keep_plots!(gui, x::_InfoItem)

Takes the plots of the current inspection or measurement of the `gui` over to the kept item `x`,
such that the next inspection or measurement does not remove them; `_forget!` removes them with the
card of `x`.
"""
function _keep_plots!(gui::LiveView, x::_BeamPoint)
    isnothing(gui.measure.inspection_plot) && return nothing
    push!(x.plots, gui.measure.inspection_plot)
    gui.measure.inspection_plot = nothing
    return nothing
end
function _keep_plots!(gui::LiveView, x::_Measurement)
    append!(x.plots, gui.measure.plots)
    empty!(gui.measure.plots)
    return nothing
end

"""Removes the plots of the kept item `x` when its card is unpinned, see `_keep_plots!`."""
function _forget!(gui::LiveView, x::_InfoItem)
    foreach(p -> delete!(gui.ax, p), x.plots)
    empty!(x.plots)
    return nothing
end
_forget!(::LiveView, _) = nothing

# The dock button of the transient card of `x` in the app layout keeps it, docked like its pin does,
# see `_keep_info!`; a kept card moves like the card of a component, see `_dock!`
function _dock!(gui::AppView, x::_InfoItem)
    c = _info_card(gui)
    !isnothing(c) && c.obj === x && return _keep_info!(gui, c)
    return invoke(_dock!, Tuple{AppView, Any}, gui, x)
end
