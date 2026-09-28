#=
Property inspector of the app layout of the live view (`live_view(...; layout = :app)`): the
"PROPERTIES" section of the right sidebar, see `_Inspector`
=#

using Makie: Button, Box, Label, Textbox, GridLayout, Observable, Point2f, RGBAf, BezierPath,
             scatter!, text!, linesegments!, lift, on, rowgap!, colgap!, rowsize!, Fixed, Auto

#=
Property list
=#

# Geometry of the property list in pixels
const _PROPERTY_ROW = 20.0f0       # row height
const _PROPERTY_FONTSIZE = 12.0f0
const _PROPERTY_MAX_ROWS = 14      # longer lists end with a row "… n more"
const _PROPERTY_GAP = 12.0f0       # minimal gap between a label and its value

"""
    _PropertyList(parent; label_color, value_color, line_color)

A list of `label  value` rows at the grid position `parent`, the label left in `label_color`, the
value right-aligned in `value_color`, rows separated by thin lines. Drawn with a constant number of
plots (one `text!` for all labels, one for all values, one `linesegments!` for the lines) in the
scene of an invisible `Box` that holds the place in the layout; the box is as high as the rows.
Values that do not fit are ellipsized, lists with more than `_PROPERTY_MAX_ROWS` rows end with
a row "… n more". Updated only by [`_set_rows!`](@ref) and on layout changes.
"""
mutable struct _PropertyList
    const box::Box
    const labels::Makie.AbstractPlot
    const values::Makie.AbstractPlot
    const lines::Makie.AbstractPlot
    rows::Vector{Tuple{String, String}}
    const font::Any
    # widths of texts in pixels, see `_text_width`
    const widths::Dict{String, Float32}
end

function _PropertyList(parent; label_color, value_color, line_color)
    box = Box(parent; visible = false, height = 0, tellwidth = false)
    scene = box.blockscene
    common = (; space = :pixel, inspectable = false, fontsize = _PROPERTY_FONTSIZE,
        markerspace = :pixel)
    labels = text!(scene, Point2f[]; common..., text = String[], color = label_color,
        align = (:left, :center))
    values = text!(scene, Point2f[]; common..., text = String[], color = value_color,
        align = (:right, :center))
    lines = linesegments!(scene, Point2f[]; space = :pixel, inspectable = false,
        color = line_color, linewidth = 1)
    list = _PropertyList(box, labels, values, lines, Tuple{String, String}[],
        _tree_font(scene, :regular), Dict{String, Float32}())
    on(_ -> _redraw!(list), scene, box.layoutobservables.computedbbox)
    return list
end

"""
    _set_rows!(list::_PropertyList, rows)

Shows the `rows` (`(label, value)` tuples of `String`s) in the `list`. The height of the list
changes only if the number of rows does, i.e. the layout is only updated then.
"""
function _set_rows!(list::_PropertyList, rows::Vector{Tuple{String, String}})
    if length(rows) > _PROPERTY_MAX_ROWS
        n = length(rows) - _PROPERTY_MAX_ROWS + 1
        rows = [rows[1:(_PROPERTY_MAX_ROWS - 1)]; ("… $n more", "")]
    end
    list.rows = rows
    height = length(rows) * _PROPERTY_ROW
    if list.box.height[] != height
        # the layout calls `_redraw!` via the new bounding box
        list.box.height[] = height
    else
        _redraw!(list)
    end
    return list
end

function _redraw!(list::_PropertyList)
    bb = list.box.layoutobservables.computedbbox[]
    (x0, y0), (w, h) = Makie.origin(bb), Makie.widths(bb)
    n = length(list.rows)
    label_pos, value_pos = Vector{Point2f}(undef, n), Vector{Point2f}(undef, n)
    label_text, value_text = Vector{String}(undef, n), Vector{String}(undef, n)
    lines = Point2f[]
    fs = _PROPERTY_FONTSIZE
    for (i, (label, value)) in enumerate(list.rows)
        y = y0 + h - (i - 0.5f0) * _PROPERTY_ROW
        # labels are short, they may take up to half the width
        label_text[i] = _fit_text(list.widths, list.font, fs, label, w / 2)
        room = w - _text_width(list.widths, list.font, fs, label_text[i]) - _PROPERTY_GAP
        value_text[i] = _fit_text(list.widths, list.font, fs, value, room)
        label_pos[i] = Point2f(x0, y)
        value_pos[i] = Point2f(x0 + w, y)
        # a line below each row but the last, on a pixel center
        if i < n
            yl = round(y - _PROPERTY_ROW / 2) + 0.5f0
            push!(lines, Point2f(x0, yl), Point2f(x0 + w, yl))
        end
    end
    Makie.update!(list.labels; arg1 = label_pos, text = label_text)
    Makie.update!(list.values; arg1 = value_pos, text = value_text)
    Makie.update!(list.lines; arg1 = lines)
    return nothing
end

#=
Formatting of `properties`: the unit in brackets at the end of a name selects the formatting of
its value, see `BMO.properties`
=#

"""
    _property_row(name, value) -> (label, text)

Returns the label (the `name` without its unit) and the formatted `value` of a property, see
`BeamletOptics.properties`: lengths `[m]` in nm, µm, mm or m, angles `[rad]` in µrad, mrad or °,
vectors in parentheses, other numbers with 4 significant digits, `Bool`s as yes or no.
"""
function _property_row(name::AbstractString, value)
    m = match(r"^(.*\S)\s*\[([^\]]+)\]$", name)
    isnothing(m) && return (String(name), _format_value(value, Val(nothing)))
    label, unit = m.captures
    return (String(label), _format_value(value, Val(Symbol(unit))))
end

# By the type of the value, then numbers and vectors by the unit
_format_value(x, ::Val) = string(x)
_format_value(s::AbstractString, ::Val) = String(s)
_format_value(x::Bool, ::Val) = x ? "yes" : "no"
_format_value(x::Real, unit::Val) = _format_number(x, unit)
_format_value(v::AbstractVector{<:Real}, unit::Val) = _format_vector(v, unit)

_format_number(x::Integer, ::Val{nothing}) = string(x)
_format_number(x, ::Val{nothing}) = _fmt_digits(x)
_format_number(x, ::Val{U}) where {U} = "$(_fmt_digits(x)) $U"
_format_number(x, ::Val{:m}) = _signed(x, _length_units)
_format_number(x, ::Val{:rad}) = _signed(x, _angle_units)

_vector_string(v, f) = "(" * join((f(x) for x in v), ", ") * ")"
# directions and other dimensionless vectors with 3 decimals
_format_vector(v, ::Val{nothing}) = _vector_string(v, x -> _fmt_digits(round(x, digits = 3)))
_format_vector(v, ::Val{U}) where {U} = _vector_string(v, _fmt_digits) * " $U"
function _format_vector(v, ::Val{:m})
    # a common unit, chosen by the largest component
    m = maximum(abs, v; init = 0.0)
    i = findfirst(((f, _),) -> round(m / f, sigdigits = 3) < 1000, _length_units)
    factor, unit = _length_units[something(i, length(_length_units))]
    return _vector_string(v, x -> _fmt_digits(x / factor, 3)) * " $unit"
end

const _length_units = (1e-9 => "nm", 1e-6 => "µm", 1e-3 => "mm", 1.0 => "m")
const _angle_units = (1e-6 => "µrad", 1e-3 => "mrad", deg2rad(1) => "°")

"""Formats `x` like `_unit_string` with its sign, `0` without a unit."""
function _signed(x::Real, units)
    iszero(x) && return "0"
    s = _unit_string(abs(x), units)
    return x < 0 ? "-" * s : s
end

"""Formats `x` with `digits` significant digits, integers without `.0`, `-0` as `0`."""
function _fmt_digits(x::Real, digits::Int = 4)
    v = round(Float64(x), sigdigits = digits)
    iszero(v) && return "0"
    return isinteger(v) && abs(v) < 1e15 ? string(Int(v)) : string(v)
end

#=
Segmented control
=#

"""
    _Segmented(parent, options; theme)

A segmented control at the grid position `parent`: a row of text buttons, one per option
`key => label`, of which the one of `selected[]` (a key) is highlighted with the accent colors of
the `theme` tokens. A click sets `selected`, which can also be set from code.
"""
struct _Segmented
    grid::GridLayout
    buttons::Vector{Button}
    keys::Vector{Symbol}
    selected::Observable{Symbol}
end

function _Segmented(parent, options::Vector{Pair{Symbol, String}}; theme, selected::Symbol = first(options).first)
    t = theme
    grid = GridLayout(parent; default_colgap = 2, halign = :left, tellwidth = false)
    sel = Observable(selected)
    ks = first.(options)
    buttons = [Button(grid[1, i]; label, fontsize = 12, padding = (7, 7, 4, 4), cornerradius = 4)
               for (i, (_, label)) in enumerate(options)]
    function look!(k)
        for (key, b) in zip(ks, buttons)
            on_ = key == k
            bg, fg = on_ ? (t.accent_soft, t.accent) : (t.field, t.text)
            b.buttoncolor[] == bg || (b.buttoncolor[] = bg)
            b.labelcolor[] == fg || (b.labelcolor[] = fg)
            b.labelcolor_hover[] = fg
            b.labelcolor_active[] = fg
            b.strokecolor[] = on_ ? t.accent : t.border
        end
        return nothing
    end
    on(look!, sel; update = true)
    for (key, b) in zip(ks, buttons)
        on(_ -> (sel[] == key || (sel[] = key)), b.clicks)
    end
    return _Segmented(grid, buttons, ks, sel)
end

#=
Inspector
=#

"""
    _Inspector

The property inspector of the app layout, the "PROPERTIES" section of the right sidebar, from top
to bottom:

- the header: the icon of the kind of the selected object (see `_tree_kind`), its name (see
  `_label`) and its type
- the pose: the pose boxes (see `_update_pose_boxes!`), the step box and the mode as a segmented
  control (`mode`), bound to the `mode` of the controls
- the properties of the object, see `BeamletOptics.properties` and `_PropertyList`; without a
  selection, a summary of the live view
- the type-dependent sections in `context`, see `_inspector_sections`

The inspector is updated on events only (see `_refresh_inspector!`): the selection, moves and
solves. The type-dependent sections are only rebuilt if the selected object needs other sections,
e.g. when a detector follows a lens; otherwise their values are updated in place.
"""
mutable struct _Inspector
    const grid::GridLayout
    const icon::Observable{BezierPath}
    const icon_color::Observable{RGBAf}
    const name::Label
    const type::Label
    const mode::_Segmented
    const list::_PropertyList
    # type-dependent sections: the registry entries they were built from, their layout and their
    # update functions `update!(obj)`
    context::GridLayout
    sections::Vector{Pair{Type, Function}}
    updates::Vector{Function}
    # the object shown (`nothing`: the summary), a flag that nothing was shown yet
    shown::Any
    fresh::Bool
    # widths of the texts of the header in pixels, per label, see `_text_width`
    const widths::NTuple{2, Dict{String, Float32}}
end

# Names of `properties` that the inspector shows elsewhere: in the header and the pose boxes
const _INSPECTOR_SKIPPED = ("Type", "Position [m]")

"""
    _build_inspector!(layout::AppLayout) -> (; pose_boxes, step_box)

Creates the "PROPERTIES" section of the right sidebar, see `_Inspector`, and registers the
built-in type-dependent sections, see `_inspector_sections`. Returns the widgets that are fields of
`LiveView`.
"""
function _build_inspector!(layout::AppLayout)
    t = layout.theme
    g = _add_sidebar_section!(layout, :right, "Properties")
    # Header: icon, name and type
    header = GridLayout(g[1, 1]; default_colgap = 8, tellwidth = false)
    icon_box = Box(header[1:2, 1]; width = 24, height = 24, visible = false)
    icon = Observable(_icon(:system))
    icon_color = Observable(RGBAf(Makie.to_color(t.muted)))
    center = lift(r -> Point2f(Makie.origin(r) .+ Makie.widths(r) ./ 2), icon_box.blockscene,
        icon_box.layoutobservables.computedbbox)
    scatter!(icon_box.blockscene, center; marker = icon, markersize = 22, color = icon_color,
        markerspace = :pixel, inspectable = false)
    name = Label(header[1, 2], "No selection"; halign = :left, font = :bold, fontsize = 14,
        tellwidth = false)
    type = Label(header[2, 2], " "; halign = :left, color = t.muted, fontsize = 12,
        tellwidth = false)
    rowgap!(header, 0)
    # Pose boxes, step and mode
    pose = GridLayout(g[2, 1]; default_colgap = 6, default_rowgap = 2, tellwidth = false)
    Label(pose[1, 1], "step"; halign = :left, fontsize = 11, color = t.muted, tellwidth = false)
    Label(pose[1, 2:3], "mode"; halign = :left, fontsize = 11, color = t.muted, tellwidth = false)
    step_box = Textbox(pose[2, 1]; placeholder = "250 nm", width = 72, halign = :left)
    mode = _Segmented(pose[2, 2:3], [:move => "Move", :rotate => "Rotate"]; theme = t)
    # Properties and type-dependent sections
    Box(g[3, 1]; height = 1, color = t.border, strokewidth = 0)
    list = _PropertyList(g[4, 1]; label_color = t.muted, value_color = t.text,
        line_color = RGBAf(Makie.to_color(t.border)))
    context = GridLayout(g[5, 1]; tellwidth = false)
    # an empty layout has no size, which would let the sidebar squeeze the inspector
    rowsize!(g, 5, Fixed(0))
    rowgap!(g, 10)
    layout.inspector = _Inspector(g, icon, icon_color, name, type, mode, list, context,
        Pair{Type, Function}[], Function[], nothing, true, (Dict{String, Float32}(), Dict{String, Float32}()))
    layout.inspectors = Pair{Type, Function}[BMO.Detector => _detector_section]
    return (; step_box)
end

"""
    _inspector_sections(gui::AppView, obj) -> Vector{Pair{Type, Function}}

Returns the entries `T => build` of the registry `gui.layout.inspectors` that apply to `obj`, in
the order of the registry. `build(gui, grid)` creates a type-dependent section of the inspector in
the `GridLayout` `grid` and returns its update function `update!(obj)`, which sets the values of
the section for the selected `obj`; it is called after each selection, move and solve.

The entries are matched with `obj isa T` in a loop instead of by dispatch, since entries are added
at runtime, e.g. for types of the user, see `_add_inspector!`; the built-in sections (detector
panel settings) are registered the same way, in `_build_inspector!`. Listeners of a section must
only observe the widgets in its `grid`, which are deleted when the inspector shows other sections.
"""
_inspector_sections(gui::AppView, obj) = filter(e -> obj isa e.first, gui.layout.inspectors)
_inspector_sections(::AppView, ::Nothing) = Pair{Type, Function}[]

"""
    _add_inspector!(gui::AppView, T::Type, build)

Registers a type-dependent section of the inspector for objects of type `T`, see
`_inspector_sections`. Shown from the next selection of such an object.
"""
function _add_inspector!(gui::AppView, T::Type, build)
    push!(gui.layout.inspectors, T => build)
    return nothing
end
_add_inspector!(gui::LiveView, _, _) = _slot_error(gui, "inspector")

"""Deletes the type-dependent sections of the `inspector` and builds the `sections` in their place."""
function _rebuild_sections!(gui::AppView, sections)
    insp = gui.layout.inspector
    old = insp.context
    foreach(delete!, _blocks!(Any[], old))
    _GLB.remove_from_gridlayout!(_GLB.gridcontent(old))
    insp.context = GridLayout(insp.grid[5, 1]; tellwidth = false, default_rowgap = 6)
    insp.sections = sections
    insp.updates = Function[build(gui, GridLayout(insp.context[i, 1]; tellwidth = false))
                            for (i, (_, build)) in enumerate(sections)]
    rowsize!(insp.grid, 5, isempty(sections) ? Fixed(0) : Auto())
    return nothing
end

"""
    _refresh_inspector!(gui::AppView)

Shows the selected object of the `gui` in the inspector (or the summary of the live view): header,
properties and type-dependent sections, which are rebuilt only if the object needs other ones. Not
called per frame, but after the selection changed, a move and a solve. A collapsed inspector is not
updated, it is refreshed when it is shown again.
"""
function _refresh_inspector!(gui::AppView; force::Bool = false)
    layout = gui.layout
    layout.right.shown || return nothing
    insp = layout.inspector
    obj = gui.controls.selected[]
    if insp.fresh || obj !== insp.shown
        insp.fresh = false
        insp.shown = obj
        _show_header!(gui, obj)
        sections = _inspector_sections(gui, obj)
        sections == insp.sections || _rebuild_sections!(gui, sections)
    end
    _set_rows!(insp.list, _inspector_rows(gui, obj))
    foreach(update! -> update!(obj), insp.updates)
    return nothing
end

"""Sets the icon, name and type of the header of the inspector for `obj` (`nothing`: no selection)."""
function _show_header!(gui::AppView, obj)
    insp, t = gui.layout.inspector, gui.layout.theme
    kind = isnothing(obj) ? :system : _tree_kind(obj)
    insp.icon[] = _icon(kind)
    insp.icon_color[] = RGBAf(Makie.to_color(isnothing(obj) ? t.muted : _tree_marker_color(t, kind)))
    name = isnothing(obj) ? "No selection" : _label(gui, obj)
    type = isnothing(obj) ? "click an object to inspect it" : string(nameof(typeof(obj)))
    # Labels do not ellipsize, the header is as wide as the sidebar minus the icon
    w = Makie.widths(insp.grid.layoutobservables.computedbbox[])[1] - 32
    font(label) = _tree_font(label.blockscene, label.font[])
    _set_text!(insp.name, _fit_text(insp.widths[1], font(insp.name), 14, name, w))
    _set_text!(insp.type, _fit_text(insp.widths[2], font(insp.type), 12, type, w))
    return nothing
end

"""Returns the rows of the property list for `obj`, see `BeamletOptics.properties`."""
function _inspector_rows(::AppView, obj)
    props = try
        BMO.properties(obj)
    catch e
        # a failing `properties` method of a user type must not break the live view
        Pair{String, Any}["Error" => sprint(showerror, e)]
    end
    return [_property_row(name, value) for (name, value) in props if !(name in _INSPECTOR_SKIPPED)]
end

"""Summary of the live view, shown without a selection."""
function _inspector_rows(gui::AppView, ::Nothing)
    objects = sum(h -> length(h.handles), gui.system_handles; init = 0)
    rows = Tuple{String, String}[
        ("Systems", string(length(gui.system_handles))), ("Objects", string(objects)),
        ("Sources", string(length(_sources(gui)))), ("Detector panels", string(length(gui.panels))),
        ("Clip planes", string(length(gui.clip_planes))),
        ("Last trace", gui.solve_time > 0 ? _ms_string(gui.solve_time) : "–")]
    return rows
end

"""Sets the mode of the controls of the `gui`, like the key `m`."""
function _set_mode!(gui::LiveView, mode::Symbol)
    ctrl = gui.controls
    ctrl.mode[] == mode && return nothing
    ctrl.mode[] = mode
    _update_selection_box!(ctrl)
    _update_help!(ctrl)
    gui.status.text[] = "$mode mode, step: $(_step_string(mode, ctrl.fine_step, ctrl.fine_angle))"
    return nothing
end

"""
Connects the inspector of the `gui`: the mode control in both directions with the `mode` of the
controls, and a refresh when the right sidebar is shown again.
"""
function _connect_inspector!(gui::AppView)
    layout = gui.layout
    insp = layout.inspector
    ctrl = gui.controls
    listeners = ctrl.listeners
    sel = insp.mode.selected
    push!(listeners, on(m -> _set_mode!(gui, m), sel))
    push!(listeners, on(m -> (sel[] == m || (sel[] = m)), ctrl.mode; update = true))
    push!(listeners, on(v -> v && _refresh_inspector!(gui), layout.collapse.right.active))
    push!(listeners, on(s -> _set_step!(gui, s), gui.step_box.stored_string))
    push!(listeners, on(_ -> _keep_keyboard!(gui), gui.step_box.focused))
    _refresh_inspector!(gui)
    return nothing
end

# The step box of the inspector takes the keyboard like the boxes of the cards, see `_typing`
_layout_boxes(gui::AppView) = (gui.step_box,)

#=
Type-dependent sections
=#

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

"""
The section "Detector panel" of a `Detector`: the mode and the color scale of its panel, see
`_set_panel_options!`. A note names the panel; without a panel, the controls have no effect.
"""
function _detector_section(gui::AppView, grid::GridLayout)
    t = gui.layout.theme
    Label(grid[1, 1:2], "Detector panel"; halign = :left, font = :bold, fontsize = 12,
        tellwidth = false)
    note = Label(grid[2, 1:2], "no panel, see the detectors kwarg"; halign = :left,
        color = t.muted, fontsize = 12, tellwidth = false)
    Label(grid[3, 1], "mode"; halign = :left, fontsize = 11, color = t.muted)
    mode = _Segmented(grid[3, 2], [:auto => "Auto", :spot => "Spot", :intensity => "Intensity"];
        theme = t)
    Label(grid[4, 1], "scale"; halign = :left, fontsize = 11, color = t.muted)
    scale = _Segmented(grid[4, 2], [:linear => "Linear", :log => "Log"]; theme = t)
    rowgap!(grid, 4)
    colgap!(grid, 8)
    panel = Ref{Any}(nothing)
    on(m -> isnothing(panel[]) || _set_panel_options!(gui, panel[]; mode = m), mode.selected)
    on(s -> isnothing(panel[]) || _set_panel_options!(gui, panel[]; colorscale = s), scale.selected)
    return function (pd)
        p = panel[] = _panel_of(gui, pd)
        _set_text!(note, isnothing(p) ? "no panel, see the detectors kwarg" : "shown as \"$(p.name)\"")
        isnothing(p) && return nothing
        mode.selected[] == p.mode || (mode.selected[] = p.mode)
        scale.selected[] == p.colorscale || (scale.selected[] = p.colorscale)
        return nothing
    end
end

#=
Properties of the objects of the live view
=#

BMO.properties(p::LiveClipPlane) = Pair{String, Any}["Type" => "Clip plane",
    "Position [m]" => collect(Float64, p.pos), "Normal" => collect(Float64, _normal(p)),
    "Size [m]" => p.size]
