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
Docked cards: the card of the selection in the inspector and the pinned cards below it, see
`_AbstractCard`
=#

"""
    _DockedCard

A card docked in the "Properties" section of the app layout, the counterpart of the floating
`_ComponentCard` (see `_AbstractCard`): the widgets declared by [`card_actions`](@ref) (in
`actions`, in the head) and [`card_rows`](@ref) (in `rows`, below the head), built by the same code
as the floating cards (see `_build_content!`), but in the colors of the `theme` of the app and with
the textboxes and sliders filling the width of the sidebar (see `_cell_attributes`). It is either

- the card of the selection (`pinned = false`), whose head is the header of the inspector, see
  `_Inspector`, or
- a card `pinned` to the object `obj`, stacked below the inspector with its own `head` (icon,
  title, pin and collapse chevron, built by the shared parts of `LiveCard.jl`), see
  `_dock_pinned!`; a `collapsed` card shows only its head and its actions.

`header` and `parent` hold the layouts of the head and of the card: the actions are placed in the
rows `actions_rows` of column 3 of the `header`, the rows in the row `rows_row` of the `parent`.
The fields `widgets` to `pose` are those of `_ComponentCard`.
"""
mutable struct _DockedCard <: _AbstractCard
    const header::GridLayout
    const parent::GridLayout
    const theme::NamedTuple
    const actions_rows::UnitRange{Int}
    const rows_row::Int
    actions::GridLayout
    rows::GridLayout
    widgets::Vector{Tuple{Any, CardWidget}}
    blocks::Vector{Any}
    textboxes::Vector{Textbox}
    listeners::Vector{Any}
    content_key::Any
    refreshing::Bool
    pose::Any
    pinned::Bool
    obj::Any
    collapsed::Bool
    # the parts of the head of a pinned card, `nothing` for the card of the selection
    head::Any
end

# Positions of the actions in the header and of the rows in the card
_docked_actions(c::_DockedCard) = _docked_actions(c.header, c.actions_rows)
_docked_actions(header::GridLayout, rows::UnitRange{Int}) = GridLayout(header[rows, 3]; default_colgap = 4)
_docked_rows_layout(parent::GridLayout, row::Int) =
    GridLayout(parent[row, 1]; default_rowgap = 4, tellwidth = false)

function _DockedCard(header::GridLayout, parent::GridLayout, theme::NamedTuple;
        actions_rows::UnitRange{Int} = 1:2, rows_row::Int = 2)
    return _DockedCard(header, parent, theme, actions_rows, rows_row,
        _docked_actions(header, actions_rows), _docked_rows_layout(parent, rows_row),
        Tuple{Any, CardWidget}[], Any[], Textbox[], Any[], nothing, false, nothing, false, nothing,
        false, nothing)
end

function _new_parts!(c::_DockedCard)
    for part in (c.actions, c.rows)
        _GLB.remove_from_gridlayout!(_GLB.gridcontent(part))
    end
    c.actions, c.rows = _docked_actions(c), _docked_rows_layout(c.parent, c.rows_row)
    return nothing
end

_card_object(gui::LiveView, c::_DockedCard) = c.pinned ? c.obj : gui.controls.selected[]
_card_boxes(c::_DockedCard) = c.textboxes
# A collapsed card shows only its actions
_declarations(c::_DockedCard, obj) = (card_actions(obj), c.collapsed ? () : card_rows(obj))

# The widgets take the theme of the figure, texts and axis colors from the tokens of the app
_card_style(c::_DockedCard, ::Type{Label}) = (; color = c.theme.text, fontsize = 12)
_card_style(::_DockedCard, ::Type{Textbox}) = (; fontsize = 12, textpadding = (5, 5, 4, 4))
_card_style(::_DockedCard, ::Type{Button}) = (; fontsize = 12, padding = (7, 7, 4, 4))
_card_style(::_DockedCard, ::Type) = (;)
_card_value(c::_DockedCard, a::_AxisColor) = c.theme.gizmo[a.k]
_row_attributes(::_DockedCard) = (; default_colgap = 6, tellwidth = false, halign = :left)

# Textboxes and sliders fill the width of the sidebar instead of the width declared for the card
_host_attributes(::_DockedCard, T::Type, attributes::NamedTuple) = _fill_width(T, attributes)
_fill_width(::Type{<:Union{Textbox, Slider}}, attributes) =
    merge(attributes, (; width = Relative(1), tellwidth = false))
_fill_width(::Type, attributes) = attributes

#=
Inspector
=#

"""
    _Inspector

The property inspector of the app layout, the "PROPERTIES" section of the right sidebar, from top
to bottom:

- the header: the icon of the kind of the selected object (see `_tree_kind`), its name (see
  `_label`) and its type, the actions of its card (e.g. hide) and the `pin` toggle, which pins a
  card to the object below the inspector (see `_toggle_pin!`); the icon, the name and the pin are
  built by the shared parts of the cards, see `_card_icon!`
- the docked `card` of the selected object: the rows of [`card_rows`](@ref), e.g. the pose, see
  `_DockedCard`
- the step box and the mode as a segmented control (`mode`), bound to the `mode` of the controls
- the properties of the object, see `BeamletOptics.properties` and `_PropertyList`; without a
  selection, a summary of the live view
- the `pinned` cards in the layout `pinned_grid`, one below the other in the order of pinning,
  each with its own head, see `_dock_pinned!`. The app layout shows no floating cards: the
  "Properties" section is formed by the card of the selection and the pinned cards.

The inspector is updated on events only (see `_refresh_inspector!`): the selection, moves, solves
and inputs. The widgets of the card are only rebuilt if the selected object declares others, e.g.
when a source follows a lens; otherwise they show the values of the new object.
"""
mutable struct _Inspector
    const grid::GridLayout
    const icon::Observable{BezierPath}
    const icon_color::Observable{RGBAf}
    const name::Label
    const type::Label
    const pin::_IconToggle
    const card::_DockedCard
    const mode::_Segmented
    const list::_PropertyList
    const pinned_grid::GridLayout
    const pinned::Vector{_DockedCard}
    # the object shown (`nothing`: the summary), a flag that nothing was shown yet
    shown::Any
    fresh::Bool
    # the pinned card that was pinned or expanded last, which stays expanded, see `_fit_pinned!`
    keep::Any
    # widths of the texts of the header in pixels, per label, see `_text_width`
    const widths::NTuple{2, Dict{String, Float32}}
end

# Names of `properties` that the inspector shows elsewhere: in the header and the pose rows
const _INSPECTOR_SKIPPED = ("Type", "Position [m]")

"""
    _build_inspector!(layout::AppLayout) -> (; step_box)

Creates the "PROPERTIES" section of the right sidebar, see `_Inspector`. Returns the widgets that
are fields of `LiveView`.
"""
function _build_inspector!(layout::AppLayout)
    t = layout.theme
    g = _add_sidebar_section!(layout, :right, "Properties")
    # Header: icon, name and type, the actions of the card and the pin
    header = GridLayout(g[1, 1]; default_colgap = 6, tellwidth = false)
    icon, icon_color = _card_icon!(header[1:2, 1]; size = 22, box = 24)
    _show_kind!(icon, icon_color, t, nothing)
    name = _card_title!(header[1, 2], t; fontsize = 14, tellwidth = false)
    name.text[] = "No selection"
    type = Label(header[2, 2], " "; halign = :left, color = t.muted, fontsize = 12,
        tellwidth = false)
    pin = _card_pin!(header[1:2, 4], t; size = 24, icon_size = 18,
        tooltip = "Pin a card below the inspector", tooltip_placement = :left)
    rowgap!(header, 0)
    colsize!(header, 2, Auto(false))
    card = _DockedCard(header, g, t)
    # Step and mode
    pose = GridLayout(g[3, 1]; default_colgap = 6, default_rowgap = 2, tellwidth = false)
    Label(pose[1, 1], "step"; halign = :left, fontsize = 11, color = t.muted, tellwidth = false)
    Label(pose[1, 2:3], "mode"; halign = :left, fontsize = 11, color = t.muted, tellwidth = false)
    step_box = _step_box!(pose[2, 1], t; placeholder = "250 nm", width = 80, halign = :left)
    mode = _Segmented(pose[2, 2:3], [:move => "Move", :rotate => "Rotate"]; theme = t)
    # Properties
    Box(g[4, 1]; height = 1, color = t.border, strokewidth = 0)
    list = _PropertyList(g[5, 1]; label_color = t.muted, value_color = t.text,
        line_color = RGBAf(Makie.to_color(t.border)))
    # Pinned cards, see `_dock_pinned!`
    pinned_grid = GridLayout(g[6, 1]; default_rowgap = 10, tellwidth = false)
    # the rows of the card are empty without a selection, see `_refresh_inspector!`, and there are
    # no pinned cards yet
    rowsize!(g, 2, Fixed(0))
    rowsize!(g, 6, Fixed(0))
    rowgap!(g, 10)
    layout.inspector = _Inspector(g, icon, icon_color, name, type, pin, card, mode, list,
        pinned_grid, _DockedCard[], nothing, true, nothing, (Dict{String, Float32}(), Dict{String, Float32}()))
    return (; step_box)
end

#=
Pinned cards of the app layout, docked below the inspector
=#

"""
    _dock_pinned!(gui::AppView, obj) -> _DockedCard

Pins a card to `obj` in the app layout: a `_DockedCard` at the end of the pinned cards of the
inspector, below a line, with its own head (icon, title, actions, pin and collapse chevron, see the
shared parts `_card_icon!`, `_card_title!`, `_card_pin!` and `_card_collapse!`) and the rows of
`obj`. The pin unpins it, see `_unpin!`; the chevron collapses it to its head. Its widgets are built
by `_refresh_inspector!`.
"""
function _dock_pinned!(gui::AppView, obj)
    insp, t = gui.layout.inspector, gui.layout.theme
    g = GridLayout(insp.pinned_grid[length(insp.pinned) + 1, 1]; default_rowgap = 4, tellwidth = false)
    line = Box(g[1, 1]; height = 1, color = t.border, strokewidth = 0)
    header = GridLayout(g[2, 1]; default_colgap = 6, tellwidth = false)
    icon, icon_color = _card_icon!(header[1, 1])
    title = _card_title!(header[1, 2], t; tellwidth = false)
    pin = _card_pin!(header[1, 4], t; active = true, tooltip_placement = :left)
    collapse = _card_collapse!(header[1, 5], t; tooltip_placement = :left)
    colgap!(header, 4, 2)
    colsize!(header, 2, Auto(false))
    c = _DockedCard(header, g, t; actions_rows = 1:1, rows_row = 3)
    c.pinned, c.obj = true, obj
    listeners = Any[
        on(v -> v || _unpin!(gui, obj), pin.active),
        on(_ -> _toggle_collapsed!(gui, c), collapse.clicks)]
    c.head = (; icon, icon_color, title, pin, collapse, line, listeners, widths = Dict{String, Float32}())
    _show_kind!(icon, icon_color, t, obj)
    _show_head!(pin, collapse, true, false)
    push!(insp.pinned, c)
    rowsize!(insp.grid, 6, Auto())
    return c
end

"""Removes the pinned card `c` of the app layout of the `gui` with its widgets and listeners."""
function _remove_pinned!(gui::AppView, c::_DockedCard)
    insp = gui.layout.inspector
    foreach(tb -> tb.focused[] && Makie.defocus!(tb), c.textboxes)
    _clear_content!(c)
    foreach(off, c.head.listeners)
    # the line and the blocks of the head: icon, title, pin and chevron
    delete!(c.head.line)
    foreach(delete!, [gc.content for gc in copy(c.header.content) if gc.content isa Makie.Block])
    _GLB.remove_from_gridlayout!(_GLB.gridcontent(c.parent))
    filter!(d -> d !== c, insp.pinned)
    insp.keep === c && (insp.keep = nothing)
    # The remaining cards move up
    for (k, d) in enumerate(insp.pinned)
        insp.pinned_grid[k, 1] = d.parent
    end
    isempty(insp.pinned) || _GLB.trim!(insp.pinned_grid)
    rowsize!(insp.grid, 6, isempty(insp.pinned) ? Fixed(0) : Auto())
    return nothing
end

"""
Collapses the pinned card `c` to its head and its actions, or expands it again; the other pinned
cards collapse if the expanded card does not fit, see `_fit_pinned!`.
"""
function _toggle_collapsed!(gui::AppView, c::_DockedCard)
    _set_collapsed!(c, !c.collapsed)
    c.collapsed || (gui.layout.inspector.keep = c)
    _refresh_inspector!(gui; force = true)
    return nothing
end

function _set_collapsed!(c::_DockedCard, collapsed::Bool)
    c.collapsed = collapsed
    _show_collapsed!(c.head.collapse, collapsed)
    # the declarations change, see `_declarations`
    c.pose = nothing
    return nothing
end

"""
    _fit_pinned!(gui::AppView; keep = nothing)

Fits the "Properties" section of the `gui` into the right sidebar, which does not scroll (see
`_overflow`): first the pinned cards collapse to their heads, the oldest first, except the card
`keep` (the one pinned or expanded last); then the property list of the selection is shortened,
see `_fit_list!`. Cards are never expanded automatically, i.e. they do not change back and forth.
"""
function _fit_pinned!(gui::AppView; keep = nothing)
    for c in gui.layout.inspector.pinned
        _overflows(gui) || break
        (c === keep || c.collapsed) && continue
        _set_collapsed!(c, true)
        _refresh_pinned!(gui, c; force = true)
    end
    _fit_list!(gui)
    return nothing
end

"""
Shortens the property list of the inspector of the `gui` by the rows that do not fit into the right
sidebar (see `_overflow`), the last row reading "… n more". The list is set again in full by the
next `_refresh_inspector!`.
"""
function _fit_list!(gui::AppView)
    excess = _overflow(gui)
    excess > 0.5 || return nothing
    list = gui.layout.inspector.list
    rows = list.rows
    n = length(rows) - ceil(Int, excess / _PROPERTY_ROW)
    n >= length(rows) && return nothing
    shown = n <= 1 ? Tuple{String, String}[] : [rows[1:(n - 1)]; ("… $(length(rows) - n + 1) more", "")]
    _set_rows!(list, shown)
    return nothing
end

"""
Returns how far [px] the parts of the right sidebar of the `gui` (the sections and their titles)
reach beyond it, above and below in total: its layout pushes the sections over the toolbar and the
dock if they do not fit.
"""
function _overflow(gui::AppView)
    stack = gui.layout.right.grid
    s = stack.layoutobservables.computedbbox[]
    isempty(stack.content) && return 0.0f0
    rs = [gc.content.layoutobservables.computedbbox[] for gc in stack.content]
    above = maximum(r -> maximum(r)[2], rs) - maximum(s)[2]
    below = minimum(s)[2] - minimum(r -> minimum(r)[2], rs)
    return Float32(max(above, 0) + max(below, 0))
end
_overflows(gui::AppView) = _overflow(gui) > 0.5

"""
    _refresh_pinned!(gui::AppView, c::_DockedCard; force = false)

Shows the object of the pinned card `c`: builds its widgets when it was pinned or collapsed (see
`_build_content!`), then its title and the values of its widgets.
"""
function _refresh_pinned!(gui::AppView, c::_DockedCard; force::Bool = false)
    if c.pose === nothing
        _build_content!(gui, c, c.obj)
        c.pose = (c.obj, nothing)
        rowsize!(c.parent, c.rows_row, isempty(c.rows.content) ? Fixed(0) : Auto())
    end
    title = c.head.title
    w = Makie.widths(gui.layout.inspector.grid.layoutobservables.computedbbox[])[1] -
        _CARD_ICON - 2 * _CARD_TOOL - 20 - _actions_width(c)
    _set_text!(title, _fit_text(c.head.widths, _tree_font(title.blockscene, title.font[]),
        _CARD_TITLE_FONTSIZE, _label(gui, c.obj), w))
    _refresh_card!(gui, c; force)
    return nothing
end

function _pin!(gui::AppView, obj)
    gui.layout.inspector.keep = _dock_pinned!(gui, obj)
    _refresh_inspector!(gui; force = true)
    return nothing
end
_is_pinned(gui::AppView, obj) = any(c -> c.obj === obj, gui.layout.inspector.pinned)
function _unpin!(gui::AppView, obj)
    foreach(c -> _remove_pinned!(gui, c), filter(c -> c.obj === obj, gui.layout.inspector.pinned))
    _on_pinned!(gui)
    return nothing
end

"""
    _refresh_inspector!(gui::AppView; force = false)

Shows the selected object of the `gui` in the inspector (or the summary of the live view): header,
docked card and properties. The widgets of the card are rebuilt only if the object declares others,
see `_build_content!`; a focused textbox of the card keeps the typed text, unless `force`. Not
called per frame, but after the selection changed, a move, a solve and an input. A collapsed
inspector is not updated, it is refreshed when it is shown again. Then the pinned cards, which
collapse (except the one pinned or expanded last) while they do not fit, see `_fit_pinned!`.
"""
function _refresh_inspector!(gui::AppView; force::Bool = false)
    layout = gui.layout
    layout.right.shown || return nothing
    insp = layout.inspector
    obj = gui.controls.selected[]
    if insp.fresh || obj !== insp.shown
        insp.fresh = false
        insp.shown = obj
        _dock_card!(gui, insp.card, obj)
        _show_header!(gui, obj)
    end
    _refresh_card!(gui, insp.card; force)
    _show_pin!(gui)
    _set_rows!(insp.list, _inspector_rows(gui, obj))
    foreach(c -> _refresh_pinned!(gui, c; force), insp.pinned)
    _fit_pinned!(gui; keep = insp.keep)
    return nothing
end

"""
    _dock_card!(gui::AppView, c::_DockedCard, obj)

Builds the widgets of the docked card `c` for `obj`, see `_build_content!`, or removes them without
a selection. The row of the card in the inspector has no height without rows.
"""
function _dock_card!(::AppView, c::_DockedCard, ::Nothing)
    _clear_content!(c)
    c.pose = nothing
    rowsize!(c.parent, 2, Fixed(0))
    return nothing
end
function _dock_card!(gui::AppView, c::_DockedCard, obj)
    # A focused box of the old object would take the keyboard, and its input the new object
    foreach(tb -> tb.focused[] && Makie.defocus!(tb), c.textboxes)
    _build_content!(gui, c, obj)
    c.pose = (obj, nothing)
    rowsize!(c.parent, 2, isempty(c.rows.content) ? Fixed(0) : Auto())
    return nothing
end

"""Sets the pin of the inspector of the `gui` to whether a card is pinned to the selected object."""
function _show_pin!(gui::AppView)
    pin = gui.layout.inspector.pin
    obj = gui.controls.selected[]
    pinned = !isnothing(obj) && _is_pinned(gui, obj)
    pin.active[] == pinned || (pin.active[] = pinned)
    visible = !isnothing(obj)
    pin.box.visible[] == visible || (pin.box.visible[] = visible)
    return nothing
end

_on_pinned!(gui::AppView) = _show_pin!(gui)

"""Sets the icon, name and type of the header of the inspector for `obj` (`nothing`: no selection)."""
function _show_header!(gui::AppView, obj)
    insp, t = gui.layout.inspector, gui.layout.theme
    _show_kind!(insp.icon, insp.icon_color, t, obj)
    name = isnothing(obj) ? "No selection" : _label(gui, obj)
    type = isnothing(obj) ? "click an object to inspect it" : string(nameof(typeof(obj)))
    # Labels do not ellipsize: the room for the texts is the column of the name, between the icon
    # and the actions of the card
    w = Makie.widths(insp.grid.layoutobservables.computedbbox[])[1] - 32 - 30 -
        _actions_width(insp.card)
    font(label) = _tree_font(label.blockscene, label.font[])
    _set_text!(insp.name, _fit_text(insp.widths[1], font(insp.name), 14, name, w))
    _set_text!(insp.type, _fit_text(insp.widths[2], font(insp.type), 12, type, w))
    return nothing
end

# Width of the actions of the docked card `c` [px], with their gap
_actions_width(c::_DockedCard) =
    isempty(c.actions.content) ? 0.0f0 : Makie.widths(c.actions.layoutobservables.computedbbox[])[1] + 6

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
controls, the pin, the step box, and a refresh when the right sidebar is shown again. The widgets of
the docked card are connected when they are built, see `_build_content!`.
"""
function _connect_inspector!(gui::AppView)
    layout = gui.layout
    insp = layout.inspector
    ctrl = gui.controls
    listeners = ctrl.listeners
    sel = insp.mode.selected
    push!(listeners, on(m -> _set_mode!(gui, m), sel))
    push!(listeners, on(m -> (sel[] == m || (sel[] = m)), ctrl.mode; update = true))
    push!(listeners, on(insp.pin.active) do v
        obj = ctrl.selected[]
        (isnothing(obj) || v == _is_pinned(gui, obj)) || _toggle_pin!(gui, obj)
        return nothing
    end)
    push!(listeners, on(v -> v && _refresh_inspector!(gui), layout.collapse.right.active))
    push!(listeners, on(s -> _set_step!(gui, s), gui.step_box.stored_string))
    push!(listeners, on(_ -> _keep_keyboard!(gui), gui.step_box.focused))
    _refresh_inspector!(gui)
    return nothing
end

# The step box and the textboxes of the docked card take the keyboard like those of the floating
# cards, see `_typing`
_layout_boxes(gui::AppView) = (gui.step_box, _card_boxes(gui.layout.inspector.card)...,
    (tb for c in gui.layout.inspector.pinned for tb in _card_boxes(c))...)

# The card of the selection is docked in the inspector, only pinned cards float in the 3D view
_selection_card_shown(::AppView) = false

# The toolbar, the sidebars, the dock and the status bar surround the 3D view: a click on their
# widgets, e.g. the pin of the inspector, neither selects nor deselects (the release of the press
# would clear the selection, and with it the inspector)
_outside_view(gui::AppView) = !Makie.is_mouseinside(gui.ax.scene)

#=
Properties of the objects of the live view
=#

BMO.properties(p::LiveClipPlane) = Pair{String, Any}["Type" => "Clip plane",
    "Position [m]" => collect(Float64, p.pos), "Normal" => collect(Float64, _normal(p)),
    "Size [m]" => p.size]
