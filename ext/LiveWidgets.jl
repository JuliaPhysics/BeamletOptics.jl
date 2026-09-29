#=
Generic widgets of the live view: the property list and the segmented control
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
    _PropertyList(parent; label_color, value_color, line_color, attributes...)

A list of `label  value` rows at the grid position `parent`, the label left in `label_color`, the
value right-aligned in `value_color`, rows separated by thin lines. Drawn with a constant number of
plots (one `text!` for all labels, one for all values, one `linesegments!` for the lines) in the
scene of an invisible `Box` that holds the place in the layout; the box is as high as the rows and
takes the `attributes`, by default `tellwidth = false` (it fills the width of its column).
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

function _PropertyList(parent; label_color, value_color, line_color, attributes...)
    box = Box(parent; visible = false, height = 0, tellwidth = false, attributes...)
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
    _Segmented(parent, options; theme, selected = first(options).first, attributes...)

A segmented control at the grid position `parent`: a row of text buttons, one per option
`key => label`, of which the one of `selected[]` (a key) is highlighted with the accent colors of
the `theme` tokens. A click sets `selected`, which can also be set from code. The `attributes` go
to the `GridLayout` of the buttons, by default `tellwidth = false`.
"""
struct _Segmented
    grid::GridLayout
    buttons::Vector{Button}
    keys::Vector{Symbol}
    selected::Observable{Symbol}
end

function _Segmented(parent, options::Vector{Pair{Symbol, String}}; theme,
        selected::Symbol = first(options).first, attributes...)
    t = theme
    grid = GridLayout(parent; default_colgap = 2, halign = :left, tellwidth = false, attributes...)
    sel = Observable(selected)
    ks = first.(options)
    # the colors of the theme also where the figure has Makie's theme, e.g. on a floating card
    buttons = [Button(grid[1, i]; label, fontsize = 12, padding = (7, 7, 4, 4), cornerradius = 4,
                   buttoncolor_hover = t.hover, buttoncolor_active = t.accent_soft, strokewidth = 1)
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

