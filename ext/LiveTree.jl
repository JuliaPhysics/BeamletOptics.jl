# Object tree of the app layout of the live view: a scrollable list of rows with expanders,
# visibility toggles (eyes), type markers and labels, drawn with a constant number of plots

using Makie: Observable, Observables, events, on, off, lift, Consume, Mouse, Point2f, Vec2f, Rect2f,
             poly!, linesegments!, lines!, scatter!, text!, RGBf, RGBAf

# Geometry of a row in pixels, from the left edge of the tree
const _TREE_PAD = 4.0f0            # top and bottom padding of the list
const _TREE_LEFT = 4.0f0           # left padding
const _TREE_INDENT = 16.0f0        # indentation per depth level
const _TREE_EXPANDER = 16.0f0      # width of the expander column
const _TREE_EYE = 20.0f0           # width of the eye column
const _TREE_MARKER = 20.0f0        # width of the type marker column
const _TREE_LABEL_GAP = 2.0f0      # gap between the type marker column and the label
const _TREE_RIGHT = 10.0f0         # right padding of the labels (room for the scroll bar)
const _TREE_SCROLLBAR = 4.0f0      # width of the scroll bar
const _TREE_SCROLL_ROWS = 3        # rows per step of the mouse wheel
const _TREE_ACCENT = 3.0f0         # width of the accent stripe of the selected row

"""
    _TreeRow(key, label, depth, kind, expandable, expanded, visible)

A row of an [`_ObjectTree`](@ref), plain data without any reference to the tree.

# Fields

- `key`: any value that identifies the row, e.g. the object it stands for. It is passed to the
  `clicked`, `eye_clicked` and `expand_clicked` observables of the tree and compared with
  `isequal` by [`_set_selected!`](@ref)
- `label`: the displayed name, ellipsized if it does not fit
- `depth`: the indentation level, `0` for top-level rows
- `kind`: the type of the row, mapped to a marker and its color by the tree, e.g. `:lens`
- `expandable`: whether the row shows an expander (e.g. a group with children)
- `expanded`: the state of the expander
- `visible`: the state of the eye, `nothing` shows no eye (e.g. for rows that cannot be hidden).
  Rows with `visible == false` are drawn muted.
"""
struct _TreeRow
    key::Any
    label::String
    depth::Int
    kind::Symbol
    expandable::Bool
    expanded::Bool
    visible::Union{Bool, Nothing}
end

"""
    _ObjectTree(parent; kwargs...)

An object tree like the feature tree of a CAD program, placed at a `GridPosition` `parent`
(e.g. `fig[1, 1]`) like a Makie block. It fills its layout cell and follows its size; with
`width` or `height` it requests a fixed size from the layout instead.

The tree only displays rows, it does not know what they stand for: the caller sets the rows with
[`_set_rows!`](@ref) (only the rows that are currently shown, i.e. the children of expanded rows)
and reacts to clicks through the observables below, e.g. by toggling the visibility of an
object and setting the rows again. [`_set_selected!`](@ref) highlights a row.

Each row shows, from left to right: the indentation by depth, an expander (if `expandable`), an
eye (if `visible !== nothing`), the type marker of its `kind` and the label, which is ellipsized
to the width of the tree. The mouse wheel scrolls the rows while the mouse is over the tree; the
scroll events are then consumed, so that a 3D scene behind the tree does not zoom. A thin scroll
bar shows the position if the rows do not fit.

The tree is drawn in a child scene with a pixel camera, clipped to the layout cell, with a
constant number of plots: only the rows in view are passed to the plots, so the cost of drawing
and scrolling does not depend on the number of rows. The tree updates its plots only on
[`_set_rows!`](@ref), [`_set_selected!`](@ref), scrolling, clicks and resizing; it has no tick or
mouse move listeners.

# Keyword arguments

- `marker = kind -> :circle`: the type marker (anything `scatter` accepts) of a row kind
- `marker_color = kind -> icon_color`: the color of the type marker of a row kind
- `eye_marker = v -> v ? :circle : :xcross`: the eye marker of the visibility state `v`
- `expand_marker = e -> e ? :dtriangle : :rtriangle`: the expander marker of the state `e`
- `row_height = 22`, `fontsize = 13`, `font = :regular`: in pixels
- `marker_size = 12`, `eye_size = 10`, `expand_size = 8`: marker sizes in pixels
- `width = nothing`, `height = nothing`: requested size, `nothing` fills the layout cell
- `guides = true`: thin vertical lines that connect the rows of a level
- colors: `background`, `text_color`, `muted_color` (text and markers of hidden rows),
  `icon_color` (eyes), `expander_color`, `selection_color` (background of the selected row),
  `accent_color` (stripe at the left of the selected row), `guide_color`, `scrollbar_color`

# Fields

- `clicked`, `eye_clicked`, `expand_clicked`: `Observable{Any}`, set to the key of the row whose
  label (or type marker), eye or expander was clicked. A click fires on the release of the left
  mouse button over the same part of the same row as the press.
- `scene`: the child scene that holds the plots, its viewport is the layout cell
- `box`: the invisible `Box` that holds the place of the tree in the layout
- `rows`: the rows set by [`_set_rows!`](@ref)
- `selected`: the key of the selected row, `nothing` if none
- `offset`: the scroll position, in pixels from the top
- `plots`: the plots of the tree, a `NamedTuple`
"""
mutable struct _ObjectTree
    const box::Makie.Box
    const scene::Makie.Scene
    const clicked::Observable{Any}
    const eye_clicked::Observable{Any}
    const expand_clicked::Observable{Any}
    rows::Vector{_TreeRow}
    selected::Any
    offset::Float32
    # Row and part (:label, :eye, :expander) of the last press of the left mouse button
    pressed::Union{Nothing, Tuple{Int, Symbol}}
    const plots::NamedTuple
    const listeners::Vector{Observables.ObserverFunction}
    # Appearance
    const marker::Function
    const marker_color::Function
    const eye_marker::Function
    const expand_marker::Function
    const row_height::Float32
    const fontsize::Float32
    const font::Any
    const guides::Bool
    const text_color::RGBAf
    const muted_color::RGBAf
    const icon_color::RGBAf
    const expander_color::RGBAf
    const guide_color::RGBAf
    # Widths of labels in pixels, see `_label_width`
    const label_widths::Dict{String, Float32}
end

function _ObjectTree(parent::Makie.GridPosition;
        row_height::Real = 22, fontsize::Real = 13, font = :regular,
        marker_size::Real = 12, eye_size::Real = 10, expand_size::Real = 8,
        width = nothing, height = nothing, guides::Bool = true,
        background = RGBf(0.985, 0.985, 0.99),
        text_color = RGBf(0.13, 0.14, 0.16),
        muted_color = RGBf(0.62, 0.63, 0.66),
        icon_color = RGBf(0.52, 0.54, 0.58),
        expander_color = RGBf(0.45, 0.47, 0.50),
        selection_color = RGBf(0.84, 0.90, 0.98),
        accent_color = RGBf(0.16, 0.45, 0.85),
        guide_color = RGBf(0.86, 0.87, 0.89),
        scrollbar_color = RGBAf(0.35, 0.37, 0.40, 0.4),
        marker = kind -> :circle,
        marker_color = kind -> icon_color,
        eye_marker = v -> v ? :circle : :xcross,
        expand_marker = e -> e ? :dtriangle : :rtriangle)
    # Holds the place in the layout, the tree is drawn in a child scene on top of it
    box = Makie.Box(parent; visible = false, width, height)
    topscene = box.blockscene
    viewport = lift(Makie.round_to_IRect2D, topscene, box.layoutobservables.computedbbox)
    scene = Makie.Scene(topscene, viewport; camera = Makie.campixel!, clear = true,
        backgroundcolor = RGBAf(Makie.to_color(background)))
    common = (; inspectable = false, space = :pixel)
    rects = [Rect2f(0, 0, 0, 0), Rect2f(0, 0, 0, 0)]
    selection = poly!(scene, rects; common..., visible = false, strokewidth = 0,
        color = [RGBAf(Makie.to_color(selection_color)), RGBAf(Makie.to_color(accent_color))])
    guide = linesegments!(scene, Point2f[]; common..., color = guide_color, linewidth = 1)
    expanders = scatter!(scene, Point2f[]; common..., marker = :rtriangle,
        markersize = expand_size, color = expander_color, strokewidth = 0)
    eyes = scatter!(scene, Point2f[]; common..., marker = :circle, markersize = eye_size,
        color = RGBAf[], strokewidth = 0)
    markers = scatter!(scene, Point2f[]; common..., marker = :circle, markersize = marker_size,
        color = RGBAf[], strokewidth = 0)
    labels = text!(scene, Point2f[]; common..., text = String[], color = RGBAf[], font,
        fontsize, align = (:left, :center), markerspace = :pixel)
    scrollbar = lines!(scene, [Point2f(0), Point2f(0)]; common..., visible = false,
        color = scrollbar_color, linewidth = _TREE_SCROLLBAR, linecap = :round)
    plots = (; selection, guide, expanders, eyes, markers, labels, scrollbar)
    tree = _ObjectTree(box, scene, Observable{Any}(nothing), Observable{Any}(nothing),
        Observable{Any}(nothing), _TreeRow[], nothing, 0.0f0, nothing, plots,
        Observables.ObserverFunction[], marker, marker_color, eye_marker, expand_marker,
        Float32(row_height), Float32(fontsize), _tree_font(scene, font), guides,
        RGBAf(Makie.to_color(text_color)), RGBAf(Makie.to_color(muted_color)),
        RGBAf(Makie.to_color(icon_color)), RGBAf(Makie.to_color(expander_color)),
        RGBAf(Makie.to_color(guide_color)), Dict{String, Float32}())
    push!(tree.listeners, on(_ -> _redraw!(tree), viewport))
    ev = events(scene)
    # Before the camera and the kinematic controls of a 3D scene (priority 200)
    push!(tree.listeners, on(ev.scroll, priority = 300) do (_, dy)
        _mouse_in_tree(tree) || return Consume(false)
        _scroll!(tree, -dy * _TREE_SCROLL_ROWS * tree.row_height)
        return Consume(true)
    end)
    push!(tree.listeners, on(ev.mousebutton, priority = 300) do event
        event.button == Mouse.left || return Consume(false)
        if event.action == Mouse.press
            _mouse_in_tree(tree) || return Consume(false)
            # Presses beside the rows are consumed, too, and so are their releases
            tree.pressed = something(_hit(tree), (0, :none))
            return Consume(true)
        elseif event.action == Mouse.release && !isnothing(tree.pressed)
            pressed, tree.pressed = tree.pressed, nothing
            hit = _mouse_in_tree(tree) ? _hit(tree) : nothing
            # The press consumed, so is the release
            (isnothing(hit) || hit != pressed) && return Consume(true)
            i, part = hit
            key = tree.rows[i].key
            part === :eye ? (tree.eye_clicked[] = key) :
            part === :expander ? (tree.expand_clicked[] = key) : (tree.clicked[] = key)
            return Consume(true)
        end
        return Consume(false)
    end)
    _redraw!(tree)
    return tree
end

"""
    _set_rows!(tree::_ObjectTree, rows)

Shows `rows` (a vector of [`_TreeRow`](@ref)) in `tree`, e.g. after a row was expanded or hidden.
The scroll position is kept (within the new content), and so is the selection, if a row with the
selected key is among `rows`.
"""
function _set_rows!(tree::_ObjectTree, rows::AbstractVector{_TreeRow})
    tree.rows = collect(rows)
    tree.pressed = nothing
    # Labels of rows that are gone are measured again when needed
    length(tree.label_widths) > 4 * length(rows) + 1000 && empty!(tree.label_widths)
    _redraw!(tree)
    return tree
end

"""
    _set_selected!(tree::_ObjectTree, key)

Highlights the row with the key `key` (compared with `isequal`) and scrolls it into view;
`nothing` clears the selection. A key that is not among the rows of the tree (e.g. an object in a
collapsed group) highlights no row, but is kept, such that the row is highlighted once it is
shown by [`_set_rows!`](@ref).
"""
function _set_selected!(tree::_ObjectTree, key)
    tree.selected = key
    i = isnothing(key) ? nothing : findfirst(r -> isequal(r.key, key), tree.rows)
    if !isnothing(i)
        h = _tree_size(tree)[2]
        top = _TREE_PAD + (i - 1) * tree.row_height
        if top < tree.offset
            tree.offset = top - _TREE_PAD
        elseif top + tree.row_height > tree.offset + h
            tree.offset = top + tree.row_height + _TREE_PAD - h
        end
    end
    _redraw!(tree)
    return tree
end

function Base.close(tree::_ObjectTree)
    foreach(off, tree.listeners)
    empty!(tree.listeners)
    isnothing(tree.scene.parent) || Makie.free(tree.scene)
    delete!(tree.box)
    return nothing
end

_tree_size(tree::_ObjectTree) = Vec2f(Makie.widths(tree.scene.viewport[]))

_content_height(tree::_ObjectTree) = 2 * _TREE_PAD + length(tree.rows) * tree.row_height

_max_offset(tree::_ObjectTree) = max(0.0f0, _content_height(tree) - _tree_size(tree)[2])

"""
    _row_y(tree, i)

Returns the vertical center of the `i`-th row in the pixel coordinates of `tree.scene`.
"""
_row_y(tree::_ObjectTree, i::Integer) =
    _tree_size(tree)[2] - _TREE_PAD - (i - 0.5f0) * tree.row_height + tree.offset

"""
    _row_columns(tree, row)

Returns the horizontal centers of the expander, the eye and the type marker of `row` and the left
edge of its label, in pixels. The columns of the expander and the eye are reserved in every row,
so that the markers and labels of a level align.
"""
function _row_columns(::_ObjectTree, row::_TreeRow)
    x = _TREE_LEFT + row.depth * _TREE_INDENT
    expander = x + _TREE_EXPANDER / 2
    eye = x + _TREE_EXPANDER + _TREE_EYE / 2
    marker = x + _TREE_EXPANDER + _TREE_EYE + _TREE_MARKER / 2
    label = x + _TREE_EXPANDER + _TREE_EYE + _TREE_MARKER + _TREE_LABEL_GAP
    return (; expander, eye, marker, label)
end

function _mouse_in_tree(tree::_ObjectTree)
    vp = tree.scene.viewport[]
    return Point2f(events(tree.scene).mouseposition[]) in Rect2f(vp)
end

"""
    _hit(tree)

Returns the row index and the part (`:expander`, `:eye` or `:label`) under the mouse, `nothing`
if the mouse is not over a row. A row without an expander or an eye has only a label part.
"""
function _hit(tree::_ObjectTree)
    vp = tree.scene.viewport[]
    p = Point2f(events(tree.scene).mouseposition[]) - Point2f(Makie.origin(vp))
    h = Makie.widths(vp)[2]
    i = floor(Int, (h - p[2] + tree.offset - _TREE_PAD) / tree.row_height) + 1
    1 <= i <= length(tree.rows) || return nothing
    row = tree.rows[i]
    c = _row_columns(tree, row)
    if row.expandable && abs(p[1] - c.expander) <= _TREE_EXPANDER / 2
        return (i, :expander)
    elseif !isnothing(row.visible) && abs(p[1] - c.eye) <= _TREE_EYE / 2
        return (i, :eye)
    end
    return (i, :label)
end

function _scroll!(tree::_ObjectTree, delta::Real)
    offset = clamp(tree.offset + Float32(delta), 0.0f0, _max_offset(tree))
    offset == tree.offset && return nothing
    tree.offset = offset
    _redraw!(tree)
    return nothing
end

"""
    _label_width(tree, s)

Returns the width of the label `s` in pixels. The widths are cached, a label is measured once.
"""
_label_width(tree::_ObjectTree, s::String) = _text_width(tree.label_widths, tree.font, tree.fontsize, s)

"""
    _fit_label(tree, s, width)

Returns `s`, or its longest prefix followed by "…" if `s` is wider than `width` pixels.
"""
_fit_label(tree::_ObjectTree, s::String, width::Real) =
    _fit_text(tree.label_widths, tree.font, tree.fontsize, s, width)

"""
    _text_width(cache, font, fontsize, s)

Returns the width of the text `s` in pixels, measured once and then taken from the `cache` (a
`Dict{String, Float32}` per font and size).
"""
function _text_width(cache::Dict{String, Float32}, font, fontsize, s::AbstractString)
    return get!(cache, s) do
        Float32(Makie.widths(Makie.text_bb(s, font, fontsize))[1])
    end
end

"""
    _fit_text(cache, font, fontsize, s, width)

Returns `s`, or its longest prefix followed by "…" if `s` is wider than `width` pixels, see
`_text_width`.
"""
function _fit_text(cache::Dict{String, Float32}, font, fontsize, s::AbstractString, width::Real)
    s = String(s)
    _text_width(cache, font, fontsize, s) <= width && return s
    n = length(s)
    lo, hi = 0, n - 1
    # Longest prefix with an ellipsis that fits, found by bisection
    while lo < hi
        m = (lo + hi + 1) ÷ 2
        if _text_width(cache, font, fontsize, first(s, m) * "…") <= width
            lo = m
        else
            hi = m - 1
        end
    end
    return lo == 0 ? "…" : rstrip(first(s, lo)) * "…"
end

"""
    _redraw!(tree)

Passes the rows in view to the plots of `tree`, after clamping the scroll position.
"""
function _redraw!(tree::_ObjectTree)
    w, h = _tree_size(tree)
    tree.offset = clamp(tree.offset, 0.0f0, _max_offset(tree))
    rows, rh = tree.rows, tree.row_height
    first_row = max(1, floor(Int, (tree.offset - _TREE_PAD) / rh) + 1)
    last_row = min(length(rows), ceil(Int, (tree.offset + h - _TREE_PAD) / rh))
    in_view = first_row:last_row
    n = length(in_view)
    label_pos = Vector{Point2f}(undef, n)
    label_text = Vector{String}(undef, n)
    label_color = Vector{RGBAf}(undef, n)
    marker_pos = Vector{Point2f}(undef, n)
    marker_shape = Vector{Any}(undef, n)
    marker_color = Vector{RGBAf}(undef, n)
    expander_pos, expander_shape = Point2f[], Any[]
    eye_pos, eye_shape, eye_color = Point2f[], Any[], RGBAf[]
    guides = Point2f[]
    scrollbar = _max_offset(tree) > 0
    label_right = w - (scrollbar ? _TREE_RIGHT : _TREE_RIGHT / 2)
    selected = nothing
    for (k, i) in enumerate(in_view)
        row = rows[i]
        y = _row_y(tree, i)
        c = _row_columns(tree, row)
        muted = row.visible === false
        label_pos[k] = Point2f(c.label, y)
        label_text[k] = _fit_label(tree, row.label, label_right - c.label)
        label_color[k] = muted ? tree.muted_color : tree.text_color
        marker_pos[k] = Point2f(c.marker, y)
        marker_shape[k] = tree.marker(row.kind)
        marker_color[k] = muted ? tree.muted_color : RGBAf(Makie.to_color(tree.marker_color(row.kind)))
        if row.expandable
            push!(expander_pos, Point2f(c.expander, y))
            push!(expander_shape, tree.expand_marker(row.expanded))
        end
        if !isnothing(row.visible)
            push!(eye_pos, Point2f(c.eye, y))
            push!(eye_shape, tree.eye_marker(row.visible))
            push!(eye_color, muted ? tree.muted_color : tree.icon_color)
        end
        # A vertical line per ancestor level, under the expanders of the ancestors
        if tree.guides
            for d in 0:(row.depth - 1)
                x = round(_TREE_LEFT + d * _TREE_INDENT + _TREE_EXPANDER / 2) + 0.5f0
                push!(guides, Point2f(x, y + rh / 2), Point2f(x, y - rh / 2))
            end
        end
        !isnothing(tree.selected) && isequal(row.key, tree.selected) && (selected = y)
    end
    p = tree.plots
    Makie.update!(p.labels; arg1 = label_pos, text = label_text, color = label_color)
    Makie.update!(p.markers; arg1 = marker_pos, marker = _tree_markers(marker_shape), color = marker_color)
    Makie.update!(p.expanders; arg1 = expander_pos, marker = _tree_markers(expander_shape))
    Makie.update!(p.eyes; arg1 = eye_pos, marker = _tree_markers(eye_shape), color = eye_color)
    Makie.update!(p.guide; arg1 = guides)
    if isnothing(selected)
        p.selection.visible = false
    else
        bottom = selected - rh / 2
        Makie.update!(p.selection; visible = true,
            arg1 = [Rect2f(0, bottom, w, rh), Rect2f(0, bottom, _TREE_ACCENT, rh)])
    end
    if scrollbar
        content = _content_height(tree)
        thumb = max(20.0f0, h * h / content)
        top = h - (h - thumb) * tree.offset / _max_offset(tree)
        x = w - _TREE_SCROLLBAR / 2 - 2
        r = _TREE_SCROLLBAR / 2
        Makie.update!(p.scrollbar; visible = true,
            arg1 = [Point2f(x, top - r - 1), Point2f(x, top - thumb + r + 1)])
    else
        p.scrollbar.visible = false
    end
    return nothing
end

# Markers as a vector with a concrete element type, symbols as their `BezierPath`s, such that
# symbols and custom markers (e.g. icons) can be mixed
_tree_marker_shape(s::Symbol) = get(Makie.default_marker_map(), s, s)
_tree_marker_shape(s) = s

function _tree_markers(shapes::Vector{Any})
    isempty(shapes) && return Makie.BezierPath[]
    return identity.(map(_tree_marker_shape, shapes))
end

# The font of the labels, a symbol names a font of the theme, e.g. `:regular`
_tree_font(scene::Makie.Scene, font::Symbol) = Makie.to_font(scene.theme[:fonts], font)
_tree_font(::Makie.Scene, font) = Makie.to_font(font)
