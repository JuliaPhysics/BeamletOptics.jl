using Makie: Figure, Observable, Point2f, Vec2f, Rect2f, GridLayout, Textbox, Button, Label, Box, Slider,
             Toggle, Menu

# Component cards of the live view: the controls of an object, next to it in the 3D view

# z translation of the scene of a card: GLMakie draws the plots in the order of this value, the
# card comes after the 3D scene, and its depth (≈ 0.005–0.02) lies in front of it. Each card gets
# its own value, `_CARD_Z0` for the first, `_CARD_DZ` more for each further one up to `_CARD_Z`,
# such that a newer card is drawn over an older one as a whole. Plus the offsets of the widgets,
# it stays within the clip range ±10000 of the pixel camera
const _CARD_Z = 9900.0f0
const _CARD_Z0 = 9600.0f0
const _CARD_DZ = 10.0f0

"""The z translation of the `k`-th card (from 1), see `_CARD_Z`."""
_card_z(k::Int) = min(_CARD_Z, _CARD_Z0 + (k - 1) * _CARD_DZ)
# Distance of the card from the bounding box of its object [px]
const _CARD_GAP = 12.0f0
# Minimum distance of the card from the edges of the 3D view and from the view cube [px]
const _CARD_MARGIN = 8.0f0
# Distance of the widgets from the edge of the background and between the parts of the card [px]
const _CARD_PADDING = 8.0f0
# Suggested bounding box of the parts that are not shown: hidden widgets still take clicks within
# their bounding box, hence they are moved far outside of the figure
const _CARD_AWAY = Rect2f(-1.0f5, -1.0f5, 0, 0)
# Opaque, unlike the progress window, such that objects behind the widgets do not shine through
const _CARD_BACKGROUND = RGBAf(0.1, 0.1, 0.12, 1)
# Line from the card to its object and the dot at the object
const _CARD_LINK_COLOR = RGBAf(0.95, 0.95, 0.95, 0.85)
# Font size of the cards, smaller than the theme, such that several cards fit into the view
const _CARD_FONTSIZE = 14

"""
    _AbstractCard

Host of the widgets that [`card_actions`](@ref) and [`card_rows`](@ref) declare for an object, see
`_build_content!`: the floating `_ComponentCard` in the 3D view, or a card docked in a sidebar,
e.g. `_DockedCard` in the inspector of the app layout. A host has the fields

- `actions::GridLayout` (the buttons of the actions, in one row) and `rows::GridLayout` (a row per
  `CardRow`), which `_new_parts!` replaces when the widgets are rebuilt
- `widgets` (each block with its `CardWidget`), `blocks`, `textboxes`, `listeners`, `content_key`,
  `refreshing` and `pose`, see `_ComponentCard`

and implements `_new_parts!(c)`, `_card_object(gui, c)`, `_card_boxes(c)` and the style of its
blocks, `_card_style(c, T)` (e.g. the colors of the dark card or of the theme of the app).
Optionally, by dispatch on the host: `_card_value(c, v)` (e.g. of an `_AxisColor`),
`_host_attributes(c, T, attributes)`, `_row_attributes(c)`, `_declarations(c, obj)`,
`_fix_caret!(c, block)` and `_on_content_built!(gui, c)`.
"""
abstract type _AbstractCard end

"""
    _ComponentCard(fig::Figure)

Card of the live view with the controls of an object, shown over the 3D view next to the bounding
box of the object and connected to it by a line, see `_update_card!`. The card of the selected
object follows the selection; a `pinned` card stays with its object `obj`, independent of the
selection. Like the progress window (`_ProgressOverlay`), it has a dark (but opaque) `background`
with light text. It consists of free layouts (a `GridLayout` with a suggested bounding box) in
`scene`, a scene with a pixel camera over the whole figure:

- `head`: the `collapse_button` ("–" or "+"), the `pin_button` ("pin" or "unpin") and the `title`,
  i.e. the label of the object
- `actions`, right of the head: the buttons of [`card_actions`](@ref) for the object
- `rows`, below the head unless `collapsed`: the rows of [`card_rows`](@ref) for the object, each in
  its own layout
- `step`, below the rows, only on the card of the selection (the keyboard steps move the selected
  object): the `step_box` of the keyboard step

The widgets of the actions and the rows are built from their declarations when the card gets an
object with other declarations, see `_build_content!`: `widgets` holds each widget with its
[`CardWidget`](@ref), `blocks` all blocks (including the texts), `textboxes` the textboxes,
`listeners` their listeners and `content_key` the layout of the declarations. While `refreshing`,
the widgets show new values and their inputs are ignored.

The `link` holds the ends of the line, the object first; its plots, like all plots of the card, are
not inspectable. For a pinned card, `corners` are the corners of the bounding box of its object in
the pose `key = (obj, P, R)` when it was pinned; they move with the object, see `_card_corners`.
`pose` is the object and its pose that the widgets show.

Parts that are not shown are moved far outside of the figure (`_CARD_AWAY`), since hidden widgets
still take clicks within their bounding box. `scene` is translated to `z` (see `_card_z`), so that
GLMakie draws the card over the 3D scene, including plots added later and transparent plots.
"""
mutable struct _ComponentCard <: _AbstractCard
    scene::Scene
    background::Box
    head::GridLayout
    actions::GridLayout
    rows::GridLayout
    step::GridLayout
    title::Label
    collapse_button::Button
    pin_button::Button
    step_box::Textbox
    link::Observable{Vector{Point2f}}
    widgets::Vector{Tuple{Any, CardWidget}}
    blocks::Vector{Any}
    textboxes::Vector{Textbox}
    listeners::Vector{Any}
    content_key::Any
    refreshing::Bool
    collapsed::Bool
    auto_collapsed::Bool
    pinned::Bool
    obj::Any
    corners::Vector{Point3f}
    key::Any
    pose::Any
end

# Content of a layout of the card `scene`, aligned at the top left corner of its suggested bounding box
function _card_part(scene::Scene)
    layout = GridLayout(; bbox = _CARD_AWAY, halign = :left, valign = :top, default_colgap = 6,
        default_rowgap = 6)
    layout.parent = scene
    return layout
end

function _ComponentCard(fig::Figure, z::Real = _CARD_Z)
    scene = Scene(fig.scene; camera = Makie.campixel!, clear = false)
    translate!(scene, 0, 0, z)
    # The line first, then the background, then the widgets: each covers the one before
    link = Observable(Point2f[])
    lines!(scene, link; color = _CARD_LINK_COLOR, linewidth = 1.5, inspectable = false)
    scatter!(scene, Makie.lift(l -> l[1:min(1, end)], link); color = _CARD_LINK_COLOR, markersize = 7,
        strokecolor = :black, strokewidth = 1, inspectable = false)
    background = Box(scene; bbox = _CARD_AWAY, color = _CARD_BACKGROUND, strokevisible = false,
        cornerradius = _PROGRESS_CORNER)
    head, step = _card_part(scene), _card_part(scene)
    collapse_button = Button(head[1, 1]; label = "–", width = 22, height = 20, padding = (0, 0, 0, 0),
        fontsize = _CARD_FONTSIZE)
    pin_button = Button(head[1, 2]; label = "pin", height = 20, padding = (6, 6, 0, 0), fontsize = _CARD_FONTSIZE)
    title = Label(head[1, 3], ""; font = :bold, halign = :left, _card_style(Label)...)
    Label(step[1, 1], "step"; halign = :right, _card_style(Label)...)
    step_box = Textbox(step[1, 2]; placeholder = "e.g. 250 nm", width = 110, _card_style(Textbox)...)
    _fix_caret!(step_box, z)
    scene.visible[] = false
    return _ComponentCard(scene, background, head, _card_part(scene), _card_part(scene), step, title,
        collapse_button, pin_button, step_box, link, Tuple{Any, CardWidget}[], Any[], Textbox[], Any[],
        nothing, false, false, false, false, nothing, Point3f[], nothing, nothing)
end

"""Returns the textboxes of the card `c`: the step box and the declared ones."""
_card_boxes(c::_ComponentCard) = (c.step_box, c.textboxes...)

"""Returns the declared widget with the `name` on the card `c` (see [`CardWidget`](@ref)), or `nothing`."""
function _card_widget(c::_AbstractCard, name::Symbol)
    i = findfirst(((_, w),) -> w.name === name, c.widgets)
    return isnothing(i) ? nothing : first(c.widgets[i])
end

#=
Widgets of the declarations, by the type of the block
=#

# Colors and sizes of the card, which is dark and compact, for the blocks of the declarations. The
# box and border colors of the textboxes are Makie's defaults, which a theme of the figure (e.g. of
# the app layout, with light boxes) must not change
_card_style(::Type{Label}) = (; color = _PROGRESS_TEXT_COLOR, fontsize = _CARD_FONTSIZE)
_card_style(::Type{Textbox}) = (; textcolor = _PROGRESS_TEXT_COLOR, boxcolor = :transparent,
    boxcolor_hover = :transparent, boxcolor_focused = :transparent,
    bordercolor = RGBAf(0.8, 0.8, 0.8, 1), fontsize = _CARD_FONTSIZE, height = 26,
    textpadding = (6, 6, 4, 4))
_card_style(::Type{Button}) = (; fontsize = _CARD_FONTSIZE, height = 24, padding = (8, 8, 2, 2))
_card_style(::Type) = (;)
_card_style(::_ComponentCard, T::Type) = _card_style(T)

"""
    _AxisColor(k)

Color of the gizmo axis `k` (1: red, 2: green, 3: blue) in a declaration, e.g. of the labels of the
rotation boxes (see `pose_card_rows`), which each host of the card resolves to its own shade, see
`_card_value`: lighter on the dark floating card, the axis colors of the theme when docked.
"""
struct _AxisColor
    k::Int
end

# Values of the attributes of a declaration on the card `c`, see `_AxisColor`
_card_value(::_AbstractCard, v) = v
_card_value(::_ComponentCard, a::_AxisColor) = _CARD_AXIS_COLORS[a.k]

"""
    _cell_attributes(c, w::CardWidget) -> NamedTuple

Attributes of the block of the declared widget `w` on the card `c`: the style of the card for the
type of the block (see `_card_style`), then the attributes of `w` with their values on the card
(see `_card_value`), which the host may adapt, see `_host_attributes`.
"""
_cell_attributes(c::_AbstractCard, w::CardWidget) = _host_attributes(c, w.type,
    (; _card_style(c, w.type)..., map(v -> _card_value(c, v), w.attributes)...))

# The attributes of a block of the type `T` on the card, as declared by default
_host_attributes(::_AbstractCard, T::Type, attributes::NamedTuple) = attributes

# Layout of a declared row of the card, see `_build_content!`
_row_attributes(::_ComponentCard) = (; halign = :left, default_colgap = 6)

# The caret and the selection of a textbox are drawn without the translation of the scene of the
# card, which is at `z`
function _fix_caret!(tb::Textbox, z::Real)
    for p in tb.editor.plots
        p isa Makie.Text || translate!(p, 0, 0, z + 5)
    end
    return nothing
end
_fix_caret!(_, ::Real) = nothing
_fix_caret!(c::_ComponentCard, b) = _fix_caret!(b, _scene_z(c))
_fix_caret!(::_AbstractCard, _) = nothing

# The z translation of the scene of the card `c`
_scene_z(c::_ComponentCard) = c.scene.transformation.translation[][3]

# The observable of a block that carries its inputs, see `CardWidget`
_widget_input(b::Slider) = b.value
_widget_input(b::Toggle) = b.active
_widget_input(b::Textbox) = b.stored_string
_widget_input(b::Button) = b.clicks
_widget_input(b::Menu) = b.selection
_widget_input(_) = nothing

"""
    _show!(block, v; force = false)

Shows the value `v` of a declared widget (see [`CardWidget`](@ref)) in its `block`. A focused
textbox keeps the typed text, unless `force`.
"""
_show!(b::Label, v; force = false) = _update!(b.text, string(v))
_show!(b::Button, v; force = false) = _update!(b.label, string(v))
_show!(b::Textbox, v; force = false) = ((b.focused[] && !force) || _set_box!(b, string(v)); nothing)
_show!(b::Slider, v; force = false) = (b.value[] == v || Makie.set_close_to!(b, v); nothing)
_show!(b::Toggle, v; force = false) = _update!(b.active, Bool(v))
_show!(_, _; force = false) = nothing

# The layout of declarations, independent of the functions `value` and `on`, see `_build_content!`
_layout_key(rows) = map(_layout_key, rows)
_layout_key(r::CardRow) = map(_layout_key, r.cells)
_layout_key(s::String) = s
_layout_key(w::CardWidget) = (w.type, w.attributes, w.name)

# The declared widgets of rows or actions, in the order of their blocks
_declared_widgets(rows) = CardWidget[w for r in rows for w in _declared_widgets(r)]
_declared_widgets(r::CardRow) = CardWidget[w for c in r.cells for w in _declared_widgets(c)]
_declared_widgets(w::CardWidget) = (w,)
_declared_widgets(::String) = ()

"""Moves the parts of the card `c` away, see `_CARD_AWAY`."""
_park_card!(c::_ComponentCard) = foreach(_park!, (c.head, c.actions, c.rows, c.step, c.background))

"""Removes the declared widgets of the card `c`, with their listeners and layouts."""
function _clear_content!(c::_AbstractCard)
    foreach(off, c.listeners)
    foreach(delete!, c.blocks)
    empty!(c.listeners)
    empty!(c.blocks)
    empty!(c.widgets)
    empty!(c.textboxes)
    # New layouts instead of the empty rows and columns of the old ones
    _new_parts!(c)
    c.content_key = nothing
    return nothing
end

"""Replaces the layouts of the actions and the rows of the card `c` by new ones, see `_clear_content!`."""
function _new_parts!(c::_ComponentCard)
    c.actions, c.rows = _card_part(c.scene), _card_part(c.scene)
    return nothing
end

"""Returns the parts of the card `c` below its head that are shown, see `_ComponentCard`."""
function _lower_parts(c::_ComponentCard)
    (c.collapsed || c.auto_collapsed) && return GridLayout[]
    rows = isempty(c.rows.content) ? GridLayout[] : [c.rows]
    return c.pinned ? rows : [rows..., c.step]
end

# Size of a layout or block [px], which does not depend on its position
_card_size(x) = Vec2f(Makie.widths(x.layoutobservables.computedbbox[]))
# Moves a layout or block such that its top left corner is at `p`, or away, see `_CARD_AWAY`
_place!(x, p::Point2f) = _update!(x.layoutobservables.suggestedbbox, Rect2f(p[1], p[2], 0, 0))
_park!(x) = _update!(x.layoutobservables.suggestedbbox, _CARD_AWAY)

# Size of the actions, empty without actions
_actions_size(c::_ComponentCard) = isempty(c.actions.content) ? Vec2f(-_CARD_PADDING, 0) : _card_size(c.actions)

"""
    _card_size(c::_ComponentCard) -> Vec2f

Size of the card `c` [px] with its actions right of the head and the parts below it (see
`_lower_parts`), including the padding of the background.
"""
function _card_size(c::_ComponentCard)
    h, a = _card_size(c.head), _actions_size(c)
    w, height = h[1] + _CARD_PADDING + a[1], max(h[2], a[2])
    for part in _lower_parts(c)
        s = _card_size(part)
        w, height = max(w, s[1]), height + _CARD_PADDING + s[2]
    end
    return Vec2f(w, height) .+ 2 * _CARD_PADDING
end

"""
    _arrange_card!(c::_ComponentCard, p::Point2f)

Moves the card `c` such that its top left corner is at the figure pixel `p`, the parts that are not
shown (see `_lower_parts`) away. Only changed positions update the layout.
"""
function _arrange_card!(c::_ComponentCard, p::Point2f)
    size = _card_size(c)
    x, y = p[1] + _CARD_PADDING, p[2] - _CARD_PADDING
    h, a = _card_size(c.head), _actions_size(c)
    line = max(h[2], a[2])
    # The head and the actions are centered vertically in the first line
    _place!(c.head, Point2f(x, y - (line - h[2]) / 2))
    _place!(c.actions, Point2f(x + h[1] + _CARD_PADDING, y - (line - a[2]) / 2))
    y -= line
    lower = _lower_parts(c)
    for part in (c.rows, c.step)
        if any(l -> l === part, lower)
            y -= _CARD_PADDING
            _place!(part, Point2f(x, y))
            y -= _card_size(part)[2]
        else
            _park!(part)
        end
    end
    _update!(c.background.layoutobservables.suggestedbbox, Rect2f(p[1], p[2] - size[2], size...))
    return nothing
end

"""Ends the input into the textboxes of the card `c`, see `Makie.defocus!`."""
function _defocus_card!(c::_ComponentCard)
    for tb in _card_boxes(c)
        tb.focused[] && Makie.defocus!(tb)
    end
    return nothing
end

"""Hides the card `c`: ends the input into its textboxes and moves all its parts away."""
function _hide_card!(c::_ComponentCard)
    c.scene.visible[] || return nothing
    _defocus_card!(c)
    _park_card!(c)
    _update!(c.link, Point2f[])
    c.scene.visible[] = false
    return nothing
end

"""Returns `true` if the card `c` is shown and the mouse of the `events` is over it."""
function _over_card(c::_ComponentCard, events::Makie.Events)
    c.scene.visible[] || return false
    return Point2f(events.mouseposition[]) in c.background.layoutobservables.computedbbox[]
end

"""Returns the 8 corners of the box `bb`."""
function _box_corners(bb)
    lo, hi = Vector{Float64}(minimum(bb)), Vector{Float64}(maximum(bb))
    return vec([Point3f(x, y, z) for x in (lo[1], hi[1]), y in (lo[2], hi[2]), z in (lo[3], hi[3])])
end

"""
    _screen_rect(scene, pts, obj) -> Rect2f

Rectangle [figure px] around the projections of the 3D points `pts` (e.g. the corners of the
bounding box of the object `obj`) in the 3D `scene`. Points behind the camera are skipped. If no
point remains, e.g. for an object behind the camera, the rectangle is the point at the edge of the
view towards `obj`, see `_screen_anchor`.
"""
function _screen_rect(scene::Scene, pts, obj)
    o = Point2f(minimum(Makie.viewport(scene)[]))
    qs = Point2f[]
    for p in pts
        (all(isfinite, p) && !_behind(scene, p)) || continue
        q = Makie.project(scene, :data, :pixel, Point3f(p))
        all(isfinite, q) && push!(qs, o + Point2f(q[1], q[2]))
    end
    isempty(qs) && push!(qs, o + _screen_anchor(scene, position(obj)))
    lo = reduce((a, b) -> min.(a, b), qs)
    hi = reduce((a, b) -> max.(a, b), qs)
    return Rect2f(lo, hi - lo)
end

"""
    _link_anchor(scene, pts, obj) -> Point2f

End of the line of a card at its object [figure px]: the projection of the center of the points
`pts` (the corners of the bounding box of `obj`), or the point at the edge of the view towards it,
see `_screen_anchor`, if it is behind the camera.
"""
function _link_anchor(scene::Scene, pts, obj)
    p = isempty(pts) ? Point3f(position(obj)) :
        (reduce((a, b) -> min.(a, b), pts) + reduce((a, b) -> max.(a, b), pts)) / 2
    o = Point2f(minimum(Makie.viewport(scene)[]))
    q = Makie.project(scene, :data, :pixel, Point3f(p))
    (_behind(scene, p) || !all(isfinite, q)) && return o + _screen_anchor(scene, p)
    return o + Point2f(q[1], q[2])
end

"""
    _card_position(sel::Rect2f, size::Vec2f, view::Rect2f) -> Point2f

Top left corner [figure px] of a card of the `size` at the screen rectangle `sel` of the bounding
box of its object in the 3D view `view`: the first of right of `sel` and top-aligned with it, left
of it, below it and above it (left-aligned) that lies inside the view with the margin
`_CARD_MARGIN`. If none does, e.g. for an object that fills the view, the card is right of `sel`,
moved into the view, where it covers the object.
"""
function _card_position(sel::Rect2f, size::Vec2f, view::Rect2f)
    lo, hi = minimum(view) .+ _CARD_MARGIN, maximum(view) .- _CARD_MARGIN
    fits(q) = q[1] >= lo[1] && q[1] + size[1] <= hi[1] && q[2] - size[2] >= lo[2] && q[2] <= hi[2]
    top = clamp(maximum(sel)[2], min(hi[2], lo[2] + size[2]), hi[2])
    left = clamp(minimum(sel)[1], lo[1], max(lo[1], hi[1] - size[1]))
    candidates = (Point2f(maximum(sel)[1] + _CARD_GAP, top), Point2f(minimum(sel)[1] - _CARD_GAP - size[1], top),
        Point2f(left, minimum(sel)[2] - _CARD_GAP), Point2f(left, maximum(sel)[2] + _CARD_GAP + size[2]))
    for q in candidates
        fits(q) && return q
    end
    return Point2f(clamp(candidates[1][1], lo[1], max(lo[1], hi[1] - size[1])), top)
end

# Rectangle of a card with the top left corner `p` and the `size`
_card_rect(p::Point2f, size::Vec2f) = Rect2f(p[1], p[2] - size[2], size...)
_overlaps(a::Rect2f, b::Rect2f) = all(minimum(a) .< maximum(b)) && all(minimum(b) .< maximum(a))
# The card with the top left corner `p` and the `size` overlaps one of the `obstacles`
_covers(p::Point2f, size::Vec2f, obstacles) = any(o -> _overlaps(_card_rect(p, size), o), obstacles)

"""Returns the screen rectangles [figure px] that the cards keep off: the view cube, if any."""
_obstacles(::Nothing) = Rect2f[]
_obstacles(cube::ViewCube) = [Rect2f(Makie.viewport(cube.scene)[])]

"""
    _avoid(p::Point2f, size::Vec2f, view::Rect2f, obstacles) -> Point2f

Moves the top left corner `p` of a card of the `size` off the `obstacles`, i.e. the view cube and
the cards placed before (see `_update_cards!`), which would otherwise take its clicks or the other
way round: to the nearest position that lies inside the `view` with the margin `_CARD_MARGIN` and
overlaps no obstacle. The candidates combine the coordinates of `p`, of the edges of the view and of
the positions next to each obstacle (left or right of it, below or above it). Keeps `p` if none is
free, i.e. if the view is too full.
"""
function _avoid(p::Point2f, size::Vec2f, view::Rect2f, obstacles)
    free(q) = !_covers(q, size, obstacles)
    free(p) && return p
    lo, hi = minimum(view) .+ _CARD_MARGIN, maximum(view) .- _CARD_MARGIN
    m = _CARD_MARGIN
    # Left edges and top edges of the card
    xs, ys = Float32[p[1], lo[1], hi[1] - size[1]], Float32[p[2], hi[2], lo[2] + size[2]]
    for o in obstacles
        push!(xs, minimum(o)[1] - m - size[1], maximum(o)[1] + m)
        push!(ys, minimum(o)[2] - m, maximum(o)[2] + m + size[2])
    end
    inside(q) = q[1] >= lo[1] - 1.0f-3 && q[1] + size[1] <= hi[1] + 1.0f-3 &&
                q[2] - size[2] >= lo[2] - 1.0f-3 && q[2] <= hi[2] + 1.0f-3
    candidates = [Point2f(x, y) for x in xs for y in ys]
    filter!(q -> inside(q) && free(q), candidates)
    isempty(candidates) && return p
    return argmin(q -> norm(q - p), candidates)
end
