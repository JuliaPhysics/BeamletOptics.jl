using Makie: Figure, Observable, Point2f, Vec2f, Rect2f, GridLayout, Textbox, Button, Label, Box

# Component cards of the live view: the controls of an object, next to it in the 3D view

# z translation of the scene of a card: GLMakie draws the plots in the order of this value, the
# card comes after the 3D scene and the progress window (`_PROGRESS_Z`), and its depth (≈ 0.005)
# lies in front of them. Plus the offsets of the widgets, it stays within the clip range ±10000 of
# the pixel camera
const _CARD_Z = 9900.0f0
# Distance of the card from the bounding box of its object [px]
const _CARD_GAP = 12.0f0
# Minimum distance of the card from the edges of the 3D view and from the view cube [px]
const _CARD_MARGIN = 8.0f0
# Distance of the widgets from the edge of the background and between the parts of the card [px]
const _CARD_PADDING = 8.0f0
# Suggested bounding box of the parts that are not shown: hidden widgets still take clicks within
# their bounding box, hence they are moved far outside of the figure
const _CARD_AWAY = Rect2f(-1.0f5, -1.0f5, 0, 0)
# Labels of the pose boxes, see `_POSE_FIELDS`, the units are given at the end of the rows
const _CARD_POSE_LABELS = ("x", "y", "z", "rx", "ry", "rv")
# Opaque, unlike the progress window, such that objects behind the widgets do not shine through
const _CARD_BACKGROUND = RGBAf(0.1, 0.1, 0.12, 1)
# Colors of the rotation labels: red, green and blue like the gizmo axes, lighter on the dark card
const _CARD_AXIS_COLORS = (RGBAf(1, 0.45, 0.45, 1), RGBAf(0.45, 0.85, 0.45, 1), RGBAf(0.55, 0.7, 1, 1))
# Line from the card to its object and the dot at the object
const _CARD_LINK_COLOR = RGBAf(0.95, 0.95, 0.95, 0.85)

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
- the actions of the object, right of the head, see `_card_actions`: `object_actions` holds the
  `hide_button` for components and sources, `plane_actions` the `flip_button` and the
  `remove_button` for clip planes
- `body`, below the head unless `collapsed`: the `pose_boxes` x, y, z [mm] and rx, ry, rv [mrad]
- `step`, below the body, only on the card of the selection (the keyboard steps move the selected
  object): the `step_box` of the keyboard step

The `link` holds the ends of the line, the object first; its plots, like all plots of the card, are
not inspectable. For a pinned card, `corners` are the corners of the bounding box of its object in
the pose `key = (obj, P, R)` when it was pinned; they move with the object, see `_card_corners`.
`pose` is the pose of the object that the pose boxes show.

The constructor adds all widgets, hidden; afterwards the card is only moved, see
`_arrange_card!`. Parts that are not shown are moved far outside of the figure (`_CARD_AWAY`),
since hidden widgets still take clicks within their bounding box. `scene` is translated to
`_CARD_Z`, so that GLMakie draws the card over the 3D scene, including plots added later and
transparent plots.
"""
mutable struct _ComponentCard
    scene::Scene
    background::Box
    head::GridLayout
    body::GridLayout
    step::GridLayout
    object_actions::GridLayout
    plane_actions::GridLayout
    title::Label
    collapse_button::Button
    pin_button::Button
    pose_boxes::Vector{Textbox}
    step_box::Textbox
    hide_button::Button
    flip_button::Button
    remove_button::Button
    link::Observable{Vector{Point2f}}
    collapsed::Bool
    pinned::Bool
    obj::Any
    corners::Vector{Point3f}
    key::Any
    pose::Any
end

function _ComponentCard(fig::Figure)
    scene = Scene(fig.scene; camera = Makie.campixel!, clear = false)
    translate!(scene, 0, 0, _CARD_Z)
    # The line first, then the background, then the widgets: each covers the one before
    link = Observable(Point2f[])
    lines!(scene, link; color = _CARD_LINK_COLOR, linewidth = 1.5, inspectable = false)
    scatter!(scene, Makie.lift(l -> l[1:min(1, end)], link); color = _CARD_LINK_COLOR, markersize = 7,
        strokecolor = :black, strokewidth = 1, inspectable = false)
    background = Box(scene; bbox = _CARD_AWAY, color = _CARD_BACKGROUND, strokevisible = false,
        cornerradius = _PROGRESS_CORNER)
    # Content of the layout, aligned at the top left corner of its suggested bounding box
    function part()
        layout = GridLayout(; bbox = _CARD_AWAY, halign = :left, valign = :top, default_colgap = 6,
            default_rowgap = 6)
        layout.parent = scene
        return layout
    end
    head, body, step, object_actions, plane_actions = part(), part(), part(), part(), part()
    text = (; color = _PROGRESS_TEXT_COLOR)
    collapse_button = Button(head[1, 1]; label = "–", width = 24, height = 22, padding = (0, 0, 0, 0))
    pin_button = Button(head[1, 2]; label = "pin", height = 22, padding = (6, 6, 0, 0))
    title = Label(head[1, 3], ""; font = :bold, halign = :left, text...)
    pose_boxes = Textbox[]
    for (k, name) in enumerate(_CARD_POSE_LABELS)
        row, col = divrem(k - 1, 3)
        color = k <= 3 ? _PROGRESS_TEXT_COLOR : _CARD_AXIS_COLORS[k - 3]
        Label(body[row + 1, 2col + 1], name; color, halign = :right)
        # Wide enough for e.g. -6869.709 (mm of a telescope), longer values scroll while typing
        push!(pose_boxes, Textbox(body[row + 1, 2col + 2]; placeholder = k <= 3 ? " " : "0", width = 80,
            textcolor = _PROGRESS_TEXT_COLOR))
    end
    Label(body[1, 7], "mm"; halign = :left, text...)
    Label(body[2, 7], "mrad"; halign = :left, text...)
    Label(step[1, 1], "step"; halign = :right, text...)
    step_box = Textbox(step[1, 2]; placeholder = "e.g. 250 nm", width = 110, textcolor = _PROGRESS_TEXT_COLOR)
    hide_button = Button(object_actions[1, 1]; label = "hide")
    flip_button = Button(plane_actions[1, 1]; label = "flip")
    remove_button = Button(plane_actions[1, 2]; label = "remove")
    # The caret and the selection of a textbox are drawn without the translation of the scene
    for tb in (pose_boxes..., step_box), p in tb.editor.plots
        p isa Makie.Text || translate!(p, 0, 0, _CARD_Z + 5)
    end
    scene.visible[] = false
    return _ComponentCard(scene, background, head, body, step, object_actions, plane_actions, title,
        collapse_button, pin_button, pose_boxes, step_box, hide_button, flip_button, remove_button,
        link, false, false, nothing, Point3f[], nothing, nothing)
end

"""Returns the actions of the card `c` for its object `obj`, see `_ComponentCard`."""
_card_actions(c::_ComponentCard, ::LiveClipPlane) = c.plane_actions
_card_actions(c::_ComponentCard, _) = c.object_actions

"""Returns the textboxes of the card `c`."""
_card_boxes(c::_ComponentCard) = (c.pose_boxes..., c.step_box)

"""Returns the parts of the card `c` below its head that are shown, see `_ComponentCard`."""
_lower_parts(c::_ComponentCard) = c.collapsed ? GridLayout[] : c.pinned ? [c.body] : [c.body, c.step]

# Size of a layout or block [px], which does not depend on its position
_card_size(x) = Vec2f(Makie.widths(x.layoutobservables.computedbbox[]))
# Moves a layout or block such that its top left corner is at `p`, or away, see `_CARD_AWAY`
_place!(x, p::Point2f) = _update!(x.layoutobservables.suggestedbbox, Rect2f(p[1], p[2], 0, 0))
_park!(x) = _update!(x.layoutobservables.suggestedbbox, _CARD_AWAY)

"""
    _card_size(c::_ComponentCard, actions) -> Vec2f

Size of the card `c` [px] with the `actions` right of the head and the parts below it (see
`_lower_parts`), including the padding of the background.
"""
function _card_size(c::_ComponentCard, actions::GridLayout)
    h, a = _card_size(c.head), _card_size(actions)
    w, height = h[1] + _CARD_PADDING + a[1], max(h[2], a[2])
    for part in _lower_parts(c)
        s = _card_size(part)
        w, height = max(w, s[1]), height + _CARD_PADDING + s[2]
    end
    return Vec2f(w, height) .+ 2 * _CARD_PADDING
end

"""
    _arrange_card!(c::_ComponentCard, actions, p::Point2f)

Moves the card `c` with the `actions` such that its top left corner is at the figure pixel `p`,
the other actions and the parts that are not shown (see `_lower_parts`) away. Only changed
positions update the layout.
"""
function _arrange_card!(c::_ComponentCard, actions::GridLayout, p::Point2f)
    size = _card_size(c, actions)
    x, y = p[1] + _CARD_PADDING, p[2] - _CARD_PADDING
    h, a = _card_size(c.head), _card_size(actions)
    line = max(h[2], a[2])
    # The head and the actions are centered vertically in the first line
    _place!(c.head, Point2f(x, y - (line - h[2]) / 2))
    _place!(actions, Point2f(x + h[1] + _CARD_PADDING, y - (line - a[2]) / 2))
    for other in (c.object_actions, c.plane_actions)
        other === actions || _park!(other)
    end
    y -= line
    lower = _lower_parts(c)
    for part in (c.body, c.step)
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
    foreach(_park!, (c.head, c.body, c.step, c.object_actions, c.plane_actions, c.background))
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

"""Returns the screen rectangles [figure px] that the cards keep off: the view cube, if any."""
_obstacles(::Nothing) = Rect2f[]
_obstacles(cube::ViewCube) = [Rect2f(Makie.viewport(cube.scene)[])]

"""
    _avoid(p::Point2f, size::Vec2f, view::Rect2f, obstacles) -> Point2f

Moves the top left corner `p` of a card of the `size` off the `obstacles`, i.e. the view cube and
the cards placed before (see `_update_cards!`), which would otherwise take its clicks or the other
way round: to the nearest position below, above, left or right of an obstacle that lies inside the
`view` with the margin `_CARD_MARGIN` and overlaps no obstacle. Keeps `p` if there is none.
"""
function _avoid(p::Point2f, size::Vec2f, view::Rect2f, obstacles)
    free(q) = !any(o -> _overlaps(_card_rect(q, size), o), obstacles)
    free(p) && return p
    lo, hi = minimum(view) .+ _CARD_MARGIN, maximum(view) .- _CARD_MARGIN
    m = _CARD_MARGIN
    candidates = Point2f[]
    for o in obstacles
        push!(candidates, Point2f(p[1], minimum(o)[2] - m), Point2f(p[1], maximum(o)[2] + m + size[2]),
            Point2f(minimum(o)[1] - m - size[1], p[2]), Point2f(maximum(o)[1] + m, p[2]))
    end
    inside(q) = q[1] >= lo[1] && q[1] + size[1] <= hi[1] && q[2] - size[2] >= lo[2] && q[2] <= hi[2]
    filter!(q -> inside(q) && free(q), candidates)
    isempty(candidates) && return p
    return argmin(q -> norm(q - p), candidates)
end
