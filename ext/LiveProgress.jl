using Makie: Observable, Point2f, Vec2f, Rect

# Progress window of the live view, anchored at a point of the 3D scene (e.g. a source marker)

# Layout in pixels relative to the anchor: the panel sits above and to the right of the anchor,
# such that it does not cover the marker at the anchor
const _PROGRESS_GAP = 16.0f0
const _PROGRESS_PANEL = Vec2f(220, 46)
# Minimum distance of the panel from the edges of the view
const _PROGRESS_MARGIN = 8.0f0
# z translation of the plots: GLMakie draws the plots in the order of this value, the window comes
# after the 3D scene, and its depth (≈ 0.01) lies in front of it, including transparent plots, which
# are drawn over plots with `overdraw`. Within the clip range ±10000 of the pixel camera, among the
# levels of the component cards, see `_CARD_Z`
const _PROGRESS_Z = 9800.0f0
const _PROGRESS_PADDING = 10.0f0
const _PROGRESS_TRACK = Vec2f(_PROGRESS_PANEL[1] - 2 * _PROGRESS_PADDING, 6)
const _PROGRESS_CORNER = 6.0f0
# Orange of the source markers, see `_live_render_source!`; the other colors come from the tokens
# of the theme of the live view, like those of the component cards
const _PROGRESS_FILL_COLOR = RGBAf(Makie.to_color(:orange))

# Left edge of the bar and vertical centers of the bar and the label
const _PROGRESS_BAR_X = _PROGRESS_GAP + _PROGRESS_PADDING
const _PROGRESS_BAR_Y = _PROGRESS_GAP + 14
const _PROGRESS_LABEL_Y = _PROGRESS_GAP + _PROGRESS_PANEL[2] - 15

"""
    _rounded_rect_marker(w, h, r)

Returns a `BezierPath` of a rectangle with rounded corners of radius `r`, centered at the origin,
of width `w` and height `h`. As a scatter marker of `markersize = w` it spans `w × h` in the
`markerspace`, if `w = 1`.
"""
function _rounded_rect_marker(w::Real, h::Real, r::Real)
    x, y = Float32(w / 2), Float32(h / 2)
    r = Float32(r)
    arc(cx, cy, a0, a1) = Makie.EllipticalArc(Point2f(cx, cy), r, r, 0.0f0, Float32(a0), Float32(a1))
    return Makie.BezierPath([
        Makie.MoveTo(Point2f(-x + r, -y)),
        Makie.LineTo(Point2f(x - r, -y)),
        arc(x - r, -y + r, -π / 2, 0),
        Makie.LineTo(Point2f(x, y - r)),
        arc(x - r, y - r, 0, π / 2),
        Makie.LineTo(Point2f(-x + r, y)),
        arc(-x + r, y - r, π / 2, π),
        Makie.LineTo(Point2f(-x, -y + r)),
        arc(-x + r, -y + r, π, 3π / 2),
        Makie.ClosePath()
    ])
end

"""
    _ProgressOverlay(ax::LScene, theme::NamedTuple = _app_theme(:light))

A small progress window of the live view in the 3D scene of `ax`: a panel with a label and a
progress bar next to an anchor point, e.g. the marker of the source that is being traced, in the
color tokens `theme` of the live view (see `_APP_THEMES`) like the component cards: the panel in
the color of the sidebars with a border, the text, the track of the bar in the color of the
borders and the bar in the orange of the source markers. The
panel has a fixed size on the screen (about 220 × 46 px, above and to the right of the anchor)
and is drawn on top of the scene. It is drawn in a child scene of the 3D scene with a pixel
camera, at the projection of the anchor, which [`_show_progress!`](@ref) updates, e.g. every frame
while the camera moves; an anchor outside of the view puts the panel at the edge towards it, see
`_screen_anchor`. Its plots come last in the drawing order of GLMakie and in front of the 3D scene
(`_PROGRESS_Z`), so that neither plots added to the 3D scene later, e.g. a sky, nor transparent
plots cover it. As plots of a child scene, they do not enter the limits of the 3D scene.

The constructor adds all plots at once, hidden, [`_show_progress!`](@ref) and
[`_hide_progress!`](@ref) only update observables. The plots are not inspectable and not clipped.

# Fields

- `scene`: the scene of the `LScene`
- `hud`: the child scene with a pixel camera that holds the plots
- `anchor`: the pixel position in `scene` the panel is attached to, `NaN` until first shown
- `fill_size`: width and height of the bar fill in pixels
- `fill_offset`: pixel offset of the center of the bar fill from the anchor
- `label`: text above the bar, e.g. "Tracing beams 42 %"
- `visible`: visibility of all plots
- `plots`: panel, bar track, bar fill and label
"""
struct _ProgressOverlay
    scene::Makie.Scene
    hud::Makie.Scene
    anchor::Observable{Point2f}
    fill_size::Observable{Vec2f}
    fill_offset::Observable{Vec2f}
    label::Observable{String}
    visible::Observable{Bool}
    plots::Vector{Makie.AbstractPlot}
end

function _ProgressOverlay(ax::LScene, theme::NamedTuple = _app_theme(:light))
    t = theme
    # Pixel coordinates of the viewport of the 3D scene, see `_screen_anchor`
    hud = Makie.campixel(ax.scene)
    anchor = Observable(Point2f(NaN))
    fill_size = Observable(Vec2f(0, _PROGRESS_TRACK[2]))
    fill_offset = Observable(Vec2f(_PROGRESS_BAR_X, _PROGRESS_BAR_Y))
    label = Observable("")
    visible = Observable(false)
    common = (; markerspace = :pixel, visible, inspectable = false,
        transparency = false, clip_planes = Makie.Plane3f[])
    w, h = _PROGRESS_PANEL
    panel = scatter!(hud, anchor; common...,
        marker = _rounded_rect_marker(1, h / w, _PROGRESS_CORNER / w), markersize = w,
        color = t.sidebar, strokecolor = t.border, strokewidth = 1,
        marker_offset = Vec2f(_PROGRESS_GAP) + _PROGRESS_PANEL / 2)
    track = scatter!(hud, anchor; common..., marker = Rect, markersize = _PROGRESS_TRACK,
        color = t.border,
        marker_offset = Vec2f(_PROGRESS_BAR_X + _PROGRESS_TRACK[1] / 2, _PROGRESS_BAR_Y))
    bar = scatter!(hud, anchor; common..., marker = Rect, markersize = fill_size,
        color = _PROGRESS_FILL_COLOR, marker_offset = fill_offset)
    # Font size of the theme, like the status line of the live view
    text = text!(hud, anchor; common..., text = label, color = t.text,
        align = (:left, :center), offset = Vec2f(_PROGRESS_BAR_X, _PROGRESS_LABEL_Y))
    plots = Makie.AbstractPlot[panel, track, bar, text]
    # In this order, after all other plots
    for (k, p) in enumerate(plots)
        translate!(p, 0, 0, _PROGRESS_Z + k)
    end
    return _ProgressOverlay(ax.scene, hud, anchor, fill_size, fill_offset, label, visible, plots)
end

"""
    _screen_anchor(scene, p) -> Point2f

Pixel position in the `scene` at which the progress window of the 3D point `p` is attached: the
projection of `p`, moved into the view such that the whole window stays visible, with a margin of
`_PROGRESS_MARGIN`. A point outside of the view thus puts the window at the edge towards it. With
the perspective projection, the projection of a point behind the camera is mirrored at the center
of the view, hence it is mirrored back and moved beyond the edge first. A point without a position
(`NaN`) puts the window at the bottom edge.
"""
function _screen_anchor(scene::Makie.Scene, p)
    w, h = Makie.widths(Makie.viewport(scene)[])
    q = Makie.project(scene, :data, :pixel, Point3f(p))
    x = Vec2f(q[1], q[2])
    c = Vec2f(w / 2, h / 2)
    if _behind(scene, p) || !all(isfinite, x)
        d = c - x
        d = all(isfinite, d) && norm(d) > 0 ? d : Vec2f(0, -1)
        x = c + d * Float32(w + h) / norm(d)
    end
    lo = Vec2f(_PROGRESS_MARGIN - _PROGRESS_GAP)
    hi = max.(lo, Vec2f(w, h) - _PROGRESS_PANEL .- (_PROGRESS_GAP + _PROGRESS_MARGIN))
    return Point2f(clamp.(x, lo, hi))
end

"""
Returns `true` if the 3D point `p` lies behind the camera of the `scene` with the perspective
projection, whose projection is mirrored at the center of the view. Always `false` with the
orthographic projection.
"""
function _behind(scene::Makie.Scene, p)
    cam = Makie.camera(scene)
    return dot(Vec3f(p) - cam.eyeposition[], cam.view_direction[]) <= 0 &&
           Makie.cameracontrols(scene).settings.projectiontype[] == Makie.Perspective
end

# Called every frame: notify an observable only if its value changes
_update!(obs::Observable, value) = (obs[] == value || (obs[] = value); nothing)

"""
    _show_progress!(o::_ProgressOverlay, anchor::AbstractVector{<:Real}, fraction::Real, label::AbstractString)

Shows the progress window `o` at the 3D point `anchor` (see `_screen_anchor`) with the bar filled
to `clamp(fraction, 0, 1)` of its width and the text `label`. Meant to be called every frame while
a solve runs, which also moves the window with the camera: only changed values update the plots,
and no plots are added.
"""
function _show_progress!(o::_ProgressOverlay, anchor::AbstractVector{<:Real}, fraction::Real,
        label::AbstractString)
    width = Float32(clamp(fraction, 0, 1)) * _PROGRESS_TRACK[1]
    _update!(o.anchor, _screen_anchor(o.scene, anchor))
    # The fill is centered at its offset, shifting it by half its width keeps it left-aligned
    _update!(o.fill_offset, Vec2f(_PROGRESS_BAR_X + width / 2, _PROGRESS_BAR_Y))
    _update!(o.fill_size, Vec2f(width, _PROGRESS_TRACK[2]))
    _update!(o.label, String(label))
    _update!(o.visible, true)
    return nothing
end

"""
    _hide_progress!(o::_ProgressOverlay)

Hides the progress window `o`, see [`_show_progress!`](@ref).
"""
_hide_progress!(o::_ProgressOverlay) = _update!(o.visible, false)
