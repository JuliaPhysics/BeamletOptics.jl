#=
Beam inspection and measuring
=#

# Screen-space pick radius of the beams [px]
const _BEAM_PICK_RADIUS = 6.0

"""
    _BeamSegment

Rendered segment of a beam from `a` to `b` [m] along the `ray`. `l0` and `opl0` are the geometric
and optical path length from the source to `a` [m]. `beamlet` is the Gaussian beamlet whose chief
ray the segment belongs to, or `nothing`.
"""
struct _BeamSegment
    a::Vector{Float64}
    b::Vector{Float64}
    ray::BMO.AbstractRay
    l0::Float64
    opl0::Float64
    beamlet::Any
end

"""
    _push_beam_segments!(segs, rays, l0, opl0, flen; beamlet = nothing)

Appends the segments of the consecutive `rays` of one beam, starting at the path lengths `l0` and
`opl0` [m], like `_push_ray_segment!`: a ray without intersection has the length `flen`.
"""
function _push_beam_segments!(segs, rays, l0, opl0, flen; beamlet = nothing)
    l, opl = Float64(l0), Float64(opl0)
    for ray in rays
        isect = BMO.intersection(ray)
        len = isnothing(isect) ? Float64(flen) : Float64(length(isect))
        a = Vector{Float64}(position(ray))
        push!(segs, _BeamSegment(a, a .+ len .* Vector{Float64}(BMO.direction(ray)), ray, l, opl, beamlet))
        l += len
        opl += len * BMO.refractive_index(ray)
    end
    return segs
end

_parent_lengths(::Nothing) = (0.0, 0.0)
_parent_lengths(p) = (Float64(length(p)), Float64(BMO.optical_path_length(p)))

function _beam_segments!(segs, ray::BMO.AbstractRay; flen)
    return _push_beam_segments!(segs, (ray,), 0.0, 0.0, flen)
end

function _beam_segments!(segs, beam::Beam; flen)
    for child in PreOrderDFS(beam)
        _push_beam_segments!(segs, BMO.rays(child), _parent_lengths(child.parent)..., flen)
    end
    return segs
end

function _beam_segments!(segs, gauss::BMO.GaussianBeamlet; flen)
    # Along the chief ray, `gauss_parameters` of a child beamlet expects the length from the source
    for child in PreOrderDFS(gauss)
        _push_beam_segments!(segs, BMO.rays(child.chief), _parent_lengths(child.parent)..., flen;
            beamlet = child)
    end
    return segs
end

function _beam_segments!(segs, agb::BMO.AstigmaticGaussianBeamlet; flen)
    for child in PreOrderDFS(agb)
        _push_beam_segments!(segs, BMO.rays(child.c), _parent_lengths(child.parent)..., flen)
    end
    return segs
end

function _beam_segments!(segs, bg::BMO.AbstractBeamGroup; flen, render_every = 1)
    bms = BMO.beams(bg)
    for i in 1:render_every:length(bms)
        _beam_segments!(segs, bms[i]; flen)
    end
    return segs
end

"""Returns the rendered segments of the beam of the render handle `h`, see `_BeamSegment`."""
function _beam_segments(h::BeamRenderHandle)
    segs = _BeamSegment[]
    if h.thing isa BMO.AbstractBeamGroup
        return _beam_segments!(segs, h.thing; flen = h.flen, render_every = h.render_every)
    end
    return _beam_segments!(segs, h.thing; flen = h.flen)
end
_beam_segments(h::GaussianRenderHandle) = _beam_segments!(_BeamSegment[], h.thing; flen = h.flen)
function _beam_segments(h::AstigmaticGroupRenderHandle)
    return _beam_segments!(_BeamSegment[], h.thing; flen = h.flen, render_every = h.render_every)
end
_beam_segments(_) = _BeamSegment[]

"""Returns the distance of the point `p` from the segment `a`-`b` in 2D."""
function _point_segment_distance(p, a, b)
    ab = (b[1] - a[1], b[2] - a[2])
    L2 = ab[1]^2 + ab[2]^2
    t = L2 > 0 ? clamp(((p[1] - a[1]) * ab[1] + (p[2] - a[2]) * ab[2]) / L2, 0, 1) : 0.0
    return hypot(p[1] - a[1] - t * ab[1], p[2] - a[2] - t * ab[2])
end

"""
    _closest_on_segment(a, b, origin, dir)

Returns the parameter `s ∈ [0, 1]` of the point `a + s (b - a)` of the segment that is closest to
the line `origin + t dir`.
"""
function _closest_on_segment(a, b, origin, dir)
    u, w = b .- a, a .- origin
    A, B, C = dot(u, u), dot(u, dir), dot(dir, dir)
    D, E = dot(u, w), dot(dir, w)
    den = A * C - B^2
    s = den > 1e-12 * A * C ? (B * E - C * D) / den : 0.0
    return clamp(s, 0.0, 1.0)
end

"""
    _inspect_beam(gui)

Returns the point of the rendered beams of the `gui` under the cursor, i.e. on the segment whose
projection is closest to the cursor within `_BEAM_PICK_RADIUS` pixels, or `nothing`. The result is
`(; point, direction, length, opl, w, R)`: the point closest to the camera ray through the cursor
and the direction of its segment, the geometric and optical path length (Σ n·L) from the source
[m], and for Gaussian beamlets the radius `w` and the curvature `R` of `gauss_parameters` at the
point, otherwise `nothing`. While a solve runs in the background, which changes the beams, nothing
is inspected.
"""
function _inspect_beam(gui::LiveView)
    # The beams are being traced by a solve in the background
    _running(gui) && return nothing
    scene = gui.ax.scene
    cursor = _px(scene)
    origin, dir = _cursor_ray(scene)
    best, dmin = nothing, _BEAM_PICK_RADIUS
    for h in gui.beam_handles, seg in _beam_segments(h)
        # Segments behind the camera are not visible
        (dot(seg.a .- origin, dir) > 0 || dot(seg.b .- origin, dir) > 0) || continue
        pa = Makie.project(scene, :data, :pixel, Point3(seg.a))
        pb = Makie.project(scene, :data, :pixel, Point3(seg.b))
        d = _point_segment_distance(cursor, pa, pb)
        if d <= dmin
            best, dmin = seg, d
        end
    end
    isnothing(best) && return nothing
    s = _closest_on_segment(best.a, best.b, origin, dir)
    point = best.a .+ s .* (best.b .- best.a)
    L = norm(point .- best.a)
    len = best.l0 + L
    opl = best.opl0 + L * BMO.refractive_index(best.ray)
    w, R = if isnothing(best.beamlet)
        nothing, nothing
    else
        BMO.gauss_parameters(best.beamlet, len)[1:2]
    end
    return (; point, direction = Vector{Float64}(BMO.direction(best.ray)), length = len, opl, w, R)
end

"""Formats a vector with 4 decimal places."""
_direction_string(v) = "(" * join((string(round(x, digits = 4)) for x in v), ", ") * ")"

"""Formats the point `p` [m] in mm with 3 decimal places."""
_point_string(p) = "(" * join((_fmt3(1e3 * x) for x in p), ", ") * ") mm"

"""Describes the inspected point `info` of a beam, see `_inspect_beam`."""
function _inspection_string(info)
    s = "beam at $(_point_string(info.point)), direction $(_direction_string(info.direction)), " *
        "path $(_fmt3(1e3 * info.length)) mm, OPL $(_fmt3(1e3 * info.opl)) mm"
    isnothing(info.w) && return s
    # `R` is the curvature, shown as the radius of curvature
    r = iszero(info.R) ? "∞" : _signed_length_string(1 / info.R)
    return s * ", w = $(_length_string(info.w)), R = $r"
end

"""Returns a marker of the points `pts` in the 3D view of the `gui`, which is never clipped."""
function _point_marker!(gui::LiveView, pts; color = :magenta)
    return scatter!(gui.ax, pts; color, markersize = 10, strokecolor = _marker_stroke(gui.layout),
        strokewidth = 1, overdraw = true, clip_planes = Plane3f[])
end

"""Removes the marker and the result of the beam inspection of the `gui`, if any."""
function _clear_inspection!(gui::LiveView)
    gui.measure.inspection = nothing
    _release_info!(gui, _BeamPoint)
    isnothing(gui.measure.inspection_plot) && return nothing
    delete!(gui.ax, gui.measure.inspection_plot)
    gui.measure.inspection_plot = nothing
    return nothing
end

"""
Shows the inspected point `info` of a beam with a marker, in the status line and on a card at the
point in the 3D view (see `_show_info!`) of the `gui`.
"""
function _show_inspection!(gui::LiveView, info)
    _clear_inspection!(gui)
    gui.measure.inspection = info
    gui.measure.inspection_plot = _point_marker!(gui, [Point3f(info.point)])
    gui.status.text[] = _inspection_string(info)
    _show_info!(gui, _BeamPoint(info))
    return nothing
end

"""Removes the points, the result, the plots and the card of the measurement of the `gui`."""
function _clear_measurement!(gui::LiveView)
    empty!(gui.measure.points)
    gui.measure.result = nothing
    _release_info!(gui, _Measurement)
    foreach(p -> delete!(gui.ax, p), gui.measure.plots)
    empty!(gui.measure.plots)
    return nothing
end

"""
    _measure(a, b)

Returns `(; distance, angle)` between the measured points `a` and `b` (`(; point, obj)`): the
distance [m] and the angle between the optical axes (local y-axes) of the objects [rad], or
`nothing` unless both are components.
"""
function _measure(a, b)
    distance = norm(b.point .- a.point)
    angle = nothing
    if a.obj isa BMO.AbstractObject && b.obj isa BMO.AbstractObject
        na, nb = _pose(a.obj)[2][:, 2], _pose(b.obj)[2][:, 2]
        angle = acos(clamp(dot(na, nb) / (norm(na) * norm(nb)), -1, 1))
    end
    return (; distance, angle)
end

"""
    _add_measure_point!(gui, point, obj)

Adds the `point` [m] of the component `obj` (or `nothing` for a point of a beam) to the measurement
of the `gui`. The second point shows the distance, and the angle between two components, in the
status line with a dashed line between the points; a third point starts a new measurement.
"""
function _add_measure_point!(gui::LiveView, point, obj)
    length(gui.measure.points) >= 2 && _clear_measurement!(gui)
    foreach(p -> delete!(gui.ax, p), gui.measure.plots)
    empty!(gui.measure.plots)
    push!(gui.measure.points, (; point = Vector{Float64}(point), obj))
    pts = [Point3f(m.point) for m in gui.measure.points]
    name(m) = isnothing(m.obj) ? "beam" : _label(gui, m.obj)
    if length(pts) == 1
        gui.status.text[] = "measure: $(name(gui.measure.points[1])) at $(_point_string(point)), " *
                            "click the second point"
    else
        a, b = gui.measure.points
        gui.measure.result = _measure(a, b)
        s = "measure: $(name(a)) to $(name(b)): distance $(round(1e3 * gui.measure.result.distance, digits = 6)) mm"
        isnothing(gui.measure.result.angle) || (s *= ", angle $(_angle_string(gui.measure.result.angle))")
        gui.status.text[] = s
        push!(gui.measure.plots, lines!(gui.ax, pts; color = :magenta, linestyle = :dash,
            linewidth = 2, overdraw = true, clip_planes = Plane3f[]))
    end
    push!(gui.measure.plots, _point_marker!(gui, pts))
    _show_info!(gui, _Measurement(copy(gui.measure.points), gui.measure.result,
        [name(m) for m in gui.measure.points]))
    return nothing
end

"""Switches measuring of the `gui` on or off, which clears the measurement."""
function _set_measuring!(gui::LiveView, on::Bool)
    _clear_measurement!(gui)
    _clear_inspection!(gui)
    gui.status.text[] = on ? "measure: click two components or beams" : "measuring off"
    return nothing
end

"""
    _on_click!(gui, obj)

Called after a click in the 3D view with the selected object `obj`, or `nothing` for a click on no
component. While measuring, the position of the component or the point of the beam under the
cursor is added to the measurement. Otherwise a click on a beam inspects it, see `_inspect_beam`,
and a click elsewhere removes the inspection. Returns `true` if a beam was clicked, then the
selection is kept.
"""
function _on_click!(gui::LiveView, obj)
    obj isa LiveClipPlane && (obj = nothing)
    info = isnothing(obj) ? _inspect_beam(gui) : nothing
    if gui.widgets.measure_toggle.active[]
        if !isnothing(obj)
            _add_measure_point!(gui, position(obj), obj)
        elseif !isnothing(info)
            _add_measure_point!(gui, info.point, nothing)
        end
    elseif isnothing(info)
        _clear_inspection!(gui)
    else
        _show_inspection!(gui, info)
    end
    return !isnothing(info)
end

"""Connects the beam inspection, esc and the measure toggle of the `gui`."""
function _connect_inspection!(gui::LiveView)
    listeners = gui.controls.listeners
    gui.controls.on_click = function (obj)
        try
            return _on_click!(gui, obj)
        catch e
            gui.last_error = _log_once(e, gui.last_error, "beam inspection")
            return false
        end
    end
    # Before the controls, which consume esc to deselect, the key is passed on
    push!(listeners, on(events(gui.ax.scene).keyboardbutton, priority = 201) do event
        (event.action == Keyboard.press && event.key == Keyboard.escape) || return Consume(false)
        gui.controls.ignore_keys() && return Consume(false)
        _clear_inspection!(gui)
        _clear_measurement!(gui)
        return Consume(false)
    end)
    push!(listeners, on(v -> _set_measuring!(gui, v), gui.widgets.measure_toggle.active))
    return nothing
end

