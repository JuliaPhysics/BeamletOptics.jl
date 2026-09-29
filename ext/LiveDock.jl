#=
Analysis dock of the app layout: a tab per panel, of which only the active one is shown, and the
lazy update of the detector panels, which computes only the panel of the active tab, see
`_DockTabs`
=#

using Makie: Aspect

# Geometry of the tab bar in pixels
const _TAB_HEIGHT = 30
const _TAB_PAD = 12.0f0            # horizontal padding of a tab
const _TAB_ICON = 16.0f0           # size of the icon in front of the title
const _TAB_ICON_GAP = 6.0f0        # gap between the icon and the title
const _TAB_UNDERLINE = 2.0f0       # height of the accent line under the active tab
const _TAB_ARROW = 22.0f0          # width of the scroll arrows, shown if the tabs do not fit

"""
    _TabBar(parent; kwargs...)

A flat tab bar like the document tabs of a CAD program, placed at a `GridPosition` `parent`
like a Makie block, with a fixed height of `_TAB_HEIGHT` pixels. Each tab shows an icon (see
`_icon`) and its title; the active tab has its title in `text_color`, its icon in `accent_color`
and an accent line below it, the other tabs are muted. A thin line in `border_color` separates
the bar from the content below it.

A left click on a tab sets `clicked` to its index (on the release over the same tab as the
press), the bar does not switch tabs by itself, see `_set_active!`. If the tabs do not fit, arrows
at both ends and the mouse wheel over the bar scroll by whole tabs.

The bar is drawn in a child scene with a pixel camera, with three plots (backgrounds and lines,
icons and arrows, titles), independent of the number of tabs. It updates its plots only on
`_add_tab!`, `_set_active!`, scrolling and resizing; it has no tick or mouse move listeners.

# Keyword arguments

- `background`, `text_color`, `muted_color`, `accent_color`, `border_color`: colors of the theme
- `fontsize = 13`, `font = :regular`

# Fields

- `clicked::Observable{Int}`: the index of the last clicked tab
- `titles`, `icons`: the titles and icon names of the tabs, see `_add_tab!`
- `active`: the index of the active tab, `0` for none
- `first`: the index of the first tab in view, changed by scrolling
"""
mutable struct _TabBar
    const box::Makie.Box
    const scene::Makie.Scene
    const clicked::Observable{Int}
    const titles::Vector{String}
    const icons::Vector{Symbol}
    active::Int
    first::Int
    # target of the last press of the left mouse button, see `_tab_hit`
    pressed::Any
    const plots::NamedTuple
    const listeners::Vector{Observables.ObserverFunction}
    const font::Any
    const fontsize::Float32
    const text_color::RGBAf
    const muted_color::RGBAf
    const accent_color::RGBAf
    const border_color::RGBAf
    # widths of the titles in pixels, measured once
    const widths::Dict{String, Float32}
end

function _TabBar(parent::Makie.GridPosition; background = RGBf(0.97, 0.97, 0.98),
        text_color = RGBf(0.12, 0.14, 0.16), muted_color = RGBf(0.41, 0.45, 0.49),
        accent_color = RGBf(0.18, 0.44, 0.86), border_color = RGBf(0.83, 0.85, 0.87),
        fontsize::Real = 13, font = :regular)
    box = Makie.Box(parent; visible = false, height = _TAB_HEIGHT, tellwidth = false)
    topscene = box.blockscene
    viewport = lift(Makie.round_to_IRect2D, topscene, box.layoutobservables.computedbbox)
    scene = Makie.Scene(topscene, viewport; camera = Makie.campixel!, clear = true,
        backgroundcolor = RGBAf(Makie.to_color(background)))
    common = (; inspectable = false, space = :pixel)
    lines = poly!(scene, [Rect2f(0, 0, 0, 0)]; common..., strokewidth = 0, color = RGBAf[_TRANSPARENT])
    # GLMakie can not draw an empty vector of markers, see `_tab_markers`
    markers = scatter!(scene, Point2f[]; common..., marker = _icon(:expand), markersize = _TAB_ICON,
        color = RGBAf[], rotation = Float32[], strokewidth = 0)
    labels = text!(scene, Point2f[]; common..., text = String[], color = RGBAf[], font, fontsize,
        align = (:left, :center), markerspace = :pixel)
    bar = _TabBar(box, scene, Observable(0), String[], Symbol[], 0, 1, nothing,
        (; lines, markers, labels), Observables.ObserverFunction[], _tree_font(scene, font),
        Float32(fontsize), RGBAf(Makie.to_color(text_color)), RGBAf(Makie.to_color(muted_color)),
        RGBAf(Makie.to_color(accent_color)), RGBAf(Makie.to_color(border_color)),
        Dict{String, Float32}())
    push!(bar.listeners, on(_ -> _redraw!(bar), viewport))
    ev = events(scene)
    # Before the camera and the kinematic controls of a 3D scene (priority 200)
    push!(bar.listeners, on(ev.scroll, priority = 300) do (dx, dy)
        _mouse_in(bar.scene) || return Consume(false)
        d = dy != 0 ? -sign(dy) : sign(dx)
        _scroll_tabs!(bar, round(Int, d))
        return Consume(true)
    end)
    push!(bar.listeners, on(ev.mousebutton, priority = 300) do event
        event.button == Mouse.left || return Consume(false)
        if event.action == Mouse.press
            _mouse_in(bar.scene) || return Consume(false)
            bar.pressed = _tab_hit(bar)
            return Consume(true)
        elseif event.action == Mouse.release && !isnothing(bar.pressed)
            pressed, bar.pressed = bar.pressed, nothing
            hit = _mouse_in(bar.scene) ? _tab_hit(bar) : nothing
            (isnothing(hit) || hit != pressed) && return Consume(true)
            if hit === :left || hit === :right
                _scroll_tabs!(bar, hit === :left ? -1 : 1)
            elseif hit isa Int
                bar.clicked[] = hit
            end
            return Consume(true)
        end
        return Consume(false)
    end)
    _redraw!(bar)
    return bar
end

"""Appends a tab with the `title` and the icon `icon` (see `_icon`) to the `bar`."""
function _add_tab!(bar::_TabBar, title::AbstractString, icon::Symbol)
    push!(bar.titles, String(title))
    push!(bar.icons, icon)
    _redraw!(bar)
    return bar
end

"""Marks the tab `i` of the `bar` as active and scrolls it into view."""
function _set_active!(bar::_TabBar, i::Int)
    bar.active = i
    if 1 <= i <= length(bar.titles)
        bar.first = min(bar.first, i)
        # Scroll right until the tab is in view
        while bar.first < i && all(s -> s.i != i, _tab_slots(bar).slots)
            bar.first += 1
        end
    end
    _redraw!(bar)
    return bar
end

_mouse_in(scene::Makie.Scene) = Point2f(events(scene).mouseposition[]) in Rect2f(scene.viewport[])

_tab_title_width(bar::_TabBar, s::String) =
    get!(() -> Float32(Makie.widths(Makie.text_bb(s, bar.font, bar.fontsize))[1]), bar.widths, s)

_tab_width(bar::_TabBar, i::Int) =
    2 * _TAB_PAD + _TAB_ICON + _TAB_ICON_GAP + _tab_title_width(bar, bar.titles[i])

"""
    _tab_slots(bar) -> (; slots, overflow)

Returns the tabs in view as `(; i, x0, x1)` (index and horizontal extent in pixels), from the tab
`bar.first` on, and whether the tabs do not all fit into the bar, i.e. the scroll arrows are shown.
"""
function _tab_slots(bar::_TabBar)
    w = Float32(Makie.widths(bar.scene.viewport[])[1])
    n = length(bar.titles)
    overflow = sum(i -> _tab_width(bar, i), 1:n; init = 0.0f0) > w
    overflow || (bar.first = 1)
    left, right = overflow ? (_TAB_ARROW, w - _TAB_ARROW) : (0.0f0, w)
    slots = @NamedTuple{i::Int, x0::Float32, x1::Float32}[]
    x = left
    for i in bar.first:n
        x1 = x + _tab_width(bar, i)
        # At least one tab, also if it is wider than the bar
        (x1 > right && !isempty(slots)) && break
        push!(slots, (; i, x0 = x, x1))
        x = x1
    end
    return (; slots, overflow)
end

"""The largest `first` of the `bar` at which the last tab is in view."""
function _last_first(bar::_TabBar)
    n = length(bar.titles)
    w = Float32(Makie.widths(bar.scene.viewport[])[1]) - 2 * _TAB_ARROW
    k, x = n, 0.0f0
    while k >= 1 && x + _tab_width(bar, k) <= w
        x += _tab_width(bar, k)
        k -= 1
    end
    return clamp(k + 1, 1, max(n, 1))
end

function _scroll_tabs!(bar::_TabBar, d::Int)
    f = clamp(bar.first + d, 1, _last_first(bar))
    f == bar.first && return nothing
    bar.first = f
    _redraw!(bar)
    return nothing
end

"""
Returns what is under the mouse in the `bar`: the index of a tab, `:left` or `:right` for the
scroll arrows, `:none` elsewhere in the bar.
"""
function _tab_hit(bar::_TabBar)
    vp = bar.scene.viewport[]
    x = Point2f(events(bar.scene).mouseposition[])[1] - Makie.origin(vp)[1]
    (; slots, overflow) = _tab_slots(bar)
    if overflow
        x < _TAB_ARROW && return :left
        x > Makie.widths(vp)[1] - _TAB_ARROW && return :right
    end
    for s in slots
        s.x0 <= x < s.x1 && return s.i
    end
    return :none
end

"""Passes the tabs in view to the plots of the `bar`."""
function _redraw!(bar::_TabBar)
    w, h = Makie.widths(bar.scene.viewport[])
    (; slots, overflow) = _tab_slots(bar)
    y = h / 2
    rects = Rect2f[Rect2f(0, 0, w, 1)]
    rect_colors = RGBAf[bar.border_color]
    label_pos, label_text, label_color = Point2f[], String[], RGBAf[]
    marker_pos, marker_shape, marker_color, marker_rot = Point2f[], Any[], RGBAf[], Float32[]
    for (; i, x0, x1) in slots
        active = i == bar.active
        push!(marker_pos, Point2f(x0 + _TAB_PAD + _TAB_ICON / 2, y))
        push!(marker_shape, _icon(bar.icons[i]))
        push!(marker_color, active ? bar.accent_color : bar.muted_color)
        push!(marker_rot, 0.0f0)
        push!(label_pos, Point2f(x0 + _TAB_PAD + _TAB_ICON + _TAB_ICON_GAP, y))
        push!(label_text, bar.titles[i])
        push!(label_color, active ? bar.text_color : bar.muted_color)
        if active
            push!(rects, Rect2f(x0, 0, x1 - x0, _TAB_UNDERLINE))
            push!(rect_colors, bar.accent_color)
        end
    end
    if overflow
        n = length(bar.titles)
        dim(c, on) = on ? c : RGBAf(c.r, c.g, c.b, 0.35f0 * c.alpha)
        # The icon `:expand` points to the right, rotated by π for the left arrow
        push!(marker_pos, Point2f(_TAB_ARROW / 2, y), Point2f(w - _TAB_ARROW / 2, y))
        push!(marker_shape, _icon(:expand), _icon(:expand))
        push!(marker_color, dim(bar.muted_color, bar.first > 1),
            dim(bar.muted_color, bar.first < _last_first(bar) && !isempty(slots) && slots[end].i < n))
        push!(marker_rot, Float32(π), 0.0f0)
    end
    p = bar.plots
    Makie.update!(p.lines; arg1 = rects, color = rect_colors)
    Makie.update!(p.markers; arg1 = marker_pos, marker = _tab_markers(marker_shape),
        color = marker_color, rotation = marker_rot)
    Makie.update!(p.labels; arg1 = label_pos, text = label_text, color = label_color)
    return nothing
end

"""
The markers of the icons of the tabs: one per tab, or a single marker without tabs, e.g. while the
dock has no panels, since GLMakie can not draw an empty vector of markers.
"""
_tab_markers(shapes::Vector{Any}) = isempty(shapes) ? _icon(:expand) : _tree_markers(shapes)

"""
    _DockTabs

State of the analysis dock of the `AppLayout`: the tab bar, per tab its content (a collapsible
part in the content cell of the dock, see `_set_shown!`) and its `DetectorPanel` (`nothing` for
other panels, see `_add_dock_panel!`), and the index of the `active` tab. Only the content of the
active tab is laid out and drawn.

# Lazy panels

Only the detector panel of the active tab is computed and drawn after a solve (see
`_computed_panels`, `_apply_panel!`), while the dock is shown. All other panels are `stale`: they
are computed once their tab becomes active or the dock is shown again (see `_refresh_dock!`), from
the hits of the last solve, if they are complete (`hits_valid`, i.e. not during or after a
cancelled solve). Panels with a history (`history = true`) are computed after every full solve to
record their metrics (see `_record_panel!`), but not drawn; their `fields` are kept, such that
showing them needs no second computation.
"""
mutable struct _DockTabs
    const bar::_TabBar
    const parts::Vector{_AppPart}
    const panels::Vector{Any}
    active::Int
    const stale::Base.IdSet{Any}
    const fields::IdDict{Any, Any}
    hits_valid::Bool
end

"""Returns the detector panel of the active tab of the dock `tabs`, `nothing` if it has none."""
_active_panel(tabs::_DockTabs) = tabs.active == 0 ? nothing : tabs.panels[tabs.active]

"""
    _build_dock!(layout::AppLayout, spec) -> Vector

Creates the tab bar of the analysis dock and a tab with a `DetectorPanel` per detector of
`spec.specs`, named by `spec.labels`, see `_build_layout`. The first tab is active. Returns the
panels.
"""
function _build_dock!(layout::AppLayout, spec)
    t = layout.theme
    g = layout.dock.grid
    bar = _TabBar(g[1, 1]; background = t.sidebar, text_color = t.text, muted_color = t.muted,
        accent_color = t.accent, border_color = t.border)
    layout.tabs = _DockTabs(bar, _AppPart[], Any[], 0, Base.IdSet{Any}(), IdDict{Any, Any}(), false)
    layout.dock_panels = Pair{String, GridLayout}[]
    rowgap!(g, 6)
    panels = Any[]
    for (i, (pd, mode, kw)) in enumerate(spec.specs)
        name = get(spec.labels, pd, "Detector $i")
        p = DetectorPanel(_add_dock_panel!(layout, name; icon = :detector)[1, 1], pd, name, mode, kw)
        _theme_panel!(p, t)
        _arrange_panel!(p)
        layout.tabs.panels[end] = p
        push!(layout.tabs.stale, p)
        push!(panels, p)
    end
    isempty(panels) || _select_tab!(layout, 1)
    return panels
end

"""
Sets the colors of the detector panel `p` that do not follow the theme of the figure to the
tokens `t` (see `_APP_THEMES`): the subtitle, the spot diagram and the power (or number of hits)
of the history in the text colors, and the x (red) and z (blue) lines of the centroid history
and of the profiles in the red and blue of `t.gizmo`. Used for the detector panels of all layouts.
"""
function _theme_panel!(p::DetectorPanel, t)
    p.ax.subtitlecolor[] = t.muted
    p.scatter_plot.color[] = t.text
    lines(ax) = filter(x -> x isa Makie.Lines, ax.scene.plots)
    x_color, z_color = t.gizmo[1], t.gizmo[3]
    if !isempty(p.history_axes)
        value_ax, centroid_ax = p.history_axes
        foreach(l -> l.color[] = t.text, lines(value_ax))
        cx, cz = lines(centroid_ax)
        cx.color[] = x_color
        cz.color[] = z_color
    end
    if !isnothing(p.profiles_ax)
        px, pz = lines(p.profiles_ax)
        px.color[] = x_color
        pz.color[] = z_color
    end
    return nothing
end

"""
Places the history and profiles axes of the panel `p` beside its axis instead of below it, since
the dock is wide and low: the axis in a column as wide as the dock is high, the history and the
profiles in columns that share the remaining width.
"""
function _arrange_panel!(p::DetectorPanel)
    grid = _GLB.gridcontent(p.ax).parent
    col = 1
    if !isempty(p.history_axes)
        col += 1
        for a in p.history_axes
            grid[1, col] = a
            a.height = nothing
        end
    end
    if !isnothing(p.profiles_ax)
        col += 1
        grid[1, col] = p.profiles_ax
        p.profiles_ax.height = nothing
    end
    Makie.trim!(grid)
    # The labels and ticks stay within the dock, i.e. under the tab bar
    grid.alignmode = Outside()
    colsize!(grid, 1, Aspect(1, 1.0))
    col > 1 && colgap!(grid, 24)
    return nothing
end

"""
    _add_dock_panel!(layout::AppLayout, title; icon = :panel_bottom) -> GridLayout

Appends a tab with the `title` (and the `icon`, see `_icon`) to the analysis dock and returns the
layout of its content. The new tab becomes the active one. The first panel switches the dock on,
later panels keep it collapsed if it is. Requires the toggle of the dock, see `_build_tools`.
"""
function _add_dock_panel!(layout::AppLayout, title::AbstractString; icon::Symbol = :panel_bottom)
    tabs = layout.tabs
    parent = layout.dock.grid
    # Placeholder box for `_set_shown!`, the content cell has the background of the dock. It is
    # created in the figure, since the layout of a collapsed dock is detached from it
    box = Box(layout.dock.box.parent; visible = false)
    parent[2, 1] = box
    part = _AppPart(parent, (2, 1), _ -> nothing, nothing, box, GridLayout(parent[2, 1]), true)
    push!(tabs.parts, part)
    push!(tabs.panels, nothing)
    push!(layout.dock_panels, String(title) => part.grid)
    _add_tab!(tabs.bar, title, icon)
    _select_tab!(layout, length(tabs.parts))
    active = layout.collapse.dock.active
    length(layout.dock_panels) == 1 && !active[] && (active[] = true)
    _update_dock!(layout)
    return part.grid
end
_add_dock_panel!(gui::AppView, title::AbstractString; kwargs...) =
    _add_dock_panel!(gui.layout, title; kwargs...)

"""
    _select_tab!(layout::AppLayout, i)

Shows the content of the tab `i` of the dock in place of the active one, which is removed from the
layout and hidden, see `_set_shown!`. Nothing is computed, see `_activate_tab!`.
"""
function _select_tab!(layout::AppLayout, i::Int)
    tabs = layout.tabs
    old = tabs.active
    tabs.active = i
    if old != i
        old == 0 || _set_shown!(tabs.parts[old], false)
        part = tabs.parts[i]
        _set_shown!(part, true)
        # The blocks of a collapsed dock stay hidden until it is shown, see `_update_dock!`
        layout.dock.shown || foreach(Makie.hide!, _blocks!(Any[], part.grid))
    end
    _set_active!(tabs.bar, i)
    return nothing
end

#=
Lazy panels, see `_DockTabs`
=#

_has_history(p::DetectorPanel) = !isempty(p.history_axes)

"""Whether the panel `p` is shown, i.e. its tab is active and the dock is not collapsed."""
_panel_shown(gui::AppView, p) = gui.layout.dock.shown && _active_panel(gui.layout.tabs) === p

_computed_panels(gui::AppView, preview::Bool) =
    filter(p -> _panel_shown(gui, p) || (!preview && _has_history(p)), gui.panels)

_shown_panels(gui::AppView) = filter(p -> _panel_shown(gui, p), gui.panels)

function _apply_panel!(gui::AppView, p::DetectorPanel, field; coarse, preview)
    tabs = gui.layout.tabs
    if _panel_shown(gui, p)
        _update_panel!(p, field; coarse, preview, record = !preview)
        delete!(tabs.stale, p)
    elseif !preview
        # Computed for its history only, see `_computed_panels`
        _record_panel!(p, field)
        coarse || (tabs.fields[p] = field)
    end
    return nothing
end

function _on_solve_started!(gui::AppView)
    tabs = gui.layout.tabs
    tabs.hits_valid = false
    foreach(p -> push!(tabs.stale, p), gui.panels)
    empty!(tabs.fields)
    return nothing
end

function _on_applied!(gui::AppView)
    gui.layout.tabs.hits_valid = true
    # e.g. the tab was switched while the solve ran in the background
    _refresh_dock!(gui)
    return nothing
end

"""
    _record_panel!(p::DetectorPanel, field)

Records the metrics of the `field` (see `_panel_field`) in the history of the panel `p` without
drawing anything, for a panel that is not shown, see `_DockTabs`. The metrics are the same as
those of `_update_panel!`: of the spot diagram, or of the computed intensity.
"""
function _record_panel!(p::DetectorPanel, field)
    try
        p.metrics = _field_metrics(p, field)
        isnothing(p.metrics) || _record_history!(p, p.metrics; draw = false)
    catch e
        p.last_error = _log_once(e, p.last_error, "update of the panel \"$(p.name)\"")
    end
    return nothing
end

_field_metrics(::DetectorPanel, ::Nothing) = nothing
_field_metrics(p::DetectorPanel, ::_SpotField) = _spot_metrics(BMO.spot_diagram(p.pd))
_field_metrics(::DetectorPanel, (x, z, I)::Tuple) = _intensity_metrics(x, z, I)
_field_metrics(::DetectorPanel, e::Exception) = throw(e)

"""
    _activate_tab!(gui::AppView, i)

Switches the dock to the tab `i` after a click on it: shows its content and computes its detector
panel if it is stale, see `_refresh_dock!`. A panel that shows a preview (coarse grid or preview
solve) becomes stale when it is left, since only shown panels are refined, see `_on_idle!`.
"""
function _activate_tab!(gui::AppView, i::Int)
    tabs = gui.layout.tabs
    (1 <= i <= length(tabs.parts) && i != tabs.active) || return nothing
    old = _active_panel(tabs)
    (gui.trace.coarse || gui.trace.preview) && !isnothing(old) && push!(tabs.stale, old)
    _select_tab!(gui.layout, i)
    _refresh_dock!(gui)
    return nothing
end

"""
    _refresh_dock!(gui::AppView)

Computes the detector panel of the active tab if it is stale and shown: from the field kept for its
history, if any, otherwise in a job like the refinement of a coarse preview (see `_start_job`), with
a progress window for long computations, cancelled by `Esc` or a change. Nothing is computed while
a solve runs or after a cancelled solve, the next solve updates the panel.
"""
function _refresh_dock!(gui::AppView)
    tabs = gui.layout.tabs
    p = _active_panel(tabs)
    (gui.layout.dock.shown && !isnothing(p) && p in tabs.stale) || return nothing
    _refresh_tab!(gui, p)
    return nothing
end

"""
    _refresh_tab!(gui::AppView, p)

Updates the stale panel `p` of the shown active tab, see `_refresh_dock!`: a `DetectorPanel` is
computed, the `update` of a `_UserPanel` of [`add_panel!`](@ref) is called, see `LiveCustom.jl`.
"""
function _refresh_tab!(gui::AppView, p::DetectorPanel)
    tabs = gui.layout.tabs
    if haskey(tabs.fields, p)
        _show_panel!(gui, p, pop!(tabs.fields, p))
        return nothing
    end
    (tabs.hits_valid && !_running(gui)) || return nothing
    msg = "computing the panel \"$(p.name)\", Esc cancels"
    status = gui.status.text[]
    # The job has the panel as its object, such that the status line is not set to "traced"
    job = _start_job(gui, r -> _show_panels!(gui, r, msg, status), p, empty(gui.pairs),
        empty(gui.beam_handles), Any[p]; coarse = gui.trace.coarse, timing = :panel_time)
    _run!(gui, job, msg)
    return nothing
end

"""
Marks the detector panel `p` stale after its options changed, see `_set_panel_options!`: a kept field
was computed with the old options and is dropped. The panel is computed once it is shown.
"""
function _refresh_panel!(gui::AppView, p::DetectorPanel)
    tabs = gui.layout.tabs
    delete!(tabs.fields, p)
    push!(tabs.stale, p)
    _refresh_dock!(gui)
    return nothing
end

"""Shows the `field` of the panel `p`, which the last solve left stale, see `_refresh_dock!`."""
function _show_panel!(gui::AppView, p::DetectorPanel, field)
    _update_panel!(p, field; coarse = gui.trace.coarse, preview = gui.trace.preview, record = false)
    _draw_history!(p)
    delete!(gui.layout.tabs.stale, p)
    return nothing
end

"""Shows the result `r` of a job of `_refresh_dock!`, the status line is restored to `status`."""
function _show_panels!(gui::AppView, r, msg, status)
    foreach((p, field) -> _show_panel!(gui, p, field), r.panels, r.fields)
    gui.status.text[] == msg && (gui.status.text[] = status)
    return nothing
end

"""Connects the tab bar and the collapse toggle of the dock of the `gui`, see `_DockTabs`."""
function _connect_dock!(gui::AppView)
    layout = gui.layout
    listeners = gui.controls.listeners
    push!(listeners, on(i -> _activate_tab!(gui, i), layout.tabs.bar.clicked))
    # After `_update_dock!`, which shows or collapses the dock
    push!(listeners, on(_ -> _refresh_dock!(gui), layout.collapse.dock.active))
    return nothing
end
