#=
App layout of `live_view(...; layout = :app)`: toolbar, sidebars, analysis dock and status bar
around the 3D view, see `AppLayout`
=#

using Makie: Box, Fixed, Auto, Outside, Rect2f, Point2f, rowsize!

const _GLB = Makie.GridLayoutBase

"""
Color tokens of the app layout by `theme`: `background` (toolbar, status bar), `sidebar`, `view`
(background of the 3D view and the detector panels), `border`, `text`, `muted` (secondary text),
`accent` and `accent_soft` (active toggles), `field` (buttons and textboxes), `hover` and
`gizmo` (labels of the rotations about the red, green and blue axes of the controls).
"""
const _APP_THEMES = Dict{Symbol, NamedTuple}(
    :light => (;
        background = Makie.to_color("#eceef1"), sidebar = Makie.to_color("#f7f8fa"),
        view = Makie.to_color("#ffffff"), border = Makie.to_color("#d3d8de"),
        text = Makie.to_color("#1f2328"), muted = Makie.to_color("#69727d"),
        accent = Makie.to_color("#2f6fdb"), accent_soft = Makie.to_color("#d7e4fa"),
        field = Makie.to_color("#ffffff"), hover = Makie.to_color("#e2e6eb"),
        gizmo = Makie.to_color.((:red, :green, :blue))),
    :dark => (;
        background = Makie.to_color("#1e2023"), sidebar = Makie.to_color("#25282c"),
        view = Makie.to_color("#34373b"), border = Makie.to_color("#3d4147"),
        text = Makie.to_color("#e3e5e8"), muted = Makie.to_color("#9aa1a9"),
        accent = Makie.to_color("#6ea4ff"), accent_soft = Makie.to_color("#2b4166"),
        field = Makie.to_color("#1a1c1f"), hover = Makie.to_color("#33373c"),
        gizmo = Makie.to_color.(("#ff6b6b", "#5ccf5c", "#6ea4ff"))))

"""Returns the color tokens of the `theme`, see `_APP_THEMES`."""
function _app_theme(theme::Symbol)
    haskey(_APP_THEMES, theme) ||
        throw(ArgumentError("theme must be one of $(Tuple(sort!(collect(keys(_APP_THEMES))))), got :$theme"))
    return _APP_THEMES[theme]
end

"""Returns the Makie theme of the widgets of the app layout with the color tokens `t`."""
function _makie_theme(t)
    return (;
        Label = (; color = t.text, fontsize = 13),
        Box = (; color = t.sidebar, strokecolor = t.border, strokewidth = 1),
        Button = (; buttoncolor = t.field, buttoncolor_hover = t.hover,
            buttoncolor_active = t.accent_soft, labelcolor = t.text, labelcolor_hover = t.text,
            labelcolor_active = t.text, strokecolor = t.border, strokewidth = 1, cornerradius = 4,
            fontsize = 13, padding = (8, 8, 5, 5)),
        Textbox = (; boxcolor = t.field, boxcolor_hover = t.field, boxcolor_focused = t.field,
            bordercolor = t.border, bordercolor_hover = t.accent, bordercolor_focused = t.accent,
            textcolor = t.text, textcolor_placeholder = t.muted, cursorcolor = t.accent,
            fontsize = 13, textpadding = (6, 6, 5, 5), cornerradius = 4),
        Menu = (; cell_color_inactive_even = t.field, cell_color_inactive_odd = t.field,
            cell_color_hover = t.hover, cell_color_active = t.accent_soft,
            selection_cell_color_inactive = t.field, textcolor = t.text, textcolor_hover = t.text,
            textcolor_active = t.text, dropdown_arrow_color = t.muted, fontsize = 13,
            textpadding = (8, 8, 5, 5)),
        Slider = (; color_inactive = t.border, color_active = t.accent,
            color_active_dimmed = t.accent_soft, linewidth = 6),
        Axis = (; backgroundcolor = t.view, titlecolor = t.text, subtitlecolor = t.muted,
            xlabelcolor = t.text, ylabelcolor = t.text, xticklabelcolor = t.muted,
            yticklabelcolor = t.muted, xtickcolor = t.border, ytickcolor = t.border,
            xgridcolor = (t.muted, 0.15), ygridcolor = (t.muted, 0.15), bottomspinecolor = t.border,
            topspinecolor = t.border, leftspinecolor = t.border, rightspinecolor = t.border,
            titlesize = 14))
end

"""
    _AppPart

Collapsible part of the app layout (a sidebar or the dock): a background `box` and the `grid` of
its content, both placed at `pos` of the `parent` layout. `resize(size)` sets the size of its
column or row in the `parent`, which is `size` while the part is shown, see `_set_shown!`.
"""
mutable struct _AppPart
    parent::GridLayout
    pos::Tuple{Int, Int}
    resize::Function
    size::Any
    box::Box
    grid::GridLayout
    shown::Bool
end

"""
    _TextToggle

Placeholder of a toolbar toggle: a text `Button` whose clicks switch `active`, shown in the accent
colors while active. Replaced by the icon toggles in S3.
"""
struct _TextToggle
    button::Button
    active::Observable{Bool}
end

function _TextToggle(pos; label::String, active::Bool, theme)
    button = Button(pos; label)
    obs = Observable(active)
    on(_ -> (obs[] = !obs[]), button.clicks)
    on(obs; update = true) do a
        button.buttoncolor[] = a ? theme.accent_soft : theme.field
        button.labelcolor[] = a ? theme.accent : theme.text
        button.labelcolor_hover[] = button.labelcolor[]
        return nothing
    end
    return _TextToggle(button, obs)
end

"""
    AppLayout

Layout of `live_view(...; layout = :app)`, an application window around the 3D view:

- a toolbar at the top, an ordered list of groups of entries, see `_add_toolbar_entry!`
- the left sidebar, a stack of titled sections ("OBJECTS", "PARAMETERS" with the sliders), and the
  right sidebar ("PROPERTIES": selection, pose boxes, step, mode), see `_add_sidebar_section!`
- the analysis dock below the 3D view with the detector panels, see `_add_dock_panel!`
- the status bar with the status line and an info label (last trace, number of rays, projection)

The sidebars and the dock are collapsed via toggles in the toolbar, then the 3D view takes their
space, see `_set_shown!`. The colors come from the tokens `theme` of the `theme` kwarg, see
`_APP_THEMES`. The fields are set by `_build_layout`.
"""
mutable struct AppLayout <: AbstractLiveLayout
    theme::NamedTuple
    # the toolbar groups `name => grid`, each entry is a column of its grid
    toolbar::GridLayout
    groups::Vector{Pair{Symbol, GridLayout}}
    # collapsible parts and the sections `title => content` of the sidebars (`:left`, `:right`)
    left::_AppPart
    right::_AppPart
    dock::_AppPart
    sections::Dict{Symbol, Vector{Pair{String, GridLayout}}}
    # the stack of a sidebar grows its sections marked `grow`, otherwise the filler at the end
    fillers::Dict{Symbol, Label}
    growing::Dict{Symbol, Bool}
    dock_panels::Vector{Pair{String, GridLayout}}
    # toolbar entries without a field of `LiveView`
    collapse::@NamedTuple{left::_TextToggle, right::_TextToggle, dock::_TextToggle}
    clip_toggle::_TextToggle
    fit_button::Button
    help_button::Button
    # inspector labels and the info label of the status bar
    selection_label::Label
    type_label::Label
    mode_label::Label
    info::Label
    AppLayout(theme::NamedTuple) = new(theme)
end

"""`LiveView` with the app layout, i.e. `live_view(...; layout = :app)`."""
const AppView = LiveView{AppLayout}

_default_size(::AppLayout) = (1600, 950)

_figure(layout::AppLayout, size) =
    Figure(; size, backgroundcolor = layout.theme.background, figure_padding = 0,
        _makie_theme(layout.theme)...)

#=
Slots
=#

"""
    _add_toolbar_entry!(layout::AppLayout, group::Symbol) -> GridPosition

Returns the position of a new entry at the end of the toolbar `group`, e.g. for a `Button`. A new
group is appended to the toolbar, after a separator.
"""
function _add_toolbar_entry!(layout::AppLayout, group::Symbol)
    i = findfirst(g -> g.first == group, layout.groups)
    grid = if isnothing(i)
        k = length(layout.groups)
        if k > 0
            Box(layout.toolbar[1, 2k]; width = 1, height = 22, color = layout.theme.border,
                strokewidth = 0)
        end
        g = GridLayout(layout.toolbar[1, 2k + 1]; default_colgap = 4)
        push!(layout.groups, group => g)
        g
    else
        layout.groups[i].second
    end
    return grid[1, length(grid.content) + 1]
end
_add_toolbar_entry!(gui::AppView, group::Symbol) = _add_toolbar_entry!(gui.layout, group)

"""
    _add_sidebar_section!(layout::AppLayout, side::Symbol, title; grow = false) -> GridLayout

Appends a section with the `title` to the left (`side = :left`) or right (`:right`) sidebar and
returns the layout of its content. A section with `grow` takes the free height of the sidebar.
"""
function _add_sidebar_section!(layout::AppLayout, side::Symbol, title::AbstractString;
        grow::Bool = false)
    haskey(layout.sections, side) || throw(ArgumentError("side must be :left or :right, got :$side"))
    part = getfield(layout, side)
    stack = part.grid
    sections = layout.sections[side]
    n = length(sections)
    t = layout.theme
    Label(stack[2n + 1, 1], uppercase(title); halign = :left, font = :bold, fontsize = 11,
        color = t.muted, tellwidth = false)
    content = GridLayout(stack[2n + 2, 1])
    push!(sections, String(title) => content)
    # The filler moves below the new section, the rows keep their sizes
    stack[2n + 3, 1] = layout.fillers[side]
    rowsize!(stack, 2n + 1, Auto())
    grow && (layout.growing[side] = true)
    rowsize!(stack, 2n + 2, grow ? Auto(false) : Auto())
    rowsize!(stack, 2n + 3, layout.growing[side] ? Fixed(0) : Auto(false))
    return content
end
_add_sidebar_section!(gui::AppView, side::Symbol, title::AbstractString; kwargs...) =
    _add_sidebar_section!(gui.layout, side, title; kwargs...)

"""
    _add_dock_panel!(layout::AppLayout, title) -> GridLayout

Appends a panel with the `title` to the analysis dock and returns its layout. The first panel
switches the dock on, later panels keep it collapsed if it is. The panels are placed in a row (tabs
follow in S4b). Requires the toolbar, see `_build_toolbar`.
"""
function _add_dock_panel!(layout::AppLayout, title::AbstractString)
    grid = GridLayout(layout.dock.grid[1, length(layout.dock_panels) + 1])
    push!(layout.dock_panels, String(title) => grid)
    active = layout.collapse.dock.active
    length(layout.dock_panels) == 1 && !active[] && (active[] = true)
    _update_dock!(layout)
    return grid
end
_add_dock_panel!(gui::AppView, title::AbstractString) = _add_dock_panel!(gui.layout, title)

#=
Collapsing
=#

"""Appends all blocks in the layout `x` to `out`, including the blocks of nested layouts."""
function _blocks!(out, gl::GridLayout)
    for c in gl.content
        _blocks!(out, c.content)
    end
    return out
end
_blocks!(out, b::Makie.Block) = push!(out, b)
function _blocks!(out, sg::SliderGrid)
    push!(out, sg)
    return _blocks!(out, sg.layout)
end

# Detached parts are laid out here, where they can not take mouse events
const _OFFSCREEN = Point2f(-1.0f5, -1.0f5)

"""
    _set_shown!(part::_AppPart, shown::Bool)

Shows or collapses the `part`. A collapsed part is removed from its parent layout, whose column
or row shrinks to zero, such that the 3D view takes the space. Since Makie keeps drawing, and
buttons keep reacting to clicks within their last bounding boxes, the blocks of the part are also
hidden and laid out off-screen. Its plots are kept, i.e. their state survives collapsing.
"""
function _set_shown!(part::_AppPart, shown::Bool)
    part.shown == shown && return nothing
    part.shown = shown
    blocks = _blocks!(Any[], part.grid)
    if shown
        part.parent[part.pos...] = part.box
        part.parent[part.pos...] = part.grid
        part.resize(part.size)
        Makie.unhide!(part.box)
        foreach(Makie.unhide!, blocks)
    else
        Makie.hide!(part.box)
        foreach(Makie.hide!, blocks)
        for x in (part.box, part.grid)
            _GLB.remove_from_gridlayout!(_GLB.gridcontent(x))
            w = GeometryBasics.widths(x.layoutobservables.computedbbox[])
            x.layoutobservables.suggestedbbox[] = Rect2f(_OFFSCREEN, w)
        end
        part.resize(Fixed(0))
    end
    return nothing
end

"""Shows the dock if its toggle is on and it has panels, otherwise collapses it."""
_update_dock!(layout::AppLayout) =
    _set_shown!(layout.dock, layout.collapse.dock.active[] && !isempty(layout.dock_panels))

#=
Construction
=#

"""Returns a collapsible part with a background box at `pos` of the `parent`, see `_AppPart`."""
function _app_part(parent::GridLayout, pos, resize, size, color; padding = 10)
    box = Box(parent[pos...]; color, cornerradius = 0)
    grid = GridLayout(parent[pos...]; alignmode = Outside(padding), default_rowgap = 8)
    resize(size)
    return _AppPart(parent, Tuple(pos), resize, size, box, grid, true)
end

"""
    _build_toolbar(layout::AppLayout, spec) -> NamedTuple

Creates the entries of the toolbar in its groups (trace, camera, display, tools, panels, help)
and returns them. The entries are placeholder text buttons and toggles, replaced by icons in S3.
"""
function _build_toolbar(layout::AppLayout, spec)
    t = layout.theme
    entry(group) = _add_toolbar_entry!(layout, group)
    toggle(group, label, active) = _TextToggle(entry(group); label, active, theme = t)
    trace_button = Button(entry(:trace); label = "Trace (t)")
    auto_trace_toggle = toggle(:trace, "Auto", spec.auto_trace)
    home_button = Button(entry(:camera); label = "Home")
    fit_button = Button(entry(:camera); label = "Fit (g)")
    views_menu = Menu(entry(:camera); options = _views_options(spec.view_specs), default = nothing,
        prompt = "Views", width = 90)
    save_view_button = Button(entry(:camera); label = "Save view")
    orthographic_toggle = toggle(:display, "Ortho", spec.orthographic)
    clip_toggle = toggle(:display, "Clip (c)", true)
    clip_beams_toggle = toggle(:display, "Beams", spec.clip_beams)
    sources_toggle = toggle(:display, "Sources (1)", spec.show_sources)
    measure_toggle = toggle(:tools, "Measure", false)
    export_button = Button(entry(:tools); label = "Export")
    collapse = (; left = toggle(:panels, "Tree", true), right = toggle(:panels, "Inspector", true),
        dock = toggle(:panels, "Dock", !isempty(spec.specs)))
    help_button = Button(entry(:help); label = "? (h)")
    return (; trace_button, auto_trace_toggle, home_button, fit_button, views_menu,
        save_view_button, orthographic_toggle, clip_toggle, clip_beams_toggle, sources_toggle,
        measure_toggle, export_button, collapse, help_button)
end

"""Creates the "PROPERTIES" section of the right sidebar: selection, pose boxes, step and mode."""
function _build_inspector!(layout::AppLayout)
    t = layout.theme
    g = _add_sidebar_section!(layout, :right, "Properties")
    layout.selection_label = Label(g[1, 1:3], "No selection"; halign = :left, font = :bold,
        fontsize = 14, tellwidth = false)
    layout.type_label = Label(g[2, 1:3], " "; halign = :left, color = t.muted, fontsize = 12,
        tellwidth = false)
    pose_boxes = Textbox[]
    for (k, field) in enumerate(_POSE_FIELDS)
        r, c = 2 * ((k - 1) ÷ 3) + 3, (k - 1) % 3 + 1
        Label(g[r, c], field; halign = :left, fontsize = 11,
            color = k <= 3 ? t.text : t.gizmo[k - 3], tellwidth = false)
        push!(pose_boxes, Textbox(g[r + 1, c]; placeholder = k <= 3 ? " " : "0", width = 72,
            halign = :left))
    end
    Label(g[7, 1], "step"; halign = :left, fontsize = 11)
    step_box = Textbox(g[8, 1:3]; placeholder = "step, e.g. 250 nm", width = 150,
        halign = :left)
    layout.mode_label = Label(g[9, 1:3], "mode: translate"; halign = :left, color = t.muted,
        fontsize = 12, tellwidth = false)
    buttons = GridLayout(g[10, 1:3]; halign = :left)
    hide_button = Button(buttons[1, 1]; label = "Hide")
    show_all_button = Button(buttons[1, 2]; label = "Show all")
    rowgap!(g, 4)
    colgap!(g, 6)
    return (; pose_boxes, step_box, hide_button, show_all_button)
end

"""Creates the placeholder of the object tree, replaced by the tree in S3."""
function _build_tree_placeholder!(layout::AppLayout)
    t = layout.theme
    g = _add_sidebar_section!(layout, :left, "Objects"; grow = true)
    Box(g[1, 1]; color = t.field, strokecolor = t.border, cornerradius = 4)
    Label(g[1, 1], "Objects"; color = t.muted, tellwidth = false, tellheight = false)
    return nothing
end

function _build_layout(layout::AppLayout, fig, spec)
    t = layout.theme
    root = fig.layout
    # Gaps of rows and columns added later, the parts are separated by their borders
    root.default_rowgap = root.default_colgap = Fixed(0)
    main = GridLayout(root[2, 1]; default_colgap = 0)
    # The 3D view first: its scene clears its area, the blocks created later draw on top of it,
    # e.g. the drop-down of the views menu
    ax = LScene(main[1, 2]; show_axis = false,
        scenekw = (; clear = true, backgroundcolor = t.view))
    colsize!(main, 2, Auto(false))
    studio_lighting!(ax; preset = spec.lighting)
    cube = spec.view_cube ? view_cube!(ax) : nothing
    # Sidebars and dock
    layout.left = _app_part(main, (1, 1), s -> colsize!(main, 1, s), Fixed(240), t.sidebar)
    layout.right = _app_part(main, (1, 3), s -> colsize!(main, 3, s), Fixed(260), t.sidebar)
    layout.dock = _app_part(root, (3, 1), s -> rowsize!(root, 3, s),
        Relative(0.3), t.sidebar; padding = 8)
    layout.sections = Dict(:left => Pair{String, GridLayout}[], :right => Pair{String, GridLayout}[])
    layout.fillers = Dict(s => Label(getfield(layout, s).grid[1, 1], ""; tellwidth = false,
        tellheight = false) for s in (:left, :right))
    layout.growing = Dict(:left => false, :right => false)
    layout.dock_panels = Pair{String, GridLayout}[]
    rowsize!(root, 2, Auto(false))
    # Toolbar
    Box(root[1, 1]; color = t.background, cornerradius = 0)
    bar = GridLayout(root[1, 1]; alignmode = Outside(8, 8, 6, 6))
    layout.toolbar = GridLayout(bar[1, 1]; halign = :left, default_colgap = 8)
    Label(bar[1, 2], ""; tellwidth = false)
    layout.groups = Pair{Symbol, GridLayout}[]
    tb = _build_toolbar(layout, spec)
    layout.collapse = tb.collapse
    layout.clip_toggle = tb.clip_toggle
    layout.fit_button = tb.fit_button
    layout.help_button = tb.help_button
    # Sidebars
    _build_tree_placeholder!(layout)
    sliders = if isempty(spec.slider_specs)
        nothing
    else
        g = _add_sidebar_section!(layout, :left, "Parameters")
        SliderGrid(g[1, 1], first.(spec.slider_specs)...; tellwidth = false)
    end
    inspector = _build_inspector!(layout)
    # Analysis dock
    panels = Any[]
    for (i, (pd, mode, kw)) in enumerate(spec.specs)
        name = get(spec.labels, pd, "Detector $i")
        p = DetectorPanel(_add_dock_panel!(layout, name)[1, 1], pd, name, mode, kw)
        # Colors of the panel that do not follow the theme
        p.ax.subtitlecolor[] = t.muted
        p.scatter_plot.color[] = t.text
        push!(panels, p)
    end
    _update_dock!(layout)
    # Status bar
    Box(root[4, 1]; color = t.background, cornerradius = 0)
    sb = GridLayout(root[4, 1]; alignmode = Outside(10, 10, 4, 4))
    status = Label(sb[1, 1], "Click on a component to select it, press h to show the controls";
        halign = :left, tellwidth = false)
    layout.info = Label(sb[1, 2], ""; halign = :right, color = t.muted)
    return (; ax, cube, panels, sliders, status, tb.trace_button, tb.auto_trace_toggle,
        tb.clip_beams_toggle, tb.orthographic_toggle, tb.sources_toggle, inspector.step_box,
        tb.export_button, inspector.hide_button, inspector.show_all_button, inspector.pose_boxes,
        tb.measure_toggle, tb.home_button, tb.save_view_button, tb.views_menu)
end

_build_menus(::AppLayout, w, _, _) = (; menu = nothing, views_menu = w.views_menu)

#=
Connections and status
=#

"""Returns the number of rays (or beams) of the source `beam`, see `_status_info`."""
_ray_count(beam::BMO.AbstractBeamGroup) = length(BMO.beams(beam))
_ray_count(_) = 1

"""Formats the duration `s` [s] of a solve in ms."""
_ms_string(s) = s < 1e-3 ? "<1 ms" : "$(round(Int, 1e3 * s)) ms"

"""Returns the text of the info label of the status bar: last solve, number of rays, projection."""
function _status_info(gui::LiveView)
    n = sum(p -> _ray_count(p.second), gui.pairs)
    traced = gui.preview ? "preview in $(_ms_string(gui.preview_time))" :
             "traced in $(_ms_string(gui.solve_time))"
    projection = gui.orthographic_toggle.active[] ? "orthographic" : "perspective"
    return "$traced · $n $(n == 1 ? "ray" : "rays") · $projection"
end

function _set_text!(label::Label, s::String)
    label.text[] == s || (label.text[] = s)
    return nothing
end

_on_solved!(gui::AppView) = _set_text!(gui.layout.info, _status_info(gui))

function _on_selected!(gui::AppView)
    obj = gui.controls.selected[]
    _set_text!(gui.layout.selection_label, isnothing(obj) ? "No selection" : _label(gui, obj))
    _set_text!(gui.layout.type_label, isnothing(obj) ? " " : string(nameof(typeof(obj))))
    return nothing
end

function _on_clipping!(gui::AppView)
    active = gui.layout.clip_toggle.active
    active[] == gui.clipping || (active[] = gui.clipping)
    return nothing
end

_mode_string(mode::Symbol) = mode == :rotate ? "mode: rotate" : "mode: translate"

"""
Connects the entries of the app layout that are not fields of `LiveView`: the clip toggle, fit,
help, the collapse toggles, and the mode and info labels. All updates are driven by events.
"""
function _connect_layout!(gui::AppView)
    layout = gui.layout
    ctrl = gui.controls
    listeners = ctrl.listeners
    push!(listeners, on(v -> v == gui.clipping || _set_clipping!(gui, v), layout.clip_toggle.active))
    push!(listeners, on(_ -> _zoom_to_selection!(gui), layout.fit_button.clicks))
    push!(listeners, on(layout.help_button.clicks) do _
        ctrl.help_shown = !ctrl.help_shown
        _update_help!(ctrl)
        return nothing
    end)
    push!(listeners, on(v -> _set_shown!(layout.left, v), layout.collapse.left.active))
    push!(listeners, on(v -> _set_shown!(layout.right, v), layout.collapse.right.active))
    push!(listeners, on(_ -> _update_dock!(layout), layout.collapse.dock.active))
    push!(listeners, on(m -> _set_text!(layout.mode_label, _mode_string(m)), ctrl.mode; update = true))
    push!(listeners, on(_ -> _on_solved!(gui), gui.orthographic_toggle.active))
    _on_clipping!(gui)
    return nothing
end
