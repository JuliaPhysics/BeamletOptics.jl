#=
App layout of `live_view(...; layout = :app)`: toolbar, sidebars, analysis dock and status bar
around the 3D view, see `AppLayout`
=#

using Makie: Box, Fixed, Auto, Outside, Rect2f, Point2f, rowsize!

const _GLB = Makie.GridLayoutBase

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
    AppLayout

Layout of `live_view(...; layout = :app)`, an application window around the 3D view:

- a toolbar at the top, an ordered list of groups of entries, see `_add_toolbar_entry!`
- the left sidebar, a stack of titled sections ("OBJECTS" with the object tree, see `_tree_rows`,
  "PARAMETERS" with the sliders), and the right sidebar ("PROPERTIES": the inspector, see
  `_Inspector`), see `_add_sidebar_section!`
- the analysis dock below the 3D view, a tab per detector panel (or other panel, see
  `_add_dock_panel!`), of which only the active one is shown and computed, see `_DockTabs`
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
    # the panels `title => content` of the dock and its tabs, a `_DockTabs`
    dock_panels::Vector{Pair{String, GridLayout}}
    tabs::Any
    # toolbar entries without a field of `LiveView`
    collapse::@NamedTuple{left::_IconToggle, right::_IconToggle, dock::_IconToggle}
    clip_toggle::_IconToggle
    fit_button::_IconButton
    views_button::_IconButton
    help_button::_IconButton
    # object tree, see `_tree_rows`: the expanded systems and groups (by key, the default is
    # expanded for systems, collapsed for groups); the names of its rows are in `gui.objects`, see
    # `_name_objects!`
    tree::_ObjectTree
    expanded::IdDict{Any, Bool}
    # the inspector with the docked card of the selection, an `_Inspector`, see LiveInspector.jl
    inspector::Any
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
group is appended to the toolbar, after a separator, but before the group `:help`, which stays
last.
"""
function _add_toolbar_entry!(layout::AppLayout, group::Symbol)
    i = findfirst(g -> g.first == group, layout.groups)
    grid = if isnothing(i)
        k = length(layout.groups)
        if k > 0
            Box(layout.toolbar[1, 2k]; width = 1, height = 22, color = layout.theme.border,
                strokewidth = 0)
        end
        # The group `:help` moves one place to the right, the new group takes its place
        j = k > 0 && layout.groups[k].first == :help ? k : k + 1
        j == k && (layout.toolbar[1, 2k + 1] = layout.groups[k].second)
        g = GridLayout(layout.toolbar[1, 2j - 1]; default_colgap = 4)
        insert!(layout.groups, j, group => g)
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

# `_add_dock_panel!` is defined with the tabs of the dock in `LiveDock.jl`

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

# Size of the icon buttons of the toolbar and of the icons on them [px]
const _TOOLBAR_BUTTON = 30
const _TOOLBAR_ICON = 20

#=
Tools: the built-in tools (see `_BUILTIN_TOOLS`) and those of `add_tool!` as flat icon buttons and
toggles in the groups of the toolbar, "show all" in the title row of the object tree
=#

# The app layout has all built-in tools
_has_tool(::AppLayout, ::Val) = true

function _tool_widget(layout::AppLayout, group::Symbol, toggle::Val, _, icon, tooltip::String,
        active::Bool)
    pos, kw = _tool_slot!(layout, Val(group))
    return _icon_tool(toggle, pos; icon, tooltip, active, kw...)
end

_icon_tool(::Val{false}, pos; active, kw...) = _IconButton(pos; kw...)
_icon_tool(::Val{true}, pos; kw...) = _IconToggle(pos; kw...)

"""
    _tool_slot!(layout::AppLayout, ::Val{group}) -> (pos, kwargs)

The position and the kwargs of the icon button or toggle of a new tool of the `group`: at the end
of the toolbar `group`, see `_add_toolbar_entry!`, or, for the group `:objects` ("show all"), at
the right of the title of the object tree, since the eyes of the tree hide and show objects.
"""
_tool_slot!(layout::AppLayout, ::Val{G}) where {G} =
    (_add_toolbar_entry!(layout, G),
        (; _icon_theme(layout.theme)..., size = _TOOLBAR_BUTTON, icon_size = _TOOLBAR_ICON))

function _tool_slot!(layout::AppLayout, ::Val{:objects})
    t = layout.theme
    return (_section_header(layout, :left, "Objects"),
        (; _icon_theme(t)..., icon_color = t.muted, size = 22, icon_size = 16, halign = :right,
            tellwidth = false, tooltip_placement = :right))
end

"""
Returns the position of the title row of the section `title` of the `side`bar, e.g. for a button
at the right of the title.
"""
function _section_header(layout::AppLayout, side::Symbol, title::AbstractString)
    i = findfirst(s -> s.first == title, layout.sections[side])
    isnothing(i) && throw(ArgumentError("the $side sidebar has no section \"$title\""))
    return getfield(layout, side).grid[2i - 1, 1]
end

"""
Creates the "OBJECTS" section of the left sidebar with the object tree (see `_tree_rows`). The
"Show all" button in its title row is a tool, see `_tool_slot!`.
"""
function _build_tree!(layout::AppLayout)
    t = layout.theme
    g = _add_sidebar_section!(layout, :left, "Objects"; grow = true)
    muted = RGBAf(Makie.to_color(t.muted))
    layout.tree = _ObjectTree(g[1, 1]; background = t.sidebar, text_color = t.text,
        muted_color = RGBAf(muted.r, muted.g, muted.b, 0.55), icon_color = t.muted,
        expander_color = t.muted, selection_color = t.accent_soft, accent_color = t.accent,
        guide_color = t.border, scrollbar_color = RGBAf(muted.r, muted.g, muted.b, 0.45),
        marker = _icon, marker_color = kind -> _tree_marker_color(t, kind), marker_size = 16,
        eye_marker = v -> _icon(v ? :eye : :eye_off), eye_size = 15,
        expand_marker = e -> _icon(e ? :collapse : :expand), expand_size = 16)
    layout.expanded = IdDict{Any, Bool}()
    return nothing
end

# The markers of sources and clip planes in the tree have the colors of their markers in the 3D
# view, see `_live_render_source!` and `_live_render_clip_plane!`
_tree_marker_color(t::NamedTuple, kind::Symbol) = _tree_marker_color(t, Val(kind))
_tree_marker_color(t::NamedTuple, ::Val) = t.muted
_tree_marker_color(::NamedTuple, ::Val{:source}) = Makie.to_color("#e8890c")
_tree_marker_color(t::NamedTuple, ::Val{:clip_plane}) = t.clip_plane_icon

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
    layout.right = _app_part(main, (1, 3), s -> colsize!(main, 3, s), Fixed(300), t.sidebar)
    layout.dock = _app_part(root, (3, 1), s -> rowsize!(root, 3, s),
        Relative(0.36), t.sidebar; padding = 8)
    layout.sections = Dict(:left => Pair{String, GridLayout}[], :right => Pair{String, GridLayout}[])
    layout.fillers = Dict(s => Label(getfield(layout, s).grid[1, 1], ""; tellwidth = false,
        tellheight = false) for s in (:left, :right))
    layout.growing = Dict(:left => false, :right => false)
    rowsize!(root, 2, Auto(false))
    # Toolbar
    Box(root[1, 1]; color = t.background, cornerradius = 0)
    bar = GridLayout(root[1, 1]; alignmode = Outside(8, 8, 6, 6))
    layout.toolbar = GridLayout(bar[1, 1]; halign = :left, default_colgap = 8)
    Label(bar[1, 2], ""; tellwidth = false)
    layout.groups = Pair{Symbol, GridLayout}[]
    # The object tree before the tools, "show all" is placed in its title row
    _build_tree!(layout)
    tb = _build_tools(layout, spec)
    layout.collapse = (; left = tb.collapse_left, right = tb.collapse_right, dock = tb.collapse_dock)
    layout.clip_toggle = tb.clip_toggle
    layout.fit_button = tb.fit_button
    layout.views_button = tb.views_button
    layout.help_button = tb.help_button
    # Sidebars
    sliders = if isempty(spec.slider_specs)
        nothing
    else
        g = _add_sidebar_section!(layout, :left, "Parameters")
        SliderGrid(g[1, 1], first.(spec.slider_specs)...; tellwidth = false)
    end
    inspector = _build_inspector!(layout)
    # Analysis dock
    panels = _build_dock!(layout, spec)
    _update_dock!(layout)
    # Status bar
    Box(root[4, 1]; color = t.background, cornerradius = 0)
    sb = GridLayout(root[4, 1]; alignmode = Outside(10, 10, 4, 4))
    status = Label(sb[1, 1], "Click on a component to select it, press h to show the controls";
        halign = :left, tellwidth = false)
    info = Label(sb[1, 2], ""; halign = :right, color = t.muted)
    return (; ax, cube, panels, sliders, status, tb.trace_button, tb.auto_trace_toggle,
        tb.clip_beams_toggle, tb.orthographic_toggle, tb.sources_toggle, inspector.step_box,
        tb.export_button, tb.show_all_button, tb.measure_toggle, tb.home_button,
        tb.save_view_button, info)
end

"""
The app layout has no component menu (the object tree replaces it); the views menu is placed in
the toolbar after the views icon, which opens it.
"""
function _build_menus(layout::AppLayout, _, _, views_options)
    views_menu = Menu(_column_after!(layout.views_button.box); options = views_options,
        default = nothing, prompt = "Views", width = 84, height = _TOOLBAR_BUTTON - 4,
        textpadding = (6, 4, 3, 3), dropdown_arrow_size = 8, valign = :center,
        # flat like the icon buttons: the closed menu has the color of the toolbar
        selection_cell_color_inactive = layout.theme.background)
    # The icon opens the menu, like a click on the menu itself
    on(_ -> (views_menu.is_open[] = !views_menu.is_open[]), layout.views_button.clicks)
    return (; menu = nothing, views_menu)
end

# The keyboard step is set in the inspector, not on a card, see `_Inspector`
_step_box(::AppLayout, w, _) = w.step_box

#=
Connections and status
=#

# The inspector follows the selection via `_update_inspector!`, called right before
function _on_selected!(gui::AppView)
    obj = gui.controls.selected[]
    _reveal!(gui, obj)
    _set_selected!(gui.layout.tree, obj)
    return nothing
end

_on_clip_planes_changed!(gui::AppView) = _update_tree!(gui)
_on_hidden!(gui::AppView) = _update_tree!(gui)
_show_hint(::AppView) = "click its eye in the object tree to show it again"

function _on_clipping!(gui::AppView)
    active = gui.layout.clip_toggle.active
    active[] == gui.clip.enabled || (active[] = gui.clip.enabled)
    return nothing
end

"""
Connects the entries of the app layout that are not fields of `LiveView`: the clip toggle, fit,
help, the collapse toggles, the info label, the object tree and the inspector. All updates are
driven by events.
"""
function _connect_layout!(gui::AppView)
    layout = gui.layout
    ctrl = gui.controls
    listeners = ctrl.listeners
    push!(listeners, on(v -> v == gui.clip.enabled || _set_clipping!(gui, v), layout.clip_toggle.active))
    push!(listeners, on(_ -> _zoom_to_selection!(gui), layout.fit_button.clicks))
    push!(listeners, on(layout.help_button.clicks) do _
        ctrl.help_shown = !ctrl.help_shown
        _update_help!(ctrl)
        return nothing
    end)
    push!(listeners, on(v -> _set_shown!(layout.left, v), layout.collapse.left.active))
    push!(listeners, on(v -> _set_shown!(layout.right, v), layout.collapse.right.active))
    push!(listeners, on(_ -> _update_dock!(layout), layout.collapse.dock.active))
    _connect_dock!(gui)
    # Object tree
    tree = layout.tree
    push!(listeners, on(key -> _tree_click!(gui, key), tree.clicked))
    push!(listeners, on(key -> _toggle_hidden!(gui, key), tree.eye_clicked))
    push!(listeners, on(key -> _toggle_expanded!(gui, key), tree.expand_clicked))
    _update_tree!(gui)
    _on_clipping!(gui)
    _connect_inspector!(gui)
    return nothing
end


#=
App layout: panels as tabs of the analysis dock, controls as sections of the left sidebar, tools
in the toolbar group `:user` before "Help"
=#

_mark_panel_stale!(gui::AppView, p::_UserPanel) = (push!(gui.layout.tabs.stale, p); nothing)
_results_valid(gui::AppView) = gui.layout.tabs.hits_valid && !_running(gui) && !gui.trace.preview

"""Hides the blocks of the collapsed `part`, e.g. those added to it since it was collapsed."""
function _hide_collapsed!(part::_AppPart)
    part.shown || foreach(Makie.hide!, _blocks!(Any[], part.grid))
    return nothing
end

"""
Adds the panel as a tab of the dock: the content is built into the shown tab, then the previously
active tab is shown again, unless `select` (or there was none). The panel is updated when its tab
is opened, see `_refresh_tab!`.
"""
function _add_user_panel!(f, gui::AppView, title::String, select::Bool)
    layout = gui.layout
    tabs = layout.tabs
    old = tabs.active
    content = _add_dock_panel!(layout, title; icon = :chart)
    i = tabs.active
    p = try
        _user_panel(f, title, content)
    finally
        (select || old == 0) || _select_tab!(layout, old)
        _hide_collapsed!(layout.dock)
    end
    tabs.panels[i] = p
    return p
end

function _refresh_tab!(gui::AppView, p::_UserPanel)
    _results_valid(gui) || return nothing
    _update_user_panel!(gui, p)
    delete!(gui.layout.tabs.stale, p)
    return nothing
end

_controls_slot!(gui::AppView, title::String) = _add_sidebar_section!(gui, :left, title)
_on_controls_added!(gui::AppView) = _hide_collapsed!(gui.layout.left)

