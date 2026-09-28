#=
App layout of `live_view(...; layout = :app)`: toolbar, sidebars, analysis dock and status bar
around the 3D view, see `AppLayout`
=#

using Makie: Box, Fixed, Auto, Outside, Rect2f, Point2f, rowsize!

const _GLB = Makie.GridLayoutBase

"""
Color tokens of the app layout by `theme`:

- chrome: `background` (toolbar, status bar), `sidebar`, `view` (background of the 3D view and the
  detector panels), `border`, `text`, `muted` (secondary text), `accent` and `accent_soft` (active
  toggles), `field` (buttons and textboxes), `hover`, `tooltip` and `tooltip_text`
- `gizmo`: the red, green and blue axes of the controls, e.g. the labels of the rotations in the
  inspector, and the x (red) and z (blue) lines of the profiles and the centroid history of the
  detector panels
- 3D view: `help` (text of the controls), `rays` (rays and beams, see `_beam_style`),
  `clip_plane` (marker of a clip plane) and `clip_plane_icon` (its icon in the object tree),
  `marker_stroke` (outline of the handles of sources, clip planes and measured points),
  `gizmo_outline` (outline of the labels of the gizmo, `nothing` for none),
  `card_border` (outline of the floating cards, `nothing` for none, see `_style_card!`) and
  `materials`, the colors of material classes that replace the colors of the render look, e.g.
  of detectors, which are too dark on a dark background, see `_theme_render!`

The light theme keeps the colors of the compact layout in the 3D view.
"""
const _APP_THEMES = Dict{Symbol, NamedTuple}(
    :light => (;
        background = Makie.to_color("#eceef1"), sidebar = Makie.to_color("#f7f8fa"),
        view = Makie.to_color("#ffffff"), border = Makie.to_color("#d3d8de"),
        text = Makie.to_color("#1f2328"), muted = Makie.to_color("#69727d"),
        accent = Makie.to_color("#2f6fdb"), accent_soft = Makie.to_color("#d7e4fa"),
        field = Makie.to_color("#ffffff"), hover = Makie.to_color("#e2e6eb"),
        tooltip = _TOOLTIP_COLOR, tooltip_text = _TOOLTIP_TEXT_COLOR,
        gizmo = Makie.to_color.((:red, :green, :blue)),
        help = Makie.to_color(:gray40), rays = Makie.to_color(:blue),
        clip_plane = Makie.to_color(:purple), clip_plane_icon = Makie.to_color("#9b40c9"),
        marker_stroke = Makie.to_color(:black), gizmo_outline = nothing, card_border = nothing,
        materials = Dict{Symbol, RGBf}()),
    :dark => (;
        background = Makie.to_color("#1e2023"), sidebar = Makie.to_color("#25282c"),
        view = Makie.to_color("#34373b"), border = Makie.to_color("#3d4147"),
        text = Makie.to_color("#e3e5e8"), muted = Makie.to_color("#9aa1a9"),
        accent = Makie.to_color("#6ea4ff"), accent_soft = Makie.to_color("#2b4166"),
        field = Makie.to_color("#1a1c1f"), hover = Makie.to_color("#33373c"),
        tooltip = Makie.to_color("#e3e5e8"), tooltip_text = Makie.to_color("#1f2328"),
        gizmo = Makie.to_color.(("#ff6b6b", "#5ccf5c", "#6ea4ff")),
        help = Makie.to_color("#a4abb3"), rays = Makie.to_color("#5b9bff"),
        clip_plane = Makie.to_color("#b48ef0"), clip_plane_icon = Makie.to_color("#b48ef0"),
        marker_stroke = Makie.to_color("#d5d9de"), gizmo_outline = Makie.to_color("#e3e5e8"),
        card_border = Makie.to_color("#5b6168"),
        materials = Dict{Symbol, RGBf}(:detector => Makie.to_color("#5d6673"),
            :polarizer => Makie.to_color("#56707d"))))

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
        Toggle = (; framecolor_inactive = t.muted, framecolor_active = t.accent,
            buttoncolor = :white),
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
    # expanded for systems, collapsed for groups), the names of unlabelled objects and the counters
    # per type of their running indices, see `_name_objects!`
    tree::_ObjectTree
    expanded::IdDict{Any, Bool}
    names::IdDict{Any, String}
    counters::Dict{String, Int}
    # the inspector with the docked card of the selection, an `_Inspector`, see LiveInspector.jl
    inspector::Any
    # the info label of the status bar
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

"""Keyword arguments of the icon buttons and toggles with the color tokens `t` of the theme."""
_icon_theme(t) = (; icon_color = t.text, hover_color = t.hover, active_color = t.accent_soft,
    active_icon_color = t.accent, tooltip_color = t.tooltip, tooltip_text_color = t.tooltip_text)

"""
    _build_toolbar(layout::AppLayout, spec) -> NamedTuple

Creates the entries of the toolbar in its groups (trace, camera, display, tools, panels, help)
and returns them: flat icon buttons (`clicks`) and toggles (`active`, initialized from `spec`)
with tooltips that name their keys, and the views menu, opened by the icon in front of it.
"""
function _build_toolbar(layout::AppLayout, spec)
    kw = (; _icon_theme(layout.theme)..., size = _TOOLBAR_BUTTON, icon_size = _TOOLBAR_ICON)
    button(group, icon, tooltip) = _IconButton(_add_toolbar_entry!(layout, group); icon, tooltip, kw...)
    toggle(group, icon, tooltip, active) =
        _IconToggle(_add_toolbar_entry!(layout, group); icon, tooltip, active, kw...)
    trace_button = button(:trace, :trace, "Trace (t)")
    auto_trace_toggle = toggle(:trace, :auto_trace, "Auto trace", spec.auto_trace)
    home_button = button(:camera, :home, "Home")
    fit_button = button(:camera, :fit, "Fit to selection (g)")
    views_button = button(:camera, :views, "Views")
    views_menu = Menu(_add_toolbar_entry!(layout, :camera); options = _views_options(spec.view_specs),
        default = nothing, prompt = "Views", width = 84, height = _TOOLBAR_BUTTON - 4,
        textpadding = (6, 4, 3, 3), dropdown_arrow_size = 8, valign = :center,
        # flat like the icon buttons: the closed menu has the color of the toolbar
        selection_cell_color_inactive = layout.theme.background)
    # The icon opens the menu, like a click on the menu itself
    on(_ -> (views_menu.is_open[] = !views_menu.is_open[]), views_button.clicks)
    save_view_button = button(:camera, :save_view, "Save view")
    orthographic_toggle = toggle(:display, :orthographic, "Orthographic", spec.orthographic)
    clip_toggle = toggle(:display, :clip, "Clipping (c)", true)
    clip_beams_toggle = toggle(:display, :clip_beams, "Clip beams", spec.clip_beams)
    sources_toggle = toggle(:display, :sources, "Source markers (1)", spec.show_sources)
    measure_toggle = toggle(:tools, :measure, "Measure", false)
    export_button = button(:tools, :export, "Export changed poses")
    collapse = (; left = toggle(:panels, :panel_left, "Object tree", true),
        right = toggle(:panels, :panel_right, "Properties", true),
        dock = toggle(:panels, :panel_bottom, "Analysis", !isempty(spec.specs)))
    help_button = button(:help, :help, "Help (h)")
    return (; trace_button, auto_trace_toggle, home_button, fit_button, views_button, views_menu,
        save_view_button, orthographic_toggle, clip_toggle, clip_beams_toggle, sources_toggle,
        measure_toggle, export_button, collapse, help_button)
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
Creates the "OBJECTS" section of the left sidebar: the object tree (see `_tree_rows`) and the
"Show all" button in its title row, which replaces the hide buttons of the compact layout together
with the eyes of the tree.
"""
function _build_tree!(layout::AppLayout)
    t = layout.theme
    g = _add_sidebar_section!(layout, :left, "Objects"; grow = true)
    show_all_button = _IconButton(_section_header(layout, :left, "Objects"); icon = :eye,
        tooltip = "Show all", _icon_theme(t)..., icon_color = t.muted, size = 22, icon_size = 16,
        halign = :right, tellwidth = false, tooltip_placement = :right)
    muted = RGBAf(Makie.to_color(t.muted))
    layout.tree = _ObjectTree(g[1, 1]; background = t.sidebar, text_color = t.text,
        muted_color = RGBAf(muted.r, muted.g, muted.b, 0.55), icon_color = t.muted,
        expander_color = t.muted, selection_color = t.accent_soft, accent_color = t.accent,
        guide_color = t.border, scrollbar_color = RGBAf(muted.r, muted.g, muted.b, 0.45),
        marker = _icon, marker_color = kind -> _tree_marker_color(t, kind), marker_size = 16,
        eye_marker = v -> _icon(v ? :eye : :eye_off), eye_size = 15,
        expand_marker = e -> _icon(e ? :collapse : :expand), expand_size = 16)
    layout.expanded = IdDict{Any, Bool}()
    layout.names = IdDict{Any, String}()
    layout.counters = Dict{String, Int}()
    return (; show_all_button)
end

# The markers of sources and clip planes in the tree have the colors of their markers in the 3D
# view, see `_live_render_source!` and `_live_render_clip_plane!`
_tree_marker_color(t::NamedTuple, kind::Symbol) = _tree_marker_color(t, Val(kind))
_tree_marker_color(t::NamedTuple, ::Val) = t.muted
_tree_marker_color(::NamedTuple, ::Val{:source}) = Makie.to_color("#e8890c")
_tree_marker_color(t::NamedTuple, ::Val{:clip_plane}) = t.clip_plane_icon

#=
Colors of the 3D view by the theme, see `_APP_THEMES`
=#

_clip_plane_color(layout::AppLayout) = layout.theme.clip_plane
_marker_stroke(layout::AppLayout) = layout.theme.marker_stroke

# The rays of sources that are drawn as lines, not the envelopes of Gaussian beamlets
_beam_style(layout::AppLayout, ::Union{BMO.AbstractRay, Beam, BMO.AbstractBeamGroup}) =
    (; color = layout.theme.rays)
_beam_style(::AppLayout, ::BMO.AstigmaticBeamGroup) = (;)

_style_card!(layout::AppLayout, c::_ComponentCard) = _style_card!(c, layout.theme.card_border)
_style_card!(::_ComponentCard, ::Nothing) = nothing
function _style_card!(c::_ComponentCard, border)
    c.background.strokecolor[] = border
    c.background.strokewidth[] = 1
    c.background.strokevisible[] = true
    return nothing
end

"""
    _theme_render!(layout::AppLayout, h::SystemRenderHandle)

Replaces the colors of the material classes of the render look by the `materials` of the theme
(see `_APP_THEMES`) in the plots of all objects of `h`, e.g. the dark detectors of the `:modern`
look on a dark background. Only plots in the color of the look are changed, i.e. not the colors
given by the user via `system_kwargs`.
"""
function _theme_render!(layout::AppLayout, h::SystemRenderHandle)
    materials = layout.theme.materials
    isempty(materials) && return nothing
    look = _materials()
    replace = Dict{RGBf, RGBf}(look[class].color => c for (class, c) in materials)
    for oh in h.handles, p in oh.plots
        haskey(p.attributes, :color) && _replace_color!(p, p.color[], replace)
    end
    return nothing
end

# Replaces a single color `c` of the plot `p` by its entry in `replace`, keeping its alpha;
# per-vertex colors, colormaps and the like are kept
function _replace_color!(p, c::Makie.Colors.Colorant, replace)
    new = get(replace, RGBf(c), nothing)
    isnothing(new) || (p.color[] = RGBAf(new, Makie.Colors.alpha(c)))
    return nothing
end
_replace_color!(_, _, _) = nothing

"""
Sets the colors of the overlays of the controls in the 3D view of the `gui` (see
`kinematic_controls!`): the help text, drawn in the block scene of the 3D view, and the outline
of the labels of the gizmo, whose colors are the fixed colors of its axes.
"""
function _theme_controls!(gui::AppView)
    foreach(p -> _theme_control!(gui, p), gui.controls.plots)
    return nothing
end

_theme_control!(::AppView, _) = nothing
function _theme_control!(gui::AppView, p::Makie.Text)
    t = gui.layout.theme
    p.parent === gui.ax.blockscene ? (p.color[] = t.help) : _outline_text!(p, t.gizmo_outline)
    return nothing
end

_outline_text!(_, ::Nothing) = nothing
function _outline_text!(p, color)
    p.strokecolor[] = color
    p.strokewidth[] = 1
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
    tb = _build_toolbar(layout, spec)
    layout.collapse = tb.collapse
    layout.clip_toggle = tb.clip_toggle
    layout.fit_button = tb.fit_button
    layout.views_button = tb.views_button
    layout.help_button = tb.help_button
    # Sidebars
    objects = _build_tree!(layout)
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
    layout.info = Label(sb[1, 2], ""; halign = :right, color = t.muted)
    return (; ax, cube, panels, sliders, status, tb.trace_button, tb.auto_trace_toggle,
        tb.clip_beams_toggle, tb.orthographic_toggle, tb.sources_toggle, inspector.step_box,
        tb.export_button, objects.show_all_button, tb.measure_toggle, tb.home_button,
        tb.save_view_button, tb.views_menu)
end

_build_menus(::AppLayout, w, _, _) = (; menu = nothing, views_menu = w.views_menu)

# The keyboard step is set in the inspector, not on a card, see `_Inspector`
_step_box(::AppLayout, w, _) = w.step_box

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

function _on_solved!(gui::AppView)
    # The inspector shows e.g. the hits of a detector after `_update_inspector!` of `_apply!`
    _set_text!(gui.layout.info, _status_info(gui))
    return nothing
end

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

"""Clip planes of the app layout are numbered, e.g. "Clip plane 2", see `_clip_plane_label`."""
function _clip_plane_label(gui::AppView)
    n = get(gui.layout.counters, "Clip plane", 0) + 1
    gui.layout.counters["Clip plane"] = n
    return "Clip plane $n"
end

"""The label of `obj` in the app layout: its entry of `labels`, else its name in the tree."""
_label(gui::AppView, obj) = get(() -> get(gui.layout.names, obj, string(nameof(typeof(obj)))),
    gui.labels, obj)

function _on_clipping!(gui::AppView)
    active = gui.layout.clip_toggle.active
    active[] == gui.clipping || (active[] = gui.clipping)
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
    _connect_dock!(gui)
    push!(listeners, on(_ -> _on_solved!(gui), gui.orthographic_toggle.active))
    # Object tree
    tree = layout.tree
    push!(listeners, on(key -> _tree_click!(gui, key), tree.clicked))
    push!(listeners, on(key -> _toggle_hidden!(gui, key), tree.eye_clicked))
    push!(listeners, on(key -> _toggle_expanded!(gui, key), tree.expand_clicked))
    _name_objects!(gui)
    _update_tree!(gui)
    _on_clipping!(gui)
    _connect_inspector!(gui)
    _theme_controls!(gui)
    return nothing
end

#=
Object tree
=#

"""
    _tree_kind(obj) -> Symbol

The kind of the row of `obj` in the object tree, which selects its icon, see `_icon`: `:lens`,
`:mirror`, `:beamsplitter`, `:polarizer`, `:detector`, `:group`, `:mesh` (objects without optical
function, e.g. housings), `:source`, `:clip_plane`, `:system` or `:object` for all other objects.
The inspector shows the same icon. Add a method for a new type of component to give it a matching
icon.
"""
_tree_kind(_) = :object
_tree_kind(::Union{BMO.Lens, BMO.DoubletLens, BMO.TripletLens}) = :lens
_tree_kind(::BMO.AbstractReflectiveOptic) = :mirror
_tree_kind(::BMO.AbstractBeamsplitter) = :beamsplitter
_tree_kind(::Union{BMO.AbstractJonesPolarizer, BMO.LinearPolarizer}) = :polarizer
_tree_kind(::BMO.AbstractDetector) = :detector
_tree_kind(::BMO.AbstractObjectGroup) = :group
_tree_kind(::Union{BMO.NonInteractableObject, BMO.IntersectableObject}) = :mesh
_tree_kind(::Union{BMO.AbstractBeam, BMO.AbstractBeamGroup}) = :source
_tree_kind(::LiveClipPlane) = :clip_plane
_tree_kind(::SystemRenderHandle) = :system

"""The rendered objects of the system of `h`, i.e. the leaves of its groups, see `_set_hidden!`."""
_leaves(h::SystemRenderHandle) = _LiveMovable[oh.obj for oh in h.handles]

"""Returns the top-level objects (outermost groups) of the system of `h`, in the order of `h`."""
function _top_levels(h::SystemRenderHandle)
    tops = _LiveMovable[]
    seen = Base.IdSet{Any}()
    for oh in h.handles
        top = _top_level(h, oh.obj)
        top in seen && continue
        push!(seen, top)
        push!(tops, top)
    end
    return tops
end

"""Returns the sources of the `gui`, i.e. the beams of its pairs, without duplicates."""
_sources(gui::LiveView) = unique(objectid, last.(gui.pairs))

"""
    _name_objects!(gui::AppView)

Names the systems ("System i") and all objects and sources without an entry in `labels` by their
type and a running index per type, e.g. "Lens 2", in the order of the tree. Names, once given, are
kept, see `_label`.
"""
function _name_objects!(gui::AppView)
    layout = gui.layout
    function name!(obj)
        (haskey(gui.labels, obj) || haskey(layout.names, obj)) && return nothing
        base = string(nameof(typeof(obj)))
        n = layout.counters[base] = get(layout.counters, base, 0) + 1
        layout.names[obj] = "$base $n"
        return nothing
    end
    for (i, h) in enumerate(gui.system_handles)
        haskey(layout.names, h) || (layout.names[h] = get(gui.labels, h.sys, "System $i"))
        foreach(top -> foreach(name!, _descendants(top)), _top_levels(h))
    end
    haskey(layout.names, gui.extras) || (layout.names[gui.extras] = "Extras")
    foreach(top -> foreach(name!, _descendants(top)), _top_levels(gui.extras))
    foreach(name!, _sources(gui))
    return nothing
end

"""
Returns the state of the eye of the row of `obj`: `nothing` if none of its `_leaves` is rendered,
i.e. it can not be hidden, otherwise whether any of them is shown.
"""
function _tree_visible(gui::LiveView, rendered, obj)
    leaves = filter(leaf -> leaf in rendered, _leaves(obj))
    isempty(leaves) && return nothing
    return !all(leaf -> leaf in gui.hidden, leaves)
end
# Clip planes are switched off, not hidden, see `_toggle_hidden!`
_tree_visible(::LiveView, _, ::LiveClipPlane) = nothing

function _push_tree_rows!(rows, gui::AppView, rendered, obj, depth)
    children = _children(obj)
    expanded = get(gui.layout.expanded, obj, false)
    push!(rows, _TreeRow(obj, _label(gui, obj), depth, _tree_kind(obj), !isempty(children),
        expanded, _tree_visible(gui, rendered, obj)))
    expanded && foreach(c -> _push_tree_rows!(rows, gui, rendered, c, depth + 1), children)
    return rows
end

"""
    _tree_rows(gui::AppView) -> Vector{_TreeRow}

Returns the rows of the object tree of the `gui`: per system a row (key: its `SystemRenderHandle`)
followed by its objects in the hierarchy of the kinematic controls (groups with their objects,
including static objects such as housings, which can be hidden, but not selected), then the
extras (see the `extras` kwarg of `live_view`) in the same way under a row "Extras" (key:
`gui.extras`, only if there are extras), the sources and the clip planes. The keys of all other
rows are the objects. Systems and the extras are expanded, and groups collapsed, by default.
Objects whose plots are all hidden (see `gui.hidden`) are muted.
"""
function _tree_rows(gui::AppView)
    rendered = Base.IdSet{Any}(oh.obj for oh in gui.controls.h.handles)
    rows = _TreeRow[]
    foreach(h -> _push_system_rows!(rows, gui, rendered, h, :system), gui.system_handles)
    isempty(gui.extras.handles) || _push_system_rows!(rows, gui, rendered, gui.extras, :group)
    foreach(src -> _push_tree_rows!(rows, gui, rendered, src, 0), _sources(gui))
    foreach(plane -> _push_tree_rows!(rows, gui, rendered, plane, 0), gui.clip_planes)
    return rows
end

# The row of the system (or the extras) of `h` with the icon of the `kind`, then its objects
function _push_system_rows!(rows, gui::AppView, rendered, h::SystemRenderHandle, kind::Symbol)
    tops = _top_levels(h)
    expanded = get(gui.layout.expanded, h, true)
    push!(rows, _TreeRow(h, _label(gui, h), 0, kind, !isempty(tops), expanded,
        _tree_visible(gui, rendered, h)))
    expanded && foreach(top -> _push_tree_rows!(rows, gui, rendered, top, 1), tops)
    return rows
end

"""Sets the rows of the object tree of the `gui`, only called on events, see `_tree_rows`."""
function _update_tree!(gui::AppView)
    _set_rows!(gui.layout.tree, _tree_rows(gui))
    return nothing
end

"""Expands or collapses the row `key` of the object tree, see `_tree_rows`."""
function _toggle_expanded!(gui::AppView, key)
    expanded = gui.layout.expanded
    expanded[key] = !get(expanded, key, key isa SystemRenderHandle)
    _update_tree!(gui)
    return nothing
end

"""
    _tree_click!(gui, key)

Handles a click on the label of a row of the object tree: selects the object like a click in the
3D view (see `_select!`), a click on a system expands or collapses it. Objects that are not
movable can not be selected.
"""
_tree_click!(gui::AppView, h::SystemRenderHandle) = _toggle_expanded!(gui, h)

function _tree_click!(gui::AppView, obj)
    if !_is_movable(gui.controls, obj)
        gui.status.text[] = "$(_label(gui, obj)) is not movable, it can not be selected"
        return nothing
    end
    _select!(gui, obj)
    return nothing
end

"""
Expands the system and the groups that contain `obj`, such that its row is shown in the object
tree. The rows are only set again if a row was expanded.
"""
_reveal!(::AppView, ::Nothing) = nothing

function _reveal!(gui::AppView, obj)
    ctrl = gui.controls
    expanded = gui.layout.expanded
    chain = _chain(ctrl, obj)
    changed = false
    for group in chain[2:end]
        get(expanded, group, false) && continue
        expanded[group] = changed = true
    end
    top = last(chain)
    for h in (gui.system_handles..., gui.extras)
        any(oh -> _top_level(h, oh.obj) === top, h.handles) || continue
        get(expanded, h, true) && break
        expanded[h] = changed = true
        break
    end
    changed && _update_tree!(gui)
    return nothing
end
