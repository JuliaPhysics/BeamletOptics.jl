#=
Compact layout of `live_view(...; layout = :compact)`, see `CompactLayout`
=#

"""
    CompactLayout

Layout of `live_view(...; layout = :compact)`: the 3D view with the detector panels on its right,
the sliders, the status row and the tool row below. The positions in `gui.fig` are fixed, e.g. the
panels are placed in `gui.fig[1, 2]`, where users may add their own axes (see
[`add_panel!`](@ref)). The widgets are Makie's buttons, toggles and menus with text labels, in the
colors of the `theme` tokens like the whole window.

The status row holds the tools of the groups `:trace` and `:display` (see `_BUILTIN_TOOLS`), the
status line and the info label; the tool row the component menu, the views menu and the tools of
all other groups, including those of [`add_tool!`](@ref), see `_tool_grid`.
"""
mutable struct CompactLayout <: AbstractLiveLayout
    # the color tokens of the `theme` kwarg, see `_APP_THEMES`
    const theme::NamedTuple
    # the grid of the detector panels in `fig[1, 2]` (`nothing` without panels), the status row
    # and the tool row, the places of the customization API, see `LiveCustom.jl`
    panels::Union{Nothing, GridLayout}
    status_row::GridLayout
    tool_row::GridLayout
    # the grids of the tools in the status row and in the tool row, see `_tool_grid`
    status_tools::GridLayout
    tools::GridLayout
    CompactLayout(theme::NamedTuple) = new(theme, nothing)
end

"""`LiveView` with the compact layout, i.e. `live_view(...; layout = :compact)`."""
const CompactView = LiveView{CompactLayout}

# The colors of the theme, with Makie's padding around the figure, unlike the app layout
_figure(layout::CompactLayout, size) =
    Figure(; size, backgroundcolor = layout.theme.background, _makie_theme(layout.theme)...)

function _build_layout(layout::CompactLayout, fig, spec)
    (; specs, slider_specs, labels, lighting, view_cube) = spec
    t = layout.theme
    ax = LScene(fig[1, 1]; show_axis = false, scenekw = (; clear = true, backgroundcolor = t.view))
    studio_lighting!(ax; preset = lighting)
    # Its click listener runs before the controls, hence clicks on the cube never select objects
    cube = view_cube ? view_cube!(ax) : nothing
    ncols = isempty(specs) ? 1 : 2

    # Detector panels in a near-square grid next to the 3D view
    panels = Any[]
    if !isempty(specs)
        grid = GridLayout(fig[1, 2])
        nc = ceil(Int, sqrt(length(specs)))
        for (i, (pd, mode, kw)) in enumerate(specs)
            parent = grid[(i - 1) ÷ nc + 1, (i - 1) % nc + 1]
            p = DetectorPanel(parent, pd, get(labels, pd, "Detector $i"), mode, kw)
            _theme_panel!(p, t)
            push!(panels, p)
        end
        colsize!(fig.layout, 1, Relative(0.6))
        layout.panels = grid
    end
    sliders = if isempty(slider_specs)
        nothing
    else
        SliderGrid(fig[2, 1:ncols], first.(slider_specs)...)
    end
    # Status row: tools, the status line and the info label at the right
    status_row = GridLayout(fig[isnothing(sliders) ? 2 : 3, 1:ncols])
    layout.status_tools = GridLayout(status_row[1, 1]; default_colgap = 12)
    status = Label(status_row[1, 2],
        "Click on a component to select it, press h to show the controls"; tellwidth = false)
    info = Label(status_row[1, 3], ""; halign = :right, color = t.muted)
    # Tool row: the component menu (added by `_build_menus`, once the movable objects are known)
    # and tools, left-aligned by the filler, which keeps the row compact enough for narrow windows
    tool_row = GridLayout(fig[isnothing(sliders) ? 3 : 4, 1:ncols])
    layout.tools = GridLayout(tool_row[1, 1]; default_colgap = 6)
    Label(tool_row[1, 2], ""; tellwidth = false)
    layout.status_row = status_row
    layout.tool_row = tool_row
    tools = _build_tools(layout, spec)
    return (; ax, cube, panels, sliders, status, info, tools...)
end

function _build_menus(layout::CompactLayout, w, options, views_options)
    tools = layout.tools
    # The component menu first, the views menu after "home"
    _GLB.insertcols!(tools, 1, 1)
    menu = Menu(tools[1, 1]; options, default = nothing, prompt = "select component", width = 150)
    views_menu = Menu(_column_after!(w.home_button); options = views_options, default = nothing,
        prompt = "views", width = 90)
    return (; menu, views_menu)
end


#=
Tools: buttons with their label, toggles with a label next to them
=#

"""
    _tool_grid(layout::CompactLayout, ::Val{group}) -> GridLayout

The grid of the tools of the `group` (see `_BUILTIN_TOOLS`): the status row for tracing and the
display, otherwise the tool row, e.g. for the camera and the tools of [`add_tool!`](@ref).
"""
_tool_grid(layout::CompactLayout, ::Val) = layout.tools
_tool_grid(layout::CompactLayout, ::Union{Val{:trace}, Val{:display}}) = layout.status_tools

"""Returns the position of a new tool at the end of the `grid`, one tool per column."""
_next_tool!(grid::GridLayout) = grid[1, length(grid.content) + 1]

_tool_widget(layout::CompactLayout, group::Symbol, ::Val{false}, label, _, _, ::Bool) =
    Button(_next_tool!(_tool_grid(layout, Val(group))); label)

function _tool_widget(layout::CompactLayout, group::Symbol, ::Val{true}, label, _, _,
        active::Bool)
    g = GridLayout(_next_tool!(_tool_grid(layout, Val(group))); default_colgap = 4)
    toggle = Toggle(g[1, 1]; active)
    Label(g[1, 2], label)
    return toggle
end

#=
Compact layout: panels below the detector panels in `fig[1, 2]`, controls above the status row,
tools in the tool row
=#

"""
Returns the grid of the panels in `fig[1, 2]` of the compact layout, created without detector
panels: then the rows below the 3D view (sliders, status and tool row, controls) span the new
column, like with detector panels.
"""
_panel_grid!(gui::CompactView) = _panel_grid!(gui, gui.layout.panels)
_panel_grid!(::CompactView, grid::GridLayout) = grid
function _panel_grid!(gui::CompactView, ::Nothing)
    root = gui.fig.layout
    for c in copy(root.content)
        c.span.rows.start > 1 && (root[c.span.rows, 1:2] = c.content)
    end
    grid = gui.layout.panels = GridLayout(root[1, 2])
    colsize!(root, 1, Relative(0.6))
    return grid
end

function _add_user_panel!(f, gui::CompactView, title::String, ::Bool)
    grid = _panel_grid!(gui)
    nrows, ncols = size(grid)
    # Below all panels, over all columns of their grid
    row = isempty(grid.content) ? 1 : nrows + 1
    box = GridLayout(grid[row, 1:ncols])
    Label(box[1, 1], title; font = :bold, halign = :left, tellwidth = false)
    layout = GridLayout(box[2, 1])
    rowgap!(box, 4)
    return _user_panel(f, title, layout)
end

function _controls_slot!(gui::CompactView, title::String)
    root = gui.fig.layout
    row = _GLB.gridcontent(gui.layout.status_row).span.rows.start
    _GLB.insertrows!(root, row, 1)
    box = GridLayout(root[row, 1:_GLB.ncols(root)])
    Label(box[1, 1], title; font = :bold)
    layout = GridLayout(box[1, 2])
    # Keeps the controls left-aligned, like the tool row
    Label(box[1, 3], ""; tellwidth = false)
    colgap!(box, 10)
    return layout
end
