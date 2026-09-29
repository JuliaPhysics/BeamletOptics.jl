#=
Compact layout of `live_view(...; layout = :compact)`, see `CompactLayout`
=#

"""
    CompactLayout

Layout of `live_view(...; layout = :compact)`: the 3D view with the detector panels on its right,
the sliders, the status row and the tool row below. The positions in `gui.fig` are fixed, e.g. the
panels are placed in `gui.fig[1, 2]`, where users may add their own axes (see
[`add_panel!`](@ref)). The widgets have Makie's look, the floating cards and the progress window
the colors of the `theme` tokens.
"""
mutable struct CompactLayout <: AbstractLiveLayout
    # the color tokens of the `theme` kwarg, see `_APP_THEMES`
    const theme::NamedTuple
    # the grid of the detector panels in `fig[1, 2]` (`nothing` without panels), the status row
    # and the tool row, the places of the customization API, see `LiveCustom.jl`
    panels::Union{Nothing, GridLayout}
    status_row::GridLayout
    tool_row::GridLayout
    CompactLayout(theme::NamedTuple) = new(theme, nothing)
end

"""`LiveView` with the compact layout, i.e. `live_view(...; layout = :compact)`."""
const CompactView = LiveView{CompactLayout}

function _build_layout(layout::CompactLayout, fig, spec)
    (; specs, slider_specs, labels, lighting, view_cube, auto_trace, clip_beams, orthographic,
        show_sources) = spec
    ax = LScene(fig[1, 1]; show_axis = false)
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
            push!(panels, DetectorPanel(parent, pd, get(labels, pd, "Detector $i"), mode, kw))
        end
        colsize!(fig.layout, 1, Relative(0.6))
        layout.panels = grid
    end
    sliders = if isempty(slider_specs)
        nothing
    else
        SliderGrid(fig[2, 1:ncols], first.(slider_specs)...)
    end
    # Status row: trace button, auto trace toggle and status line
    status_row = GridLayout(fig[isnothing(sliders) ? 2 : 3, 1:ncols])
    trace_button = Button(status_row[1, 1]; label = "Trace (t)")
    auto_trace_toggle = Toggle(status_row[1, 2]; active = auto_trace)
    Label(status_row[1, 3], "auto trace")
    clip_beams_toggle = Toggle(status_row[1, 4]; active = clip_beams)
    Label(status_row[1, 5], "clip beams")
    orthographic_toggle = Toggle(status_row[1, 6]; active = orthographic)
    Label(status_row[1, 7], "orthographic")
    sources_toggle = Toggle(status_row[1, 8]; active = show_sources)
    Label(status_row[1, 9], "sources (1)")
    export_button = Button(status_row[1, 10]; label = "Export")
    status = Label(status_row[1, 11],
        "Click on a component to select it, press h to show the controls"; tellwidth = false)
    # Tool row: component menu (added by `_build_menus`, once the movable objects are known) and
    # "show all"
    tool_row = GridLayout(fig[isnothing(sliders) ? 3 : 4, 1:ncols])
    show_all_button = Button(tool_row[1, 2]; label = "show all")
    # Measuring and camera tools, the views menu is added with the component menu
    measure_toggle = Toggle(tool_row[1, 3]; active = false)
    Label(tool_row[1, 4], "measure")
    home_button = Button(tool_row[1, 5]; label = "home")
    save_view_button = Button(tool_row[1, 7]; label = "save view")
    # Keeps the tool row left-aligned and compact enough for narrow windows
    Label(tool_row[1, 8], ""; tellwidth = false)
    colgap!(tool_row, 6)
    layout.status_row = status_row
    layout.tool_row = tool_row
    return (; ax, cube, panels, sliders, status, trace_button, auto_trace_toggle, clip_beams_toggle,
        orthographic_toggle, sources_toggle, export_button, show_all_button, measure_toggle,
        home_button, save_view_button, tool_row)
end

function _build_menus(::CompactLayout, w, options, views_options)
    menu = Menu(w.tool_row[1, 1]; options, default = nothing, prompt = "select component",
        width = 150)
    views_menu = Menu(w.tool_row[1, 6]; options = views_options, default = nothing,
        prompt = "views", width = 90)
    return (; menu, views_menu)
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

"""Returns `n` positions at the end of the tool row of the compact `layout`, before its filler."""
function _tool_slots!(layout::CompactLayout, n::Int)
    row = layout.tool_row
    k = _GLB.ncols(row)
    row[1, k + n] = only(_GLB.contents(row[1, k]))
    return [row[1, k + i] for i in 0:(n - 1)]
end

_tool_widget(gui::CompactView, ::Val{false}, label, _, _) =
    Button(only(_tool_slots!(gui.layout, 1)); label)

function _tool_widget(gui::CompactView, ::Val{true}, label, _, _)
    toggle_pos, label_pos = _tool_slots!(gui.layout, 2)
    toggle = Toggle(toggle_pos; active = false)
    Label(label_pos, label)
    return toggle
end

