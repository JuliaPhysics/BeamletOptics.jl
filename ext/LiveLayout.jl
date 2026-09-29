#=
Layouts, see `AbstractLiveLayout`
=#

"""
    _live_layout(layout::Symbol, theme::Symbol) -> AbstractLiveLayout

Returns the layout of the `layout` and `theme` kwargs of `live_view`, the only place where the
names are mapped to the layout types. The tokens of the `theme` color the floating cards and the
progress window of all layouts, and the whole window of the app layout.
"""
function _live_layout(layout::Symbol, theme::Symbol)
    tokens = _app_theme(theme)
    layout == :compact && return CompactLayout(tokens)
    layout == :app && return AppLayout(tokens)
    throw(ArgumentError("layout must be :compact or :app, got :$layout"))
end

_default_size(::AbstractLiveLayout) = (1400, 800)
_figure(::AbstractLiveLayout, size) = Figure(; size)

# Optional parts of the layout interface, see `AbstractLiveLayout`
_connect_layout!(::LiveView) = nothing
_on_solved!(::LiveView) = nothing
_on_selected!(::LiveView) = nothing
_on_clipping!(::LiveView) = nothing
_on_clip_planes_changed!(::LiveView) = nothing
_on_hidden!(::LiveView) = nothing

# Colors of the 3D view, which a layout may adapt to its theme, see `_APP_THEMES`
_clip_plane_color(::AbstractLiveLayout) = :purple
_marker_stroke(::AbstractLiveLayout) = :black
"""Default kwargs of `live_render!` of the source `beam` in the `layout`, e.g. the color of rays."""
_beam_style(::AbstractLiveLayout, _) = (;)
"""Styles the floating card `c` (see `_ComponentCard`) for the `layout`, e.g. with an outline."""
_theme_render!(::AbstractLiveLayout, _) = nothing

function _slot_error(gui::LiveView, what)
    throw(ArgumentError("the $(nameof(typeof(gui.layout))) of this live view has no $what"))
end
_add_toolbar_entry!(gui::LiveView, _) = _slot_error(gui, "toolbar")
_add_sidebar_section!(gui::LiveView, _, _) = _slot_error(gui, "sidebars")
_add_dock_panel!(gui::LiveView, _) = _slot_error(gui, "dock")

"""
    _step_box(layout, w, card) -> Textbox

The box of the keyboard step of the `layout` with the widgets `w` (see `_build_layout`): the step
box of the `card` of the selection by default, see `_ComponentCard`.
"""
_step_box(::AbstractLiveLayout, _, card) = card.step_box
