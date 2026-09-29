#=
Extras of the live view, i.e. objects that are rendered but not traced, and the opacity of
mechanical objects
=#

"""Objects without optical function, e.g. housings and mounts from CAD files, see `card_rows`."""
const _Mechanics = Union{BMO.NonInteractableObject, BMO.IntersectableObject}

#=
Extras, see the `extras` kwarg of `live_view`
=#

_extra_spec(obj::BMO.AbstractObject) = (obj, (;))
_extra_spec(p::Pair{<:BMO.AbstractObject, <:NamedTuple}) = (p.first, p.second)
function _extra_spec(x)
    throw(ArgumentError("invalid extra $(repr(x)), use `obj` or `obj => (; kwargs...)` with an `AbstractObject` `obj`"))
end

"""
    _extra_specs(extras, systems) -> Vector{Tuple{AbstractObject, NamedTuple}}

Returns `(obj, render_kwargs)` of each entry of the `extras` kwarg of `live_view`. An extra must
not be an object of the `systems` (it would be rendered twice) nor be given twice.
"""
function _extra_specs(extras::Union{AbstractVector, Tuple}, systems)
    specs = [_extra_spec(e) for e in extras]
    traced = Base.IdSet{Any}(o for sys in systems for o in BMO.objects(sys))
    seen = Base.IdSet{Any}()
    for (obj, _) in specs, leaf in _leaves(obj)
        leaf in traced && throw(ArgumentError("the extra $(nameof(typeof(leaf))) is an object of a system, extras are not traced"))
        leaf in seen && throw(ArgumentError("the extra $(nameof(typeof(leaf))) is given twice"))
        push!(seen, leaf)
    end
    return specs
end
_extra_specs(extras, _) = throw(ArgumentError("extras must be a vector of `obj` or `obj => (; kwargs...)`, got $(repr(extras))"))

"""
    _live_render_extras!(ax, specs) -> SystemRenderHandle

Live-renders the extras `specs` (see `_extra_specs`), each object with its own kwargs, like the
objects of a system (groups per object, see `live_render!`). The handle holds a `System` of the
extras, which is never solved, and is listed in the object tree of the app layout as "Extras".
"""
function _live_render_extras!(ax::_RenderEnv, specs)
    handles = ObjectRenderHandle[]
    parent = IdDict{BMO.AbstractObject, BMO.AbstractObject}()
    for (obj, kw) in specs
        h = live_render!(ax, BMO.System(obj); kw...)
        append!(handles, h.handles)
        merge!(parent, h.parent)
    end
    return SystemRenderHandle(ax, BMO.System(BMO.AbstractObject[first.(specs)...]), handles, parent)
end

"""Returns `true` if `obj` is an extra of the `gui` or belongs to one, see `_live_render_extras!`."""
function _is_extra(gui::LiveView, obj)
    top = _top_level(gui.extras, obj)
    return any(o -> o === top, gui.extras.sys.objects)
end

"""
    _on_moved!(gui, obj)

Called after `obj` was moved (or changed via its card with `solve = true`): solves the systems (see
`_on_change!`), except for a clip plane, whose planes are applied again, and an extra, which is
not part of any system.
"""
_on_moved!(gui::LiveView, plane::LiveClipPlane) = _on_clip_change!(gui, plane)
function _on_moved!(gui::LiveView, obj)
    _is_extra(gui, obj) || return _on_change!(gui, obj)
    gui.status.text[] = _pose_string(gui, obj)
    return nothing
end

#=
Opacity, see `_set_opacity!`
=#

"""
    _Opacity

Opacity of an object of a live view set via its card, see `_set_opacity!`: the current `value`, the
`initial` opacity as rendered and, per plot, its `alpha`, `transparency` and opacity as rendered
(`base`).
"""
mutable struct _Opacity
    value::Float64
    initial::Float64
    base::IdDict{AbstractPlot, Tuple{Float32, Bool, Float32}}
end

# Alpha of a single color, per-vertex colors (and none) count as opaque
_color_alpha(c::Union{Makie.Colorant, Symbol, AbstractString, Tuple}) = Float32(Makie.to_color(c).alpha)
_color_alpha(_) = 1.0f0

_plot_alpha(p::AbstractPlot) = haskey(p, :alpha) ? Float32(p.alpha[]) : 1.0f0
# The opacity of a plot: its `alpha` times the alpha of its color
_plot_opacity(p::AbstractPlot) = _plot_alpha(p) * (haskey(p, :color) ? _color_alpha(p.color[]) : 1.0f0)
_plot_base(p::AbstractPlot) = (_plot_alpha(p), Bool(p.transparency[]), _plot_opacity(p))

"""The opacity of an object with the `plots` as rendered: the largest opacity of its plots."""
_rendered_opacity(plots) = isempty(plots) ? 1.0 : Float64(maximum(_plot_opacity, plots))

"""Returns the opacity of `obj` in the `gui` (0 to 1), see `_set_opacity!`."""
_opacity(gui::LiveView, obj) = _opacity(gui, obj, get(gui.objects.opacity, obj, nothing))
_opacity(gui::LiveView, obj, ::Nothing) = _rendered_opacity(_object_plots(gui.controls.h, obj))
_opacity(::LiveView, _, o::_Opacity) = o.value

"""
    _set_opacity!(gui, obj, o)

Sets the opacity of the rendered object `obj` of the `gui` to `o` (0 to 1): the opacity of its most
opaque plot, see `_plot_opacity`; the other plots, e.g. the feature edges, are scaled with it. The
opacity as rendered and the attributes of the plots are stored in `gui.objects.opacity` on the first call
(see `_Opacity`), such that the opacity survives hiding and showing. Only the plots of `obj` change:
plots that are not opaque become `transparency = true` (order independent transparency of GLMakie),
an opaque object keeps the cheaper opaque rendering, see `_set_transparency!`. At 0 the object is
hidden like via the "hide" action, a larger opacity shows it again.
"""
function _set_opacity!(gui::LiveView, obj, o::Real)
    o = clamp(Float64(o), 0.0, 1.0)
    plots = _object_plots(gui.controls.h, obj)
    rec = get!(() -> _Opacity(o, _rendered_opacity(plots), IdDict{AbstractPlot, Tuple{Float32, Bool, Float32}}()),
        gui.objects.opacity, obj)
    rec.value = o
    s = rec.initial > 0 ? o / rec.initial : o
    for p in plots
        _apply_opacity!(p, get!(() -> _plot_base(p), rec.base, p)..., s)
    end
    hidden = _all_hidden(gui, obj)
    if (o <= 0) != hidden
        _set_hidden!(gui, obj, !hidden)
        _on_hidden!(gui)
    end
    gui.status.text[] = "$(_label(gui, obj)): opacity $(round(Int, 100 * o)) %" * (o <= 0 ? ", hidden" : "")
    return nothing
end

"""
    _apply_opacity!(p, alpha, transparency, opacity, s)

Scales the opacity of the plot `p` with the `alpha`, `transparency` and `opacity` as rendered by
`s`, at most to 1. `s = 1` restores the plot as rendered.
"""
function _apply_opacity!(p::AbstractPlot, alpha, transparency, opacity, s)
    f = opacity > 0 ? min(s, 1 / opacity) : s
    a = Float32(alpha * f)
    haskey(p, :alpha) && (p.alpha[] == a || (p.alpha[] = a))
    _set_transparency!(p, s == 1 ? transparency : (transparency || s < 1) && opacity * f < 1)
    return nothing
end

"""
    _set_transparency!(p, transparency)

Sets the `transparency` of the plot `p`. GLMakie builds the render object of a plot for its
transparency (the shader writes into the buffers of the order independent transparency) and ignores
later changes, hence the render object is rebuilt in each screen of the plot, like Makie moves a
plot. The geometry stays on the CPU side of the plot; for a mesh with 1 M triangles, the rebuild
takes a few ms.
"""
function _set_transparency!(p::AbstractPlot, transparency::Bool)
    p.transparency[] == transparency && return nothing
    p.transparency[] = transparency
    scene = Makie.parent_scene(p)
    for screen in scene.current_screens
        delete!(screen, scene, p)
        insert!(screen, scene, p)
    end
    return nothing
end

"""
Opacity below which a click in the 3D view does not select mechanics (see `_Mechanics`): a click
into the empty space inside or over a transparent housing is meant for what lies behind it, see
`_pickable`. They stay selectable via the object tree and the component menu.
"""
const _MIN_PICK_OPACITY = 0.5

# Mechanics that are (almost) see-through are not selected by a click, see `_MIN_PICK_OPACITY`. The
# opacity is that of the plots, i.e. as rendered or as set via `_set_opacity!`
_pickable(::KinematicController, ::_Mechanics, plots) =
    !_is_hidden(plots) && _rendered_opacity(plots) >= _MIN_PICK_OPACITY

"""Gives an object of the `gui` shown again at opacity 0 its initial opacity, see `_set_hidden!`."""
_restore_opacity!(::LiveView, _, ::Nothing) = nothing
function _restore_opacity!(gui::LiveView, leaf, o::_Opacity)
    o.value > 0 && return nothing
    _set_opacity!(gui, leaf, o.initial > 0 ? o.initial : 1.0)
    return nothing
end

#=
Card rows of the mechanics
=#

_opacity_percent(gui::LiveView, obj) = round(Int, 100 * _opacity(gui, obj))

"""The row of the opacity of an object on its card: a slider 0-100 %, see `_set_opacity!`."""
_opacity_row() = CardRow("opacity",
    CardWidget(Slider; name = :opacity, range = 0:100, width = 150, value = _opacity_percent,
        on = (gui, obj, v) -> _set_opacity!(gui, obj, v / 100)),
    CardWidget(Label; name = :opacity_value, width = 40, halign = :right,
        value = (gui, obj) -> "$(_opacity_percent(gui, obj)) %"))

# Mechanics, e.g. a housing: the pose and the opacity, which does not change the optics
card_rows(obj::_Mechanics) = (pose_card_rows(obj)..., _opacity_row())
