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
_tree_kind(::Union{SystemRenderHandle, BMO.AbstractSystem}) = :system

"""The rendered objects of the system of `h`, i.e. the leaves of its groups, see `_set_hidden!`."""
_leaves(h::SystemRenderHandle) = _LiveMovable[oh.obj for oh in h.handles]

"""
The objects of the system `sys` as they are rendered, i.e. the leaves of its groups (see
`live_render!`), e.g. for the card of an inspected system, see `_inspect!`.
"""
_leaves(sys::BMO.AbstractSystem) = reduce(vcat, (_leaves(obj) for obj in sys.objects); init = _LiveMovable[])

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
    _name_objects!(gui::LiveView)

Names the systems ("System i") and all objects and sources without an entry in `labels` by their
type and a running index per type, e.g. "Lens 2", in the order of the tree. A system has its name
twice, as its `SystemRenderHandle` (the key of its row and menu entry) and as the system itself,
which its card shows (see `_inspect!`). Names, once given, are kept, see `_label`. The component
menu, if any, shows the names.
"""
function _name_objects!(gui::LiveView)
    layout = gui.objects
    function name!(obj)
        (haskey(gui.labels, obj) || haskey(layout.names, obj)) && return nothing
        base = string(nameof(typeof(obj)))
        n = layout.counters[base] = get(layout.counters, base, 0) + 1
        layout.names[obj] = "$base $n"
        return nothing
    end
    for (i, h) in enumerate(gui.system_handles)
        haskey(layout.names, h) || (layout.names[h] = get(gui.labels, h.sys, "System $i"))
        haskey(layout.names, h.sys) || (layout.names[h.sys] = layout.names[h])
        foreach(top -> foreach(name!, _descendants(top)), _top_levels(h))
    end
    haskey(layout.names, gui.extras) || (layout.names[gui.extras] = "Extras")
    haskey(layout.names, gui.extras.sys) || (layout.names[gui.extras.sys] = "Extras")
    foreach(top -> foreach(name!, _descendants(top)), _top_levels(gui.extras))
    foreach(name!, _sources(gui))
    _refresh_menu_options!(gui, gui.widgets.menu)
    return nothing
end

"""
Returns the state of the eye of the row of `obj`: `nothing` if none of its `_leaves` is rendered,
i.e. it can not be hidden, otherwise whether any of them is shown.
"""
function _tree_visible(gui::LiveView, rendered, obj)
    leaves = filter(leaf -> leaf in rendered, _leaves(obj))
    isempty(leaves) && return nothing
    return !all(leaf -> leaf in gui.objects.hidden, leaves)
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
Objects whose plots are all hidden (see `gui.objects.hidden`) are muted.
"""
function _tree_rows(gui::AppView)
    rendered = Base.IdSet{Any}(oh.obj for oh in gui.controls.h.handles)
    rows = _TreeRow[]
    foreach(h -> _push_system_rows!(rows, gui, rendered, h, :system), gui.system_handles)
    isempty(gui.extras.handles) || _push_system_rows!(rows, gui, rendered, gui.extras, :group)
    foreach(src -> _push_tree_rows!(rows, gui, rendered, src, 0), _sources(gui))
    foreach(plane -> _push_tree_rows!(rows, gui, rendered, plane, 0), gui.clip.planes)
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
3D view (see `_select!`). A click on a system (or the extras) and on an object that is not movable
inspects it instead, i.e. shows its card in the inspector without selecting it, see `_inspect!`;
the expander of a system expands or collapses it.
"""
_tree_click!(gui::AppView, h::SystemRenderHandle) = _inspect!(gui, h)

function _tree_click!(gui::AppView, obj)
    _is_movable(gui.controls, obj) ? _select!(gui, obj) : _inspect!(gui, obj)
    return nothing
end

# The row of the inspected object is highlighted like the selected one, see `_on_selected!`
_show_inspected!(::AppView, ::Nothing) = nothing
function _show_inspected!(gui::AppView, obj)
    _reveal!(gui, obj)
    _set_selected!(gui.layout.tree, _row_key(gui, obj))
    return nothing
end

"""
Expands the system and the groups that contain `obj`, such that its row is shown in the object
tree. The rows are only set again if a row was expanded.
"""
_reveal!(::AppView, ::Nothing) = nothing
# The rows of the systems are always shown
_reveal!(::AppView, ::BMO.AbstractSystem) = nothing

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
