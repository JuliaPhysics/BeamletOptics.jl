using Makie: translate!, rotate!, Quaternion, AbstractPlot

"""Values that can be live-rendered and moved by the controls: objects and sources."""
const _LiveMovable = Union{BMO.AbstractObject, BMO.AbstractBeam, BMO.AbstractBeamGroup}

"""
    ObjectRenderHandle <: AbstractRenderHandle

Live rendering handle of an `AbstractObject`, see [`live_render!`](@ref), or of the marker of a
source, see `_live_render_source!`.

The plots are generated once by the `draw` function in the reference pose `P0`, `R0` of the object.
Afterwards, the rigid transformation from the reference pose to the current pose is applied as the
model matrix of the plots, i.e. `x ↦ R * R0' * (x - P0) + P`.
"""
mutable struct ObjectRenderHandle{O <: _LiveMovable} <: AbstractRenderHandle
    ax::_RenderEnv
    obj::O
    plots::Vector{AbstractPlot}
    # renders the object in its current pose
    draw::Function
    # reference pose
    P0::Point3{Float64}
    R0::Matrix{Float64}
    # reference poses of the MultiShape subparts
    subP0::Vector{Point3{Float64}}
    subR0::Vector{Matrix{Float64}}
    # last rendered pose
    P::Point3{Float64}
    R::Matrix{Float64}
end

function Base.show(io::IO, h::ObjectRenderHandle)
    print(io, "ObjectRenderHandle(", nameof(typeof(h.obj)), ", ", length(h.plots), " plots)")
end

_pose(x) = _pose(BMO.kinematic_trait_of(x), x)
_pose(::Any, x) = (Point3{Float64}(position(x)), Matrix{Float64}(orientation(x)))
# Beams only have a direction, which is the local y-axis of the frame
_pose(::BMO.Movable{BMO.Directed}, x) = (Point3{Float64}(position(x)), _direction_frame(BMO.direction(x)))

"""Returns a right-handed orthonormal frame with the local y-axis along `d`."""
function _direction_frame(d)
    y = normalize(Vector{Float64}(d))
    x = normalize(cross(y, abs(y[3]) < 0.9 ? [0.0, 0, 1] : [1.0, 0, 0]))
    return hcat(x, y, cross(x, y))
end

function _subposes(obj)
    (obj isa BMO.AbstractObject && BMO.shape_trait_of(obj) isa BMO.MultiShape) ||
        return Point3{Float64}[], Matrix{Float64}[]
    parts = BMO.shape(obj)
    return [_pose(p)[1] for p in parts], [_pose(p)[2] for p in parts]
end

"""Returns the plots that are added to the `ax` by `f()`."""
function _capture_new_plots(f, ax::_RenderEnv)
    before = Set(objectid(p) for p in ax.scene.plots)
    f()
    return AbstractPlot[p for p in ax.scene.plots if objectid(p) ∉ before]
end

"""
    live_render!(ax, obj::AbstractObject; kwargs...)

Renders the `obj` via [`render!`](@ref) and returns an `ObjectRenderHandle`. Kinematic changes of the
`obj` are applied to the plots via [`update_render!`](@ref) without regenerating the geometry.
"""
function live_render!(ax::_RenderEnv, obj::BMO.AbstractObject; kwargs...)
    return _live_render_movable!(ax, obj, () -> render!(ax, obj; kwargs...))
end

function _live_render_movable!(ax::_RenderEnv, obj, draw)
    plots = _capture_new_plots(draw, ax)
    P0, R0 = _pose(obj)
    subP0, subR0 = _subposes(obj)
    return ObjectRenderHandle(ax, obj, plots, draw, P0, R0, subP0, subR0, P0, R0)
end

"""
    _live_render_source!(ax, src; size, color = :orange)

Renders a marker of the source `src` (a beam or beam group), i.e. an arrow of length `size` along
its direction and a sphere at its position, and returns an `ObjectRenderHandle`. The marker allows
selecting and moving the source with the [`kinematic_controls!`](@ref).
"""
function _live_render_source!(ax::_RenderEnv, src; size::Real, color = :orange)
    draw = function ()
        p, d = Point3f(position(src)), Vec3f(size * normalize(BMO.direction(src)))
        arrows3d!(ax, [p], [d]; color, shaftradius = 0.05, tipradius = 0.15, tiplength = 0.35)
        mesh!(ax, GeometryBasics.Sphere(p, Float32(size / 5)); color)
        # Constant size on the screen, such that the source is visible in the overview as well
        scatter!(ax, [p]; color, markersize = 12, strokecolor = :black, strokewidth = 1)
    end
    return _live_render_movable!(ax, src, draw)
end

"""
    _quat_from_rotmatrix(R)

Converts the rotation matrix `R` into a `Makie.Quaternion`.
"""
function _quat_from_rotmatrix(R::AbstractMatrix{T}) where {T}
    tr = R[1, 1] + R[2, 2] + R[3, 3]
    if tr > 0
        S = sqrt(tr + 1) * 2
        w = T(0.25) * S
        x = (R[3, 2] - R[2, 3]) / S
        y = (R[1, 3] - R[3, 1]) / S
        z = (R[2, 1] - R[1, 2]) / S
    elseif R[1, 1] > R[2, 2] && R[1, 1] > R[3, 3]
        S = sqrt(1 + R[1, 1] - R[2, 2] - R[3, 3]) * 2
        w = (R[3, 2] - R[2, 3]) / S
        x = T(0.25) * S
        y = (R[1, 2] + R[2, 1]) / S
        z = (R[1, 3] + R[3, 1]) / S
    elseif R[2, 2] > R[3, 3]
        S = sqrt(1 + R[2, 2] - R[1, 1] - R[3, 3]) * 2
        w = (R[1, 3] - R[3, 1]) / S
        x = (R[1, 2] + R[2, 1]) / S
        y = T(0.25) * S
        z = (R[2, 3] + R[3, 2]) / S
    else
        S = sqrt(1 + R[3, 3] - R[1, 1] - R[2, 2]) * 2
        w = (R[2, 1] - R[1, 2]) / S
        x = (R[1, 3] + R[3, 1]) / S
        y = (R[2, 3] + R[3, 2]) / S
        z = T(0.25) * S
    end
    return Quaternion(x, y, z, w)
end

"""
    _is_rigid(h::ObjectRenderHandle; atol = 1e-6, rtol = 1e-6)

Checks if all subparts of a `MultiShape` object still match a rigid motion of the whole object.
"""
function _is_rigid(h::ObjectRenderHandle; atol = 1e-6, rtol = 1e-6)
    isempty(h.subP0) && return true
    parts = BMO.shape(h.obj)
    length(parts) != length(h.subP0) && return false
    P, R = _pose(h.obj)
    Rd = R * h.R0'
    for (i, part) in enumerate(parts)
        p, r = _pose(part)
        isapprox(Rd * (h.subP0[i] - h.P0) + P, p; atol, rtol) || return false
        isapprox(Rd * h.subR0[i], r; atol, rtol) || return false
    end
    return true
end

function _rerender!(h::ObjectRenderHandle)
    for plot in h.plots
        delete!(h.ax, plot)
    end
    h.plots = _capture_new_plots(h.draw, h.ax)
    h.P0, h.R0 = _pose(h.obj)
    h.P, h.R = h.P0, h.R0
    h.subP0, h.subR0 = _subposes(h.obj)
    return h
end

"""
    update_render!(h::ObjectRenderHandle)

Applies the current pose of the object to its plots. If the subparts of a `MultiShape` object have
not been moved rigidly, the object is rendered again.
"""
function update_render!(h::ObjectRenderHandle)
    # The pose of a MultiShape object is the pose of its first part, hence check the other parts first
    _is_rigid(h) || return _rerender!(h)
    P, R = _pose(h.obj)
    (P == h.P && R == h.R) && return h
    # Makie model matrix: x ↦ Rd * x + t, since scale and origin are not used
    Rd = R * h.R0'
    t = P - Rd * h.P0
    q = _quat_from_rotmatrix(Rd)
    for plot in h.plots
        translate!(plot, t...)
        rotate!(plot, q)
    end
    h.P, h.R = P, R
    return h
end

function remove_render!(h::ObjectRenderHandle)
    for plot in h.plots
        delete!(h.ax, plot)
    end
    empty!(h.plots)
    return nothing
end

"""Walks up the parent plots of `plot` until `isowner` is true, since picking returns primitive child plots."""
function _walk_to_owner(isowner, plot)
    p = plot
    while p isa AbstractPlot
        isowner(p) && return p
        p = p.parent
    end
    return nothing
end

function pick_object(h::ObjectRenderHandle, plot)
    owner = _walk_to_owner(p -> any(q -> q === p, h.plots), plot)
    return isnothing(owner) ? nothing : h.obj
end

"""
    SystemRenderHandle <: AbstractRenderHandle

Live rendering handle of an `AbstractSystem`, which holds one `ObjectRenderHandle` per object.
Object groups are rendered per object, their hierarchy is stored in `parent`, which maps each
object of a group to the enclosing group. Top-level objects have no entry.
"""
mutable struct SystemRenderHandle{S <: BMO.AbstractSystem} <: AbstractRenderHandle
    ax::_RenderEnv
    sys::S
    handles::Vector{ObjectRenderHandle}
    parent::IdDict{BMO.AbstractObject, BMO.AbstractObject}
end

function Base.show(io::IO, h::SystemRenderHandle)
    n = sum(length(oh.plots) for oh in h.handles; init = 0)
    print(io, "SystemRenderHandle(", length(h.handles), " objects, ", n, " plots)")
end

"""
    live_render!(ax, sys::AbstractSystem; kwargs...)

Live-renders all objects of the `sys`tem, see [`live_render!`](@ref). Object groups are rendered per
object, such that each object of a group can be moved on its own without rendering the group again.
[`pick_object`](@ref) returns the top-level object of the system, i.e. the outermost group.
"""
function live_render!(ax::_RenderEnv, sys::BMO.AbstractSystem; kwargs...)
    handles = ObjectRenderHandle[]
    parent = IdDict{BMO.AbstractObject, BMO.AbstractObject}()
    function render_obj!(obj)
        if obj isa BMO.AbstractObjectGroup
            for child in BMO.shape(obj)
                parent[child] = obj
                render_obj!(child)
            end
            return nothing
        end
        push!(handles, live_render!(ax, obj; kwargs...))
        return nothing
    end
    # Avoid use of objects(sys), which flattens the groups
    foreach(render_obj!, sys.objects)
    return SystemRenderHandle(ax, sys, handles, parent)
end

"""Returns the top-level object of `obj` in the hierarchy of `h`, i.e. the outermost group."""
function _top_level(h::SystemRenderHandle, obj)
    while haskey(h.parent, obj)
        obj = h.parent[obj]
    end
    return obj
end

function update_render!(h::SystemRenderHandle)
    foreach(update_render!, h.handles)
    return h
end

function remove_render!(h::SystemRenderHandle)
    foreach(remove_render!, h.handles)
    return nothing
end

"""
    _pick_leaf(h::SystemRenderHandle, plot)

Returns the rendered (leaf) object of `plot`, i.e. an object of a group, or `nothing`. The owner is
looked up in the handles, since a fallback rerender replaces the plots of a handle.
"""
function _pick_leaf(h::SystemRenderHandle, plot)
    p = plot
    while p isa AbstractPlot
        for oh in h.handles
            any(q -> q === p, oh.plots) && return oh.obj
        end
        p = p.parent
    end
    return nothing
end

function pick_object(h::SystemRenderHandle, plot)
    leaf = _pick_leaf(h, plot)
    return isnothing(leaf) ? nothing : _top_level(h, leaf)
end
