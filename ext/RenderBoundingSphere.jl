# for docs refer to src/Render/RenderBoundingSphere.jl

"""
    _great_circles(center, radius; segments = 64)

Returns the points of the three great circles of the sphere of the `radius` around the `center`
that are parallel to the xy-, yz- and zx-plane, as one polyline for a `lines` plot: each circle is
closed and has the given number of `segments`, the circles are separated by a `NaN` point.
"""
function _great_circles(center, radius; segments::Int = 64)
    c = Point3d(center)
    r = Float64(radius)
    φ = range(0, 2π; length = segments + 1)
    points = Point3d[]
    sizehint!(points, 3 * (segments + 1) + 2)
    # the circle in the plane of the axes (i, j)
    for (k, (i, j)) in enumerate(((1, 2), (2, 3), (3, 1)))
        k > 1 && push!(points, Point3d(NaN))
        for a in φ
            p = zeros(3)
            p[i] = r * cos(a)
            p[j] = r * sin(a)
            push!(points, c + Point3d(p))
        end
    end
    return points
end

# A sphere
render!(::_RenderEnv, ::BMO.NoBoundingSphere; kwargs...) = nothing

function render!(ax::_RenderEnv, sphere::BMO.AbstractBoundingSphere; color = :magenta, linewidth = 1, kwargs...)
    lines!(ax, _great_circles(position(sphere), BMO.radius(sphere)); color, linewidth, kwargs...)
    return nothing
end

"""
    render_bounding_spheres!(ax, x; color = :magenta, main_color = :orange, linewidth = 1, kwargs...)

Draws the bounding spheres that the solver tests for `x`, a shape, an object or an object group,
via `render!(ax, sphere)`, for the `show_bounding_sphere` option of `render!`:

- shape: its sphere
- object of one shape: the sphere of the object, none for an object without a sphere (e.g. a
  `NonInteractableObject`)
- object of several parts (e.g. a doublet lens, a cube beamsplitter) and object group: the sphere of
  each part, since the solver tests each part on its own, and then the main sphere around all parts
  in `main_color`, which the solver tests first. An object without a main sphere, because a part has
  none, gets none. `main_color` also applies to the main spheres of nested parts and groups.

The other `kwargs` are passed on to `render!(ax, sphere)`. Returns `nothing`.
"""
render_bounding_spheres!(ax::_RenderEnv, shape::BMO.AbstractShape; main_color = :orange, kwargs...) =
    render!(ax, BMO.bounding_sphere_of(shape); kwargs...)

render_bounding_spheres!(ax::_RenderEnv, obj::BMO.AbstractObject; kwargs...) =
    render_bounding_spheres!(ax, BMO.shape_trait_of(obj), obj; kwargs...)

# `bounding_sphere_of` of the object, e.g. none for a `NonInteractableObject` whatever its shape
render_bounding_spheres!(ax::_RenderEnv, ::BMO.SingleShape, obj; main_color = :orange, kwargs...) =
    render!(ax, BMO.bounding_sphere_of(obj); kwargs...)

# The solver tests the parts one by one, hence one sphere per part, and the main sphere of the object
# first, which is drawn in `main_color`
render_bounding_spheres!(ax::_RenderEnv, ::BMO.MultiShape, obj; kwargs...) =
    _render_part_spheres!(ax, obj, BMO.shape(obj); kwargs...)

# An object group: the spheres of its objects, then the main sphere of the group
render_bounding_spheres!(ax::_RenderEnv, group::BMO.AbstractObjectGroup; kwargs...) =
    _render_part_spheres!(ax, group, BMO.objects(group); kwargs...)

function _render_part_spheres!(ax::_RenderEnv, whole, parts; main_color = :orange, kwargs...)
    for part in parts
        render_bounding_spheres!(ax, part; main_color, kwargs...)
    end
    render!(ax, BMO.bounding_sphere_of(whole); kwargs..., color = main_color)
    return nothing
end
