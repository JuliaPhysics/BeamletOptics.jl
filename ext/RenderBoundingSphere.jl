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

# Spheres. `main_color` is consumed by objects and groups, a sphere value has no use for it.
render_bounding_sphere!(::_RenderEnv, ::BMO.NoBoundingSphere; kwargs...) = nothing

function render_bounding_sphere!(
        ax::_RenderEnv, sphere::BMO.AbstractBoundingSphere; color = :magenta, linewidth = 1, main_color = nothing, kwargs...
    )
    lines!(ax, _great_circles(position(sphere), BMO.radius(sphere)); color, linewidth, kwargs...)
    return nothing
end

# A shape
render_bounding_sphere!(ax::_RenderEnv, shape::BMO.AbstractShape; kwargs...) =
    render_bounding_sphere!(ax, BMO.bounding_sphere_of(shape); kwargs...)

# An object: the sphere of its shape, or of each part and then the main sphere of the object
render_bounding_sphere!(ax::_RenderEnv, obj::BMO.AbstractObject; kwargs...) =
    render_bounding_sphere!(ax, BMO.shape_trait_of(obj), obj; kwargs...)

# `bounding_sphere_of` of the object, e.g. none for a `NonInteractableObject` whatever its shape
render_bounding_sphere!(ax::_RenderEnv, ::BMO.SingleShape, obj; kwargs...) =
    render_bounding_sphere!(ax, BMO.bounding_sphere_of(obj); kwargs...)

# The solver tests the parts one by one, hence one sphere per part, and the main sphere of the object
# first, which is drawn in `main_color`
render_bounding_sphere!(ax::_RenderEnv, ::BMO.MultiShape, obj; kwargs...) =
    _render_part_spheres!(ax, obj, BMO.shape(obj); kwargs...)

# An object group: the spheres of its objects, then the main sphere of the group
render_bounding_sphere!(ax::_RenderEnv, group::BMO.AbstractObjectGroup; kwargs...) =
    _render_part_spheres!(ax, group, BMO.objects(group); kwargs...)

function _render_part_spheres!(ax::_RenderEnv, whole, parts; main_color = :orange, kwargs...)
    for part in parts
        render_bounding_sphere!(ax, part; main_color, kwargs...)
    end
    render_bounding_sphere!(ax, BMO.bounding_sphere_of(whole); kwargs..., color = main_color)
    return nothing
end

function render_bounding_sphere!(ax::_RenderEnv, sys::BMO.AbstractSystem; kwargs...)
    # Avoid use of objects(sys)
    for obj in sys.objects
        render_bounding_sphere!(ax, obj; kwargs...)
    end
    return nothing
end
