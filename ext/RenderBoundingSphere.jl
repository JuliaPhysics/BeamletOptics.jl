# for docs refer to Render.jl

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

render_bounding_sphere!(ax::_RenderEnv, shape::BMO.AbstractShape; kwargs...) =
    render_bounding_sphere!(ax, BMO.world_bounding_sphere(shape); kwargs...)

# A shape without a bounding sphere
render_bounding_sphere!(::_RenderEnv, ::Nothing; kwargs...) = nothing

# A sphere `(center, radius)` in world coordinates
function render_bounding_sphere!(ax::_RenderEnv, sphere::Tuple; color = :magenta, linewidth = 1, kwargs...)
    center, radius = sphere
    lines!(ax, _great_circles(center, radius); color, linewidth, kwargs...)
    return nothing
end

render_bounding_sphere!(ax::_RenderEnv, obj::BMO.AbstractObject; kwargs...) =
    render_bounding_sphere!(ax, BMO.shape_trait_of(obj), obj; kwargs...)

render_bounding_sphere!(ax::_RenderEnv, ::BMO.SingleShape, obj; kwargs...) =
    render_bounding_sphere!(ax, BMO.shape(obj); kwargs...)

# The solver tests the parts one by one, hence one sphere per part
function render_bounding_sphere!(ax::_RenderEnv, ::BMO.MultiShape, obj; kwargs...)
    for part in BMO.shape(obj)
        render_bounding_sphere!(ax, part; kwargs...)
    end
    return nothing
end

function render_bounding_sphere!(ax::_RenderEnv, sys::BMO.AbstractSystem; kwargs...)
    # Avoid use of objects(sys)
    for obj in sys.objects
        render_bounding_sphere!(ax, obj; kwargs...)
    end
    return nothing
end
