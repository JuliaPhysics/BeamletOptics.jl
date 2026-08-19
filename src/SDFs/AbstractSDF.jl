
"""
    AbstractSDF <: AbstractShape

Provides a shape function based on signed distance functions. See https://iquilezles.org/articles/distfunctions/ for more information.

# Implementation reqs.

Subtypes of `AbstractSDF` should implement all reqs. of `AbstractShape` as well as the following:

# Functions

- `sdf(::AbstractSDF, point)`: a function that returns the signed distance for a point in 3D space
"""
abstract type AbstractSDF{T} <: AbstractShape{T} end

function orientation!(sdf::AbstractSDF, dir)
    sdf.dir = dir
    transposed_orientation!(sdf, copy(transpose(dir)))
end

transposed_orientation(sdf::AbstractSDF) = sdf.transposed_dir

transposed_orientation!(sdf::AbstractSDF, tdir) = (sdf.transposed_dir = tdir)

"""
    _world_to_sdf(sdf, point)

Transforms the coordinates of `point` into a reference frame where the `sdf` lies at the origin. Useful to represent translation and rotation.
If rotations are applied, the rotation is applied around the local sdf coordinate system.
"""
function _world_to_sdf(sdf::AbstractSDF, point)
    # transforms world coords to sdf coords
    T = transposed_orientation(sdf)
    # rotates around local xyz system
    return T * (point - position(sdf))
end

"""
    bounding_sphere(sdf)

Returns `nothing` or a `center` point and the `radius` of a sphere which encloses
the shape of the SDF. This function is currently only used for rendering SDFs but
might be used in the future to optimize the raymarching algorithm, by tracing against
the bounding sphere of and SDF first, instead of calling the more costly complex SDF.
"""
bounding_sphere(::AbstractSDF) = nothing

function bounding_box(s::AbstractSDF)
    if bounding_sphere(s) === nothing
        # fallback behavior, use the SDF itself to find a bounding box
        xmin = sdf(s, Point3(-1000, 0, 0)) - 1000
        ymin = sdf(s, Point3(0, -1000, 0)) - 1000
        zmin = sdf(s, Point3(0, 0, -1000)) - 1000
        xmax = 1000 - sdf(s, Point3(1000, 0, 0))
        ymax = 1000 - sdf(s, Point3(0, 1000, 0))
        zmax = 1000 - sdf(s, Point3(0, 0, 1000))
    else
        c_local, r = bounding_sphere(s)
        # Transform center from local SDF coordinates to world coordinates
        center = orientation(s) * c_local + position(s)
        xmin = center[1] - r
        xmax = center[1] + r
        ymin = center[2] - r
        ymax = center[2] + r
        zmin = center[3] - r
        zmax = center[3] + r
    end

    return xmin, xmax, ymin, ymax, zmin, zmax
end

"""
    normal3d(s::AbstractSDF, pos)

Computes the normal vector of `s` at `pos`.
"""
normal3d(s::AbstractSDF, pos) = normal_fd(s, pos)

"""
    surface_tag(sdf::AbstractSDF, point)
    surface_tag(sdf::AbstractSDF, point, normal)

Returns a symbolic surface tag (e.g. `:front`, `:back`, `:side`, `:top`, `:bottom`) for a hit point on the SDF.
Defaults to `:unknown` if no tag is implemented for the given SDF.
"""
surface_tag(s::AbstractSDF, point) = surface_tag(s, _world_to_sdf(s, point), transposed_orientation(s) * normal3d(s, point))
surface_tag(s::AbstractSDF, point, normal) = :unknown

function numeric_gradient(s::AbstractSDF{S}, pos::AbstractArray{R}; h = promote_type(S, R)(1e-6)) where {S, R}
    T = promote_type(S, R)
    p = Point3{T}(pos)
    h_T = T(h)

    # 4-point tetrahedron technique (Inigo Quilez) - only 4 SDF evaluations, isotropic
    k1 = Point3{T}( 1, -1, -1)
    k2 = Point3{T}(-1, -1,  1)
    k3 = Point3{T}(-1,  1, -1)
    k4 = Point3{T}( 1,  1,  1)

    grad = k1 * sdf(s, p + h_T * k1) +
           k2 * sdf(s, p + h_T * k2) +
           k3 * sdf(s, p + h_T * k3) +
           k4 * sdf(s, p + h_T * k4)

    return normalize(grad)
end
numeric_gradient(s::AbstractSDF, pos) = numeric_gradient(s, Point3(pos))

function normal_fd(s::AbstractSDF, p)
    try
        normal = normalize(gradient(x -> sdf(s, x), p))
        if all(!isnan, normal) && norm(normal) > 0
            return normal
        end
    catch
    end
    # fallback
    return numeric_gradient(s, p)
end

"""
    _refine_root_1d(shape, p0, dir, ta, tb; max_iter=6, tol=1e-13)

Refines the 1D intersection distance along the ray `p(t) = p0 + t * dir` within `[ta, tb]` using the secant method.
"""
function _refine_root_1d(shape::AbstractSDF, p0::Point3{T}, dir::Point3{T}, ta::T, tb::T;
        max_iter = 6, tol = Config.get_sdf_raymarch_eps()) where T
    ga = sdf(shape, p0 + ta * dir)
    gb = sdf(shape, p0 + tb * dir)

    abs(ga) <= tol && return ta
    abs(gb) <= tol && return tb

    for _ in 1:max_iter
        denom = gb - ga
        abs(denom) < eps(T) && break

        # Secant step
        t_next = tb - gb * (tb - ta) / denom

        # Guard against wild extrapolation outside [min(ta, tb), max(ta, tb)]
        t_min = min(ta, tb)
        t_max = max(ta, tb)
        if !(t_min <= t_next <= t_max)
            t_next = (ta + tb) / 2
        end

        g_next = sdf(shape, p0 + t_next * dir)
        if abs(g_next) <= tol
            return t_next
        end

        if ga * g_next <= 0
            tb = t_next
            gb = g_next
        else
            ta = tb
            ga = gb
            tb = t_next
            gb = g_next
        end
    end
    return tb
end

"""
    _raymarch_outside(shape::AbstractSDF, pos, dir; num_iter=1000, eps=Config.get_sdf_raymarch_eps(), t_start=0.0)

Perform the ray marching algorithm if the starting pos is outside of `shape`.
Uses adaptive sphere tracing combined with 1D secant root refinement.
"""
function _raymarch_outside(shape::AbstractSDF{S},
        pos::AbstractArray{R},
        dir::AbstractArray{R},
        num_iter = 1000,
        eps = Config.get_sdf_raymarch_eps();
        t_start = zero(promote_type(S, R))) where {S, R}
    T = promote_type(S, R)
    p0 = Point3{T}(pos)
    d = Point3{T}(dir)

    t0 = T(t_start)
    current_pos = p0 + t0 * d
    dist = sdf(shape, current_pos)

    escaped = dist > eps
    prev_t = t0

    i = 1
    while i <= num_iter
        # When near the surface, use a small step; cap the blind escape growth to avoid tunneling
        min_step = escaped ? T(eps) : min(T(1e-4), T(eps) + t0 * T(0.001))
        step_size = max(dist, min_step)

        prev_t = t0
        t0 += step_size
        current_pos = p0 + t0 * d
        dist = sdf(shape, current_pos)
        i += 1

        if dist > eps
            escaped = true
        elseif escaped
            # Hit detected or zero-crossing; refine intersection distance with 1D secant method
            t_hit = _refine_root_1d(shape, p0, d, prev_t, t0; tol = eps)
            hit_pos = p0 + t_hit * d
            normal = normal3d(shape, hit_pos)

            # Filter out false positive hits caused by numerical noise when leaving an SDF
            if !(dot(d, normal) > eps)
                return ShapeIntersection(t_hit, normal, shape)
            end
        end
    end
    return nothing
end

"""
    _raymarch_inside(object::AbstractSDF, pos, dir; num_iter=1000, dl=Config.get_sdf_inside_step())

Perform the ray marching algorithm if the starting pos is inside of `object`.
Marches forward along the ray direction and refines the exit boundary point.
"""
function _raymarch_inside(object::AbstractSDF{S},
        pos::AbstractArray{R},
        dir::AbstractArray{R},
        num_iter = 1000,
        dl = Config.get_sdf_inside_step()) where {S, R}
    T = promote_type(S, R)
    p0 = Point3{T}(pos)
    d = Point3{T}(dir)

    t0 = zero(T)
    current_pos = p0
    dist = sdf(object, current_pos)
    eps = Config.get_sdf_raymarch_eps()

    # If starting on the surface, take a small step inside
    if abs(dist) <= eps
        t0 += T(eps)
        current_pos = p0 + t0 * d
        dist = sdf(object, current_pos)
    end

    prev_t = t0
    i = 1
    while i <= num_iter
        # For exact Euclidean SDFs, abs(dist) is the distance to boundary. dl acts as upper bound.
        step_size = dist < 0 ? min(abs(dist), T(dl)) : T(eps)
        step_size = max(step_size, T(eps))

        prev_t = t0
        t0 += step_size
        current_pos = p0 + t0 * d
        dist = sdf(object, current_pos)

        # Reached or crossed the exit boundary
        if dist > eps
            t_hit = _refine_root_1d(object, p0, d, prev_t, t0; tol = eps)
            hit_pos = p0 + t_hit * d
            normal = normal3d(object, hit_pos)
            # Ensure the ray is actually exiting the surface (not crossing an internal seam)
            if dot(d, normal) > -eps
                return ShapeIntersection(t_hit, normal, object)
            end
        end
        i += 1
    end
    return nothing
end

"""
    bounding_sphere_entry(object::AbstractSDF, ray::AbstractRay)

Computes if the ray intersects the bounding sphere of `object`, and if so, returns `(true, t_enter)`.
If `object` has no bounding sphere, returns `(true, 0.0)`.
"""
function bounding_sphere_entry(object::AbstractSDF, ray::AbstractRay)
    bs = bounding_sphere(object)
    if bs === nothing
        return (true, 0.0) # fallback: no bounding sphere
    end
    c_local, r = bs

    # Transform ray to local SDF coordinates
    local_pos = _world_to_sdf(object, position(ray))
    local_dir = transposed_orientation(object) * direction(ray)

    # Ray-sphere intersection
    v = local_pos - c_local
    b = dot(v, local_dir)
    c = dot(v, v) - r^2

    # If origin is outside (c > 0) and ray points away (b > 0), it misses
    if c > 0 && b > 0
        return (false, 0.0)
    end

    disc = b^2 - c
    if disc < 0
        return (false, 0.0)
    end

    if c > 0
        t_enter = -b - sqrt(disc)
        return (true, max(0.0, t_enter))
    else
        return (true, 0.0)
    end
end

function intersects_bounding_sphere(object::AbstractSDF, ray::AbstractRay)
    return bounding_sphere_entry(object, ray)[1]
end

"""
    intersect3d(object::AbstractSDF, ray::AbstractRay)

Intersection algorithm for sdf based shapes.
"""
function intersect3d(object::AbstractSDF, ray::AbstractRay)
    hit_bs, t_enter = bounding_sphere_entry(object, ray)
    if !hit_bs
        return nothing
    end

    pos = position(ray)
    dir = direction(ray)
    d = sdf(object, pos)

    # Test if outside of sdf, else inside
    if d > Config.get_sdf_surface_threshold()
        # If ray origin is outside bounding sphere, jump empty air up to t_enter (with margin)
        t_start = t_enter > 0 ? max(0.0, t_enter - 100 * Config.get_sdf_raymarch_eps()) : 0.0
        return _raymarch_outside(object, pos, dir; t_start = t_start)
    end

    # Test if normal and ray dir oppose or align to determine if ray enters or exits object
    n = normal3d(object, pos)
    if dot(dir, n) <= 0
        return _raymarch_inside(object, pos, dir)
    else
        return _raymarch_outside(object, pos, dir)
    end
end

# generic SDF transformations
"""
    op_revolve_z(p, sdf2d::Function, offset)

Calculates the SDF at point `p` for the given 2D-SDF function with `offset` by revolving
the 2D shape around the z-axis.

"""
function op_revolve_z(p::Point3{T}, sdf2d::Function, offset = zero(T)) where {T <: Real}
    q = Point2(norm(Point2(p[1], p[2])) - offset, p[3])
    return sdf2d(q)
end

"""
    op_revolve_y(p, sdf2d::Function, offset)

Calculates the SDF at point `p` for the given 2D-SDF function with `offset` by revolving
the 2D shape around the y-axis.

"""
function op_revolve_y(p::Point3{T}, sdf2d::Function, offset = zero(T)) where {T <: Real}
    q = Point2(norm(Point2(p[1], p[3])) - offset, p[2])
    return sdf2d(q)
end

"""
    op_extrude_z(p, sdf2d::Function, height)

Calculates the SDF at point `p` for the given 2D-SDF function and extrudes the shape to
`height` along the z-axis.

"""
function op_extrude_z(p::Point3{T}, sdf2d::Function, height::Real) where {T <: Real}
    d = sdf2d(Point2(p[1], p[2]))
    w = Point2(d, abs(p[3]) - height)

    return min(max(w[1], w[2]), zero(T)) + norm(max.(w, zero(T)))
end

"""
    op_extrude_x(p, sdf2d::Function, height)

Calculates the SDF at point `p` for the given 2D-SDF function and extrudes the shape to
`height` along the x-axis.

"""
function op_extrude_x(p::Point3{T}, sdf2d::Function, height::Real) where {T <: Real}
    d = sdf2d(Point2(p[2], p[3]))
    w = Point2(d, abs(p[1]) - height)

    return min(max(w[1], w[2]), zero(T)) + norm(max.(w, zero(T)))
end
