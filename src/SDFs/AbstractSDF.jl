
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
    bounding_box(s::AbstractSDF)

Returns the limits `(xmin, xmax, ymin, ymax, zmin, zmax)` of an axis-aligned box around `s` in world
coordinates, e.g. for the sampling grid of the rendering. The SDF is probed from ±1000 m along the
axes, which gives the extent of the shape if the SDF is exact outside of it. A type whose SDF is
not exact far away returns the box around its [`bounding_sphere_of`](@ref) instead.
"""
function bounding_box(s::AbstractSDF)
    xmin = sdf(s, Point3(-1000, 0, 0)) - 1000
    ymin = sdf(s, Point3(0, -1000, 0)) - 1000
    zmin = sdf(s, Point3(0, 0, -1000)) - 1000
    xmax = 1000 - sdf(s, Point3(1000, 0, 0))
    ymax = 1000 - sdf(s, Point3(0, 1000, 0))
    zmax = 1000 - sdf(s, Point3(0, 0, 1000))
    return xmin, xmax, ymin, ymax, zmin, zmax
end

"""
    normal3d(s::AbstractSDF, pos)

Computes the normal vector of `s` at `pos`.
"""
normal3d(s::AbstractSDF, pos) = normal_fd(s, pos)

function numeric_gradient(s::AbstractSDF, pos)
    # approximate ∇ of s at pos
    eps = 1e-8
    norm = Point3(sdf(s, pos + Point3(eps, 0, 0)) - sdf(s, pos - Point3(eps, 0, 0)),
        sdf(s, pos + Point3(0, eps, 0)) - sdf(s, pos - Point3(0, eps, 0)),
        sdf(s, pos + Point3(0, 0, eps)) - sdf(s, pos - Point3(0, 0, eps)))
    return normalize(norm)
end

function normal_fd(s::AbstractSDF, p)
    normal = normalize(gradient(x -> sdf(s, x), p))
    all(!isnan, normal) && return normal
    # fallback
    return numeric_gradient(s, p)
end

"""
    SDF_MISS_DISTANCE

Path length in [m] after which a ray that is marched towards an `AbstractSDF` counts as a
miss. A shape that lies farther away than this from the start of the ray is not hit.
"""
const SDF_MISS_DISTANCE = 1e15

"""
    _raymarch_outside(shape::AbstractSDF, pos, dir, t_max=SDF_MISS_DISTANCE, num_iter=1000, eps=1e-10)

Perform the ray marching algorithm if the starting pos is outside of `shape`. The ray is a miss
once it has travelled the path length `t_max` in [m] without a hit.
"""
function _raymarch_outside(shape::AbstractSDF{S},
        pos::AbstractArray{R},
        dir::AbstractArray{R},
        t_max = SDF_MISS_DISTANCE,
        num_iter = 1000,
        eps = Config.get_sdf_raymarch_eps()) where {S, R}
    T = promote_type(S, R)
    dist = sdf(shape, pos)
    t0 = zero(T)

    # `escaped` tracks if the ray has definitively moved away from the surface
    escaped = dist > eps
    i = 1
    while i <= num_iter
        # When trapped in the surface noise floor (dist < eps), we slowly accelerate the minimum step
        # proportionally to the distance traveled (t0 * 0.01). 
        # This logarithmic escape prevents exhausting num_iter on highly inaccurate SDFs, 
        # while bounding the blind step to 1% of the traveled distance to prevent tunneling.
        min_step = escaped ? eps : (eps + t0 * 0.01)
        step_size = max(dist, min_step)
        pos = pos + step_size * dir
        t0 += step_size
        dist = sdf(shape, pos)
        i += 1

        # A ray escaping the scene makes `dist` track the true remaining distance, so `t0`
        # doubles with every step. Up to the overflow to `Inf` that takes about 1000 steps,
        # i.e. all of `num_iter`, whereas `SDF_MISS_DISTANCE` is passed after about 60, and the
        # exit of a bounding sphere after a few.
        t0 > t_max && return nothing
        # A non-finite probe, e.g. a `NaN` returned by `sdf`, fails `dist > eps` below and
        # would be misread as re-entering the surface. It is a miss, not a hit.
        (isfinite(dist) && isfinite(t0)) || return nothing

        if dist > eps
            escaped = true
        elseif escaped
            normal = normal3d(shape, pos)
            # Filter out false positive hits caused by numerical noise when leaving an SDF.
            if !(dot(dir, normal) > eps)
                return Intersection(t0, normal, shape)
            end
        end
    end
    return nothing
end

"""
    _raymarch_inside(object::AbstractSDF, pos, dir; num_iter=1000, dl=0.1)

Perform the ray marching algorithm if the starting pos is inside of `object`.
"""
function _raymarch_inside(object::AbstractSDF{S},
        pos::AbstractArray{R},
        dir::AbstractArray{R},
        num_iter = 1000,
        dl = Config.get_sdf_inside_step()) where {S, R}
    # this method assumes semi-concave objects, i.e. might fail depending on the choice of dl
    T = promote_type(S, R)
    t0::T = 0
    i = 1
    # march the ray a fixed distance dl until position is outside of sdf, since some sdfs are not exact on the inside
    while i <= num_iter
        pos = pos + dl * dir
        t0 += dl
        dist = sdf(object, pos)
        # once outside the sdf, fall back to _raymarch_outside
        if dist > 0
            intersection = _raymarch_outside(object, pos, -dir)
            if intersection === nothing
                break
            end
            intersection.t = t0 - intersection.t
            return intersection
        end
        i += 1
    end
    # return no intersection if too many iterations or actual miss occurs
    return nothing
end

"""
    intersect3d(object::AbstractSDF, ray::AbstractRay, t_max = SDF_MISS_DISTANCE)

Intersection algorithm for sdf based shapes. A ray that starts outside of the `object` is a miss
once it has travelled the path length `t_max` in [m] without a hit.
"""
function intersect3d(object::AbstractSDF, ray::AbstractRay, t_max = SDF_MISS_DISTANCE)
    pos = position(ray)
    dir = direction(ray)
    d = sdf(object, pos)
    h = Config.get_sdf_surface_threshold()
    # Test if outside of sdf, else inside
    if d > h
        return _raymarch_outside(object, pos, dir, t_max)
    end
    # The one-sided difference along the ray is used to determine inside/outside instead of the normal
    if sdf(object, pos + h * dir) ≤ d
        return _raymarch_inside(object, pos, dir)
    else
        return _raymarch_outside(object, pos, dir, t_max)
    end
    # Return no intersection else
    return nothing
end

"""
    intersect3d(sphere::SingleBoundingSphere, object::AbstractSDF, ray::AbstractRay)

Tests the `ray` against the bounding `sphere` of the `object`, see
`intersect3d(sphere, shape, ray)`. If the ray passes through the sphere, it is marched as in
`intersect3d(object, ray)`, but is a miss as soon as it has left the sphere, since the `object` lies
within. A ray that leaves the `object` is thus a miss after a few steps.

The radius is enlarged by the surface threshold of the marching algorithm
(`Config.get_sdf_surface_threshold()`), below which a ray counts as starting on the surface.
"""
function intersect3d(sphere::SingleBoundingSphere, object::AbstractSDF, ray::AbstractRay)
    t_out = _sphere_exit(sphere, ray, Config.get_sdf_surface_threshold())
    isnothing(t_out) && return nothing
    return intersect3d(object, ray, t_out)
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
