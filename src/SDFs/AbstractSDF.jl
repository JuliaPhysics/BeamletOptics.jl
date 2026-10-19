
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
    _raymarch_inside(shape::AbstractSDF, pos, dir, t_max=SDF_MISS_DISTANCE, num_iter=10_000, eps=1e-10, h=1e-9)

Perform the ray marching algorithm if the starting pos is inside of `shape`, or on its surface with
`dir` pointing inside, and return the intersection where the ray leaves `shape` first. The ray is
marched by the distance `-sdf`, which is not larger than the distance to the surface for an SDF that
is exact or a lower bound inside, hence no step passes the exit. A point where the SDF is about 0 is
an exit only if the ray leaves the object there, not at an internal seam of a composite SDF. If a
step ends outside, e.g. for an SDF that overestimates the distance inside or after a minimum step,
the exit within this step is found within the bracket of this step, see `_bracketed_exit`. The ray is a miss once it has travelled the path length `t_max`
in [m] without leaving `shape`. A ray that grazes the surface of a thin part from inside takes
steps of about its depth, hence the many iterations. A ray that starts on the surface may start up to the surface
threshold `h` outside of it, see `intersect3d(object::AbstractSDF, ray)`.
"""
function _raymarch_inside(shape::AbstractSDF{S},
        pos::AbstractArray{R},
        dir::AbstractArray{R},
        t_max = SDF_MISS_DISTANCE,
        num_iter = 10_000,
        eps = Config.get_sdf_raymarch_eps(),
        h = Config.get_sdf_surface_threshold()) where {S, R}
    T = promote_type(S, R)
    t0 = zero(T)
    dist = sdf(shape, pos)
    # A ray that starts on the surface counts as inside once it is deeper than `eps`, before that a
    # point of the starting surface is no exit.
    entered = dist < -eps
    ω = T(1.6)
    for _ in 1:num_iter
        # The minimum step leaves the starting surface and the internal seams of composite SDFs, e.g.
        # the face between two parts of a `UnionSDF`, where the SDF is about 0 within the object. A
        # step that thereby passes the exit is caught below.
        t_prev, dist_prev = t0, dist
        step = max(-dist, eps + t0 * T(1e-3))
        # Over-relaxed step (Keinert et al., 2014): it is kept if the balls of radius -sdf around
        # both ends cover it, i.e. no surface lies in between, else the plain step is taken. The
        # factor grows while steps are kept, e.g. when the ray leaves the starting surface or a seam.
        t0 = t_prev + ω * step
        dist = sdf(shape, pos + t0 * dir)
        if dist ≤ 0 && -dist_prev - dist ≥ t0 - t_prev
            ω = min(2ω, T(32))
        else
            ω = T(1.6)
            t0 = t_prev + step
            dist = sdf(shape, pos + t0 * dir)
        end
        t0 > t_max && return nothing
        p = pos + t0 * dir
        # A non-finite probe, e.g. a `NaN` returned by `sdf`, is a miss
        (isfinite(dist) && isfinite(t0)) || return nothing
        # The normal at the exit, see `_exit_normal`
        if dist > eps && entered
            # The step ended outside: the exit lies within it
            t_exit = _bracketed_exit(shape, pos, dir, t_prev, dist_prev, t0, dist, eps)
            return Intersection(t_exit, _exit_normal(shape, pos, dir, t_exit, eps, h), shape)
        elseif dist > h
            # The ray left before it got deeper than `eps`, i.e. it grazes the surface
            return Intersection(t0, normal3d(shape, p), shape)
        elseif dist < -eps
            entered = true
        elseif entered
            dist_ahead = sdf(shape, p + h * dir)
            if dist_ahead > 0
                # On the surface and leaving it, in contrast to an internal seam. The SDF value is
                # not the distance for every SDF, hence the exit is refined within (t0, t0 + h].
                t_exit = _bracketed_exit(shape, pos, dir, t0, dist, t0 + h, dist_ahead, eps)
                return Intersection(t_exit, _exit_normal(shape, pos, dir, t_exit, eps, h), shape)
            end
        end
    end
    # return no intersection if too many iterations
    return nothing
end

# The outward normal where the ray from `pos` along `dir` leaves `shape` at the path length `t`. It is
# evaluated outside, `eps` behind the exit, like the normal of a hit of `_raymarch_outside`:
# - within a composite SDF, e.g. where the parts of a union overlap, the gradient may stem from a part
#   other than the one whose surface is left
# - exactly on the surface, the gradient of some SDFs is not defined, e.g. of a distance that is 0
# - as close as possible, since near thin features another surface may lie close behind the exit
# A normal that points against the ray, which leaves, stems from a part that touches the exit point
# from outside, e.g. at the vertex of a concave face on the back of a plate. It is then evaluated
# farther outside.
function _exit_normal(shape::AbstractSDF, pos, dir, t, eps, h)
    normal = normal3d(shape, pos + (t + eps) * dir)
    dot(normal, dir) > 0 && return normal
    return normal3d(shape, pos + (t + h) * dir)
end

# The path length at which the ray from `pos` along `dir` leaves `shape`, for the SDF values
# fa = sdf(a) ≤ 0 < fb = sdf(b) or fa ≤ eps < fb. The Illinois variant of the regula falsi keeps the
# bracket, i.e. the result is reliable also where the SDF is not exact, and converges in a few steps
# where the SDF is about linear, as near a surface. The result is the secant root of the final bracket
# of width ≤ `eps`: it changes continuously with the parameters of the shape and the ray, e.g. for
# finite differences and automatic differentiation, and it lies within `eps` of the surface, i.e. a
# ray that starts there and leaves the shape is not taken for a ray inside.
function _bracketed_exit(shape::AbstractSDF, pos, dir, a, fa, b, fb, eps)
    # the SDF values at the ends, `fa` and `fb` are scaled by the Illinois steps
    sa, sb = fa, fb
    side = 0
    for _ in 1:100
        b - a > eps || break
        m = (a * fb - b * fa) / (fb - fa)
        # fall back to bisection if the secant leaves the bracket, e.g. for a NaN
        (a < m < b) || (m = (a + b) / 2)
        fm = sdf(shape, pos + m * dir)
        if fm > 0
            b, fb, sb = m, fm, fm
            side == -1 && (fa /= 2)
            side = -1
        else
            a, fa, sa = m, fm, fm
            side == 1 && (fb /= 2)
            side = 1
        end
    end
    m = (a * sb - b * sa) / (sb - sa)
    return a ≤ m ≤ b ? m : (a + b) / 2
end

"""
    intersect3d(object::AbstractSDF, ray::AbstractRay, t_max = SDF_MISS_DISTANCE)

Intersection algorithm for sdf based shapes. A ray that starts outside of the `object` is a miss
once it has travelled the path length `t_max` in [m] without a hit. A ray that starts inside, or on
the surface pointing inside, hits the surface where it leaves the `object` first, see
`_raymarch_inside`.
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
    # On the surface, the one-sided difference along the ray determines whether the ray points
    # inside or outside, instead of the normal. Below the surface the ray is inside in any direction.
    if d < -h || sdf(object, pos + h * dir) ≤ d
        # within the object, the exit lies within any bounding sphere, hence no `t_max`
        return _raymarch_inside(object, pos, dir)
    else
        return _raymarch_outside(object, pos, dir, t_max)
    end
    # Return no intersection else
    return nothing
end

"""
    intersect3d(sphere::AbstractBoundingSphere, object::AbstractSDF, ray::AbstractRay)

Tests the `ray` against the bounding `sphere` of the `object`, see
`intersect3d(sphere, shape, ray)`. If the ray passes through the sphere, it is marched as in
`intersect3d(object, ray)`, but is a miss as soon as it has left the sphere, since the `object` lies
within. A ray that leaves the `object` is thus a miss after a few steps.

The radius is enlarged by the surface threshold of the marching algorithm
(`Config.get_sdf_surface_threshold()`), below which a ray counts as starting on the surface.
"""
function intersect3d(sphere::AbstractBoundingSphere, object::AbstractSDF, ray::AbstractRay)
    t_out = _sphere_exit(sphere, ray, Config.get_sdf_surface_threshold())
    isnothing(t_out) && return nothing
    return intersect3d(object, ray, t_out)
end

# Without a sphere the ray is marched up to the miss distance
intersect3d(::NoBoundingSphere, object::AbstractSDF, ray::AbstractRay) = intersect3d(object, ray)

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
