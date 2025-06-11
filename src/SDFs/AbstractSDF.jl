const eps_srf = 1e-9     # helps to decide if on SDF surface
const eps_ray = 1e-10     # helps to terminate outside ray marching routine
const eps_ins = 1e-0      # internal ray marching length step

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
        center, r = bounding_sphere(s)
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
normal3d(s::AbstractSDF, pos) = numeric_gradient(s, pos)

function numeric_gradient(s::AbstractSDF, pos)
    # approximate ∇ of s at pos
    eps = 1e-8
    norm = Point3(sdf(s, pos + Point3(eps, 0, 0)) - sdf(s, pos - Point3(eps, 0, 0)),
        sdf(s, pos + Point3(0, eps, 0)) - sdf(s, pos - Point3(0, eps, 0)),
        sdf(s, pos + Point3(0, 0, eps)) - sdf(s, pos - Point3(0, 0, eps)))
    return normalize(norm)
end

"""
    _raymarch_outside(shape::AbstractSDF, pos, dir; num_iter=1000, eps=1e-10)

Perform the ray marching algorithm if the starting pos is outside of `shape`.
"""
function _raymarch_outside(shape::AbstractSDF{S},
        pos::AbstractArray{R},
        dir::AbstractArray{R},
        num_iter = 1000,
        eps = eps_ray) where {S, R}
    T = promote_type(S, R)
    dist = sdf(shape, pos)
    t0 = dist
    i = 1
    # march the ray based on the last returned distance
    while i <= num_iter
        pos = pos + dist * dir
        dist = sdf(shape, pos)
        t0 += dist
        i += 1
        # surface has been reached if distance is less than tolerance (highly convex surfaces can require many iterations)
        if dist < eps
            normal = normal3d(shape, pos)
            return Intersection(t0, normal, shape)
        end
    end
    # return no intersection if tolerance has not been reached or actual miss occurs
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
        dl = eps_ins) where {S, R}
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
    intersect3d(sphere::AbstractSphere, ray::Ray)

Intersection algorithm for sdf based shapes.
"""
function intersect3d(object::AbstractSDF, ray::AbstractRay)
    pos = position(ray)
    dir = direction(ray)
    d = sdf(object, pos)
    # Test if outside of sdf, else inside
    if d > eps_srf
        return _raymarch_outside(object, pos, dir)
    end
    # Test if normal and ray dir oppose or align to determine if ray exits object
    n = normal3d(object, pos)
    if dot(dir, n) ≤ 0
        return _raymarch_inside(object, pos, dir)
    end
    # Return no intersection else
    return nothing
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
