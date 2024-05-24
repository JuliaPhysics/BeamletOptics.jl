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
    render_object!(axis, s::AbstractSDF)

Render the surface of `s` based on the marching cubes algorithm.
"""
function render_object!(axis, s::AbstractSDF)
    # Calculate object limits from far away
    xmin = sdf(s, [-1000, 0, 0]) - 1000
    ymin = sdf(s, [0, -1000, 0]) - 1000
    zmin = sdf(s, [0, 0, -1000]) - 1000
    xmax = 1000 - sdf(s, [1000, 0, 0])
    ymax = 1000 - sdf(s, [0, 1000, 0])
    zmax = 1000 - sdf(s, [0, 0, 1000])
    x = LinRange(xmin - 1e-5, xmax + 1e-5, 100)
    y = LinRange(ymin - 1e-5, ymax + 1e-5, 100)
    z = LinRange(zmin - 1e-5, zmax + 1e-5, 100)
    sdf_values = Float32.([sdf(s, [i, j, k]) for i in x, j in y, k in z])
    mc = MC(sdf_values; x = Float32.(x), y = Float32.(y), z = Float32.(z))
    march(mc)
    vertices = transpose(reinterpret(reshape, Float32, mc.vertices))
    faces = transpose(reinterpret(reshape, Int64, mc.triangles))
    render_sdf_mesh!(axis, vertices, faces, transparency = true)
    return nothing
end
render_sdf_mesh!(::Any, vertices, faces; transparency = true) = nothing

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
        dir::AbstractArray{R};
        num_iter = 1000,
        eps = 1e-10) where {S, R}
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
        dir::AbstractArray{R};
        num_iter = 1000,
        dl = 0.1) where {S, R}
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
            if isnothing(intersection)
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
    eps_srf = 1e-9     # helps to decide if on SDF surface
    eps_ray = 1e-10     # helps to terminate outside ray marching routine
    eps_ins = 1e-0      # internal ray marching length step
    dist = sdf(object, position(ray))
    # Test if outside of sdf, else inside
    if dist > eps_srf
        return _raymarch_outside(object, position(ray), direction(ray), eps = eps_ray)
    end
    # Test if normal and ray dir oppose or align to determine if ray exits object
    test = dot(direction(ray), normal3d(object, position(ray))) > 0
    if test ≤ 0
        return _raymarch_inside(object, position(ray), direction(ray), dl = eps_ins)
    end
    # Return no intersection else
    return nothing
end
