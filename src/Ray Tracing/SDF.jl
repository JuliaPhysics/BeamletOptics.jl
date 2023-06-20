"""
    AbstractSDF <: AbstractShape

Provides a shape function based on signed distance functions. See https://iquilezles.org/articles/distfunctions/ for more information.

# Implementation reqs.
Subtypes of `AbstractSDF` should implement all reqs. of `AbstractShape` as well as the following:

# Functions
- `sdf(::AbstractSDF, point)`: a function that returns the signed distance for a point in 3D space
"""
abstract type AbstractSDF{T} <: AbstractShape{T} end

"""
    _world_to_sdf(sdf, point)

Transforms the coordinates of `point` into a reference frame where the `sdf` lies at the origin. Useful to represent translation and rotation.
If rotations are applied, the rotation is applied around the local sdf coordinate system.
"""
function _world_to_sdf(sdf::AbstractSDF, point)
    # transforms world coords to sdf coords
    T = transpose(orientation(sdf))
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
    x = LinRange(xmin - 1e-5, xmax + 1e-5, 50)
    y = LinRange(ymin - 1e-5, ymax + 1e-5, 50)
    z = LinRange(zmin - 1e-5, zmax + 1e-5, 50)
    sdf_values = Float32.([sdf(s, [i, j, k]) for i in x, j in y, k in z])
    mc = MC(sdf_values; x=Float32.(x), y=Float32.(y), z=Float32.(z))
    march(mc)
    vertices = transpose(reinterpret(reshape, Float32, mc.vertices))
    faces = transpose(reinterpret(reshape, Int64, mc.triangles))
    mesh!(axis, vertices, faces, transparency=true)
    return nothing
end

"""
    normal3d(s::AbstractSDF, pos)

Computes the normal vector of `s` at `pos`.
"""
function normal3d(s::AbstractSDF, pos)
    # approximate ∇ of s at pos
    eps = 1e-7
    norm = [
        sdf(s, pos + [eps, 0, 0]) - sdf(s, pos - [eps, 0, 0])
        sdf(s, pos + [0, eps, 0]) - sdf(s, pos - [0, eps, 0])
        sdf(s, pos + [0, 0, eps]) - sdf(s, pos - [0, 0, eps])
    ]
    normalize3d!(norm)
    return norm
end

"""
    _raymarch_outside(object::AbstractSDF, pos, dir; num_iter=1000, eps=1e-10)

Perform the ray marching algorithm if the starting pos is outside of `object`.
"""
function _raymarch_outside(object::AbstractSDF{S}, pos::AbstractArray{R}, dir::AbstractArray{R}; num_iter=1000, eps=1e-10) where {S,R}
    T = promote_type(S, R)
    dist = sdf(object, pos)
    t0 = dist
    i = 1
    # march the ray based on the last returned distance
    while i <= num_iter
        pos = pos + dist * dir
        dist = sdf(object, pos)
        t0 += dist
        i += 1
        # surface has been reached if distance is less than tolerance (highly convex surfaces can require many iterations)
        if dist < eps
            normal = normal3d(object, pos)
            return Intersection{T}(t0, normal, nothing)
        end
    end
    # return no intersection if tolerance has not been reached or actual miss occurs
    return nothing
end

"""
    _raymarch_inside(object::AbstractSDF, pos, dir; num_iter=1000, dl=0.1)

Perform the ray marching algorithm if the starting pos is inside of `object`.
"""
function _raymarch_inside(object::AbstractSDF{S}, pos::AbstractArray{R}, dir::AbstractArray{R}; num_iter=1000, dl=0.1) where {S,R}
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
    eps_srf = 1e-5      # helps to decide if on SDF surface
    eps_ray = 1e-10     # helps to terminate outside ray marching routine
    eps_ins = 1e-0      # internal ray marching length step
    dist = sdf(object, position(ray))
    # Test if outside of sdf, else inside
    if dist > eps_srf
        return _raymarch_outside(object, position(ray), direction(ray), eps=eps_ray)
    end
    # Test if normal and ray dir oppose or align to determine if ray exits object
    test = fast_dot3d(direction(ray), normal3d(object, position(ray))) > 0
    if test ≤ 0
        return _raymarch_inside(object, position(ray), direction(ray), dl=eps_ins)
    end
    # Return no intersection else
    return nothing
end

"""
    SphereSDF <: AbstractSDF

Implements sphere SDF. Orientation is fixed to unity matrix.
"""
struct SphereSDF{T} <: AbstractSDF{T}
    id::UUID
    pos::Vector{T}
    radius::T
end

orientation(::SphereSDF{T}) where {T} = Matrix{T}(I, 3, 3)
orientation!(::SphereSDF, ::Any) = nothing

function sdf(sphere::SphereSDF, point)
    p = _world_to_sdf(sphere, point)
    return norm3d(p) - sphere.radius
end

"""
    CylinderSDF <: AbstractSDF

Implements cylinder SDF. Cylinder is initially orientated along the y-axis and symmetrical in x-z.
"""
struct CylinderSDF{T} <: AbstractSDF{T}
    id::UUID
    dir::Matrix{T}
    pos::Vector{T}
    radius::T
    height::T
end

function sdf(cylinder::CylinderSDF, point)
    p = _world_to_sdf(cylinder, point)
    d = abs.([norm2d([p[1], p[3]]), p[2]]) - [cylinder.radius, cylinder.height]
    return min(maximum(d), 0) + norm2d(max.(d, 0))
end

"""
    CutSphereSDF <: AbstractSDF

Implements SDF of a sphere which is cut off in the x-z-plane at some point along the y-axis.
"""
struct CutSphereSDF{T} <: AbstractSDF{T}
    id::UUID
    dir::Matrix{T}
    pos::Vector{T}
    radius::T
    height::T
    w::T
end

"""
    CutSphereSDF(pos, radius, height)

Constructs a sphere with `radius` which is cut off along the y-axis at `height`.
"""
function CutSphereSDF(pos::Vector{V}, radius::R, height::H) where {V,R,H}
    if abs(height) ≥ radius
        error("Cut off height must be smaller than radius")
    end
    T = promote_type(V, R, H)
    w = sqrt(radius^2 - height^2)
    return CutSphereSDF{T}(uuid4(), Matrix{T}(I, 3, 3), pos, radius, height, w)
end

function sdf(cs::CutSphereSDF, point)
    p = _world_to_sdf(cs, point)
    q = [norm2d([p[1], p[3]]), p[2]]
    s = max((cs.height - cs.radius) * q[1]^2 + cs.w^2 * (cs.height + cs.radius - 2 * q[2]), cs.height * q[1] - cs.w * q[2])
    if s < 0
        return norm2d(q) - cs.radius
    elseif q[1] < cs.w
        return cs.height - q[2]
    else
        return norm2d(q - [cs.w, cs.height])
    end
end

"""
    ThinLensSDF <: AbstractSDF

Implements SDF of a spherical lens with which is composed of two joint cut spheres.
"""
struct ThinLensSDF{T} <: AbstractSDF{T}
    id::UUID
    dir::Matrix{T}
    pos::Vector{T}
    front::CutSphereSDF{T}
    back::CutSphereSDF{T}
end

sag(r::Real, l::Real) = r - sqrt(r^2 - 0.25 * l^2)

"""
    ThinLensSDF(r1, r2, d)

Constructs a thin lens with two spherical surfaces (`r1` and `r2`) of diameter `d`. Radii must be ≥ 0 and smaller than d/2!
"""
function ThinLensSDF(r1::L, r2::M, d::O) where {L,M,O}
    if 2r1 < d || 2r2 < d
        error("r1, r2 ≥ d/2 !")
    end
    T = promote_type(L, M, O)
    s1 = r1 - sag(r1, d)
    s2 = r2 - sag(r2, d)
    # Create cut spheres and shift by surface radius - sag
    front = CutSphereSDF([0, 0, 0.0], r1, s1)
    back = CutSphereSDF([0, 0, 0.0], r2, s2)
    translate3d!(front, [0, -s1, 0])
    # Rotate and move back cut sphere
    zrotate3d!(back, deg2rad(180))
    translate3d!(back, [0, s2, 0])
    return ThinLensSDF{T}(uuid4(), Matrix{T}(I, 3, 3), zeros(T, 3), front, back)
end

function sdf(tl::ThinLensSDF, pos)
    p = _world_to_sdf(tl, pos)
    return min(sdf(tl.front, p), sdf(tl.back, p))
end
