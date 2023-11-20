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
function normal3d(s::AbstractSDF, pos)
    # approximate ∇ of s at pos
    eps = 1e-7
    norm = Point3(sdf(s, pos + Point3(eps, 0, 0)) - sdf(s, pos - Point3(eps, 0, 0)),
        sdf(s, pos + Point3(0, eps, 0)) - sdf(s, pos - Point3(0, eps, 0)),
        sdf(s, pos + Point3(0, 0, eps)) - sdf(s, pos - Point3(0, 0, eps)))
    return normalize(norm)
end

"""
    _raymarch_outside(object::AbstractSDF, pos, dir; num_iter=1000, eps=1e-10)

Perform the ray marching algorithm if the starting pos is outside of `object`.
"""
function _raymarch_outside(object::AbstractSDF{S},
        pos::AbstractArray{R},
        dir::AbstractArray{R};
        num_iter = 1000,
        eps = 1e-10) where {S, R}
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

"""
    SphereSDF <: AbstractSDF

Implements sphere SDF. Orientation is fixed to unity matrix.
"""
mutable struct SphereSDF{T} <: AbstractSDF{T}
    const id::UUID
    pos::Point3{T}
    radius::T
end

SphereSDF(r::T) where {T} = SphereSDF{T}(uuid4(), zeros(T, 3), r)

orientation(::SphereSDF{T}) where {T} = Matrix{T}(I, 3, 3)
orientation!(::SphereSDF, ::Any) = nothing

function sdf(sphere::SphereSDF, point)
    p = _world_to_sdf(sphere, point)
    return norm(p) - sphere.radius
end

"""
    CylinderSDF <: AbstractSDF

Implements cylinder SDF. Cylinder is initially orientated along the y-axis and symmetrical in x-z.
"""
mutable struct CylinderSDF{T} <: AbstractSDF{T}
    const id::UUID
    dir::Matrix{T}
    pos::Point3{T}
    radius::T
    height::T
end

function CylinderSDF(r::R, h::H) where {R, H}
    T = promote_type(R, H)
    return CylinderSDF{T}(uuid4(), Matrix{T}(I, 3, 3), Point3{T}(0), r, h)
end

function sdf(cylinder::CylinderSDF, point)
    p = _world_to_sdf(cylinder, point)
    d = abs.(Point2(norm(Point2(p[1], p[3])), p[2])) -
        Point2(cylinder.radius, cylinder.height)
    return min(maximum(d), 0) + norm(max.(d, 0))
end

"""
    CutSphereSDF <: AbstractSDF

Implements SDF of a sphere which is cut off in the x-z-plane at some point along the y-axis.
"""
mutable struct CutSphereSDF{T} <: AbstractSDF{T}
    const id::UUID
    dir::Matrix{T}
    pos::Point3{T}
    radius::T
    height::T
    w::T
end

"""
    CutSphereSDF(pos, radius, height)

Constructs a sphere with `radius` which is cut off along the y-axis at `height`.
"""
function CutSphereSDF(radius::R, height::H) where {R, H}
    if abs(height) ≥ radius
        error("Cut off height must be smaller than radius")
    end
    T = promote_type(R, H)
    w = sqrt(radius^2 - height^2)
    return CutSphereSDF{T}(uuid4(), Matrix{T}(I, 3, 3), zeros(T, 3), radius, height, w)
end

function sdf(cs::CutSphereSDF, point)
    p = _world_to_sdf(cs, point)
    q = Point2(norm([p[1], p[3]]), p[2])
    s = max((cs.height - cs.radius) * q[1]^2 + cs.w^2 * (cs.height + cs.radius - 2 * q[2]),
        cs.height * q[1] - cs.w * q[2])
    if s < 0
        return norm(q) - cs.radius
    elseif q[1] < cs.w
        return cs.height - q[2]
    else
        return norm(q - Point2(cs.w, cs.height))
    end
end

"""
    AbstractSphericalLensSDF{T} <: AbstractSDF{T}

An abstract type for SDFs which represent spherical lenses, i.e. biconvex or plano-concave.
"""
abstract type AbstractSphericalLensSDF{T} <: AbstractSDF{T} end

"""
    sag(r::Real, l::Real)

Calculates the sag of a cut circle with radius `r` and chord length `l`
"""
sag(r::Real, l::Real) = r - sqrt(r^2 - 0.25 * l^2)

check_sag(r, d) = 2 * abs(r) < d ? error("r=$(r) must be ≥ d/2, d=$(d)!") : nothing

"""
    ThinLensSDF <: AbstractSphericalLensSDF

Implements the SDF of an ideal spherical lens which is composed of two joint cut spheres.
"""
mutable struct ThinLensSDF{T} <: AbstractSphericalLensSDF{T}
    const id::UUID
    dir::Matrix{T}
    pos::Point3{T}
    front::CutSphereSDF{T}
    back::CutSphereSDF{T}
end

"""
    ThinLensSDF(r1, r2, d=1inch)

Constructs a cylindrical bi-convex thin lens SDF with:

- `r1 > 0`: radius of convex front
- `r2 > 0`: radius of convex back
- `d`: lens diameter, default value is one inch

The spherical surfaces are constructed flush.
"""
function ThinLensSDF(r1::L, r2::M, d::O = 1inch) where {L, M, O}
    check_sag(r1, d)
    check_sag(r2, d)
    T = promote_type(L, M, O)
    s1 = r1 - sag(r1, d)
    s2 = r2 - sag(r2, d)
    # Create cut spheres and shift by surface radius - sag
    front = CutSphereSDF(r1, s1)
    back = CutSphereSDF(r2, s2)
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

"""
    BiConvexLensSDF <: AbstractSphericalLensSDF

Implements a cylindrical lens SDF with two convex surfaces and a cylindrical mid section.
"""
mutable struct BiConvexLensSDF{T} <: AbstractSphericalLensSDF{T}
    const id::UUID
    dir::Matrix{T}
    pos::Point3{T}
    front::CutSphereSDF{T}
    back::CutSphereSDF{T}
    mid::CylinderSDF{T}
end

"""
    BiConvexLensSDF(r1, r2, l, d=1inch)

Constructs a cylindrical bi-convex lens SDF with:

- `r1` > 0: radius of convex front
- `r2` > 0: radius of convex back
- `l`: lens thickness
- `d`: lens diameter, default value is one inch

The spherical surfaces are constructed flush with the cylinder surface.
"""
function BiConvexLensSDF(r1::L, r2::M, l::N, d::O = 1inch) where {L, M, N, O}
    T = promote_type(L, M, N, O)
    check_sag(r1, d)
    check_sag(r2, d)
    s1 = sag(r1, d)
    s2 = sag(r2, d)
    # Calculate length of cylindrical section
    l = l - (s1 + s2)
    s1 = r1 - s1
    s2 = r2 - s2
    front = CutSphereSDF(r1, s1)
    back = CutSphereSDF(r2, s2)
    mid = CylinderSDF(d / 2, l / 2)
    # Shift and rotate cut spheres into position
    translate3d!(front, [0, -s1 + l / 2, 0])
    zrotate3d!(back, deg2rad(180))
    translate3d!(back, [0, s2 - l / 2, 0])
    return BiConvexLensSDF{T}(uuid4(), Matrix{T}(I, 3, 3), zeros(T, 3), front, back, mid)
end

function sdf(bcl::BiConvexLensSDF, pos)
    p = _world_to_sdf(bcl, pos)
    return min(sdf(bcl.front, p), sdf(bcl.mid, p), sdf(bcl.back, p))
end

"""
    BiConcaveLensSDF <: AbstractSphericalLensSDF

Implements a cylindrical lens SDF with two concave surfaces and a cylindrical mid section.
"""
mutable struct BiConcaveLensSDF{T} <: AbstractSphericalLensSDF{T}
    const id::UUID
    dir::Matrix{T}
    pos::Point3{T}
    front::SphereSDF{T}
    back::SphereSDF{T}
    mid::CylinderSDF{T}
end

"""
    BiConcaveLensSDF(r1, r2, l, d=1inch)

Constructs a bi-concave lens SDF with:

- `r1` > 0: radius of concave front
- `r2` > 0: radius of convex back
- `l`: lens thickness
- `d`: lens diameter, default value is one inch

The spherical surfaces are constructed flush with the cylinder surface.
"""
function BiConcaveLensSDF(r1::L, r2::M, l::N, d::O = 1inch) where {L, M, N, O}
    T = promote_type(L, M, N, O)
    check_sag(r1, d)
    check_sag(r2, d)
    s1 = sag(r1, d)
    s2 = sag(r2, d)
    # Calculate length of cylindrical section
    l = l + s1 + s2
    if s1 + s2 ≥ l
        error("r1=$(r1), r2=$(r2) results in hollow lens")
    end
    front = SphereSDF(r1)
    back = SphereSDF(r2)
    mid = CylinderSDF(d / 2, l / 2)
    # Shift and rotate subtraction spheres into position
    translate3d!(front, [0, (r1 + l / 2 - s1), 0])
    translate3d!(back, [0, -(r2 + l / 2 - s2), 0])
    return BiConcaveLensSDF{T}(uuid4(), Matrix{T}(I, 3, 3), zeros(T, 3), front, back, mid)
end

function sdf(bcl::BiConcaveLensSDF, pos)
    p = _world_to_sdf(bcl, pos)
    return max(-sdf(bcl.front, p), sdf(bcl.mid, p), -sdf(bcl.back, p))
end

"""
    PlanoConvexLensSDF <: AbstractSphericalLensSDF

Implements a cylindrical lens SDF with one convex and one planar surface.
"""
mutable struct PlanoConvexLensSDF{T} <: AbstractSphericalLensSDF{T}
    const id::UUID
    dir::Matrix{T}
    pos::Point3{T}
    front::CutSphereSDF{T}
    back::CylinderSDF{T}
end

"""
    PlanoConvexLensSDF(r, l, d=1inch)

Constructs a plano-convex lens SDF with:

- `r` > 0: front radius
- `l`: lens thickness
- `d`: lens diameter, default value is one inch

The spherical surface is constructed flush with the cylinder surface.
"""
function PlanoConvexLensSDF(r::R, l::L, d::D = 1inch) where {R, L, D}
    T = promote_type(R, L, D)
    check_sag(r, d)
    s = sag(r, d)
    # Calculate length of cylindrical section
    l = l - s
    s = r - s
    front = CutSphereSDF(r, s)
    back = CylinderSDF(d / 2, l / 2)
    # Shift and rotate cut spheres into position
    translate3d!(front, [0, -s + l / 2, 0])
    return PlanoConvexLensSDF{T}(uuid4(), Matrix{T}(I, 3, 3), zeros(T, 3), front, back)
end

function sdf(pcl::PlanoConvexLensSDF, pos)
    p = _world_to_sdf(pcl, pos)
    return min(sdf(pcl.front, p), sdf(pcl.back, p))
end

"""
    PlanoConcaveLensSDF <: AbstractSphericalLensSDF

Implements a cylindrical lens SDF with one concave and one planar surface.
"""
mutable struct PlanoConcaveLensSDF{T} <: AbstractSphericalLensSDF{T}
    const id::UUID
    dir::Matrix{T}
    pos::Point3{T}
    front::SphereSDF{T}
    back::CylinderSDF{T}
end

"""
    PlanoConcaveLensSDF(r, l, d=1inch)

Constructs a plano-concave lens SDF with:

- `r` > 0: front radius
- `l`: lens thickness
- `d`: lens diameter, default value is one inch

The spherical surface is constructed flush with the cylinder surface.
"""
function PlanoConcaveLensSDF(r::R, l::L, d::D = 1inch) where {R, L, D}
    T = promote_type(R, L, D)
    check_sag(r, d)
    s = sag(r, d)
    # Calculate length of cylindrical section
    l = l + s
    if s ≥ l
        error("r=$(r) results in hollow lens")
    end
    front = SphereSDF(r)
    back = CylinderSDF(d / 2, l / 2)
    # Shift and rotate subtraction sphere into position
    translate3d!(front, [0, (r + l / 2 - s), 0])
    return PlanoConcaveLensSDF{T}(uuid4(), Matrix{T}(I, 3, 3), zeros(T, 3), front, back)
end

function sdf(pcl::PlanoConcaveLensSDF, pos)
    p = _world_to_sdf(pcl, pos)
    return max(-sdf(pcl.front, p), sdf(pcl.back, p))
end

"""
    ConvexConcaveLensSDF <: AbstractSphericalLensSDF

Implements a cylindrical lens SDF with one convex and one concave surface, as well as a mid section.
"""
mutable struct ConvexConcaveLensSDF{T} <: AbstractSphericalLensSDF{T}
    const id::UUID
    dir::Matrix{T}
    pos::Point3{T}
    front::CutSphereSDF{T}
    back::SphereSDF{T}
    mid::CylinderSDF{T}
end

"""
    ConvexConcaveLensSDF(r1, r2, l, d=1inch)

Constructs a positive/negative meniscus lens SDF with:

- `r1` > 0: radius of convex front
- `r2` > 0: radius of concave back
- `l`: lens thickness
- `d`: lens diameter, default value is one inch

The spherical surface is constructed flush with the cylinder surface.
"""
function ConvexConcaveLensSDF(r1::L, r2::M, l::N, d::O = 1inch) where {L, M, N, O}
    T = promote_type(L, M, N, O)
    check_sag(r1, d)
    check_sag(r2, d)
    # Calculate convex and concave sag
    s1 = sag(r1, d)
    s2 = sag(r2, d)
    # Calculate length of cylindrical section
    l = l - s1 + s2
    s1 = r1 - s1
    if s2 ≥ l + s1
        error("r1=$(r1), r2=$(r2) results in hollow lens")
    end
    front = CutSphereSDF(r1, s1)
    back = SphereSDF(r2)
    mid = CylinderSDF(d / 2, l / 2)
    # Shift and rotate subtraction spheres into position
    translate3d!(front, [0, (-s1 + l / 2), 0])
    translate3d!(back, [0, -(r2 + l / 2 - s2), 0])
    return ConvexConcaveLensSDF{T}(uuid4(),
        Matrix{T}(I, 3, 3),
        zeros(T, 3),
        front,
        back,
        mid)
end

function sdf(ccl::ConvexConcaveLensSDF, pos)
    p = _world_to_sdf(ccl, pos)
    return max(min(sdf(ccl.front, p), sdf(ccl.mid, p)), -sdf(ccl.back, p))
end
