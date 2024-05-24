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

"""
    SphereSDF <: AbstractSDF

Implements sphere SDF. Orientation is fixed to unity matrix.
"""
mutable struct SphereSDF{T} <: AbstractSDF{T}
    pos::Point3{T}
    radius::T
end

SphereSDF(r::T) where {T} = SphereSDF{T}(zeros(T, 3), r)

orientation(::SphereSDF{T}) where {T} = SArray{Tuple{3,3}, T}(I)
transposed_orientation(::SphereSDF{T}) where {T} = SArray{Tuple{3,3}, T}(I)
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
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    radius::T
    height::T
end

function CylinderSDF(r::R, h::H) where {R, H}
    T = promote_type(R, H)
    return CylinderSDF{T}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        Point3{T}(0),
        r,
        h)
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
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
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
    return CutSphereSDF{T}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        zeros(T, 3),
        radius,
        height,
        w)
end

function sdf(cs::CutSphereSDF, point)
    p = _world_to_sdf(cs, point)
    q = Point2(norm(Point2(p[1], p[3])), p[2])
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
    RingSDF <: AbstractSDF

Implements the SDF of a ring in the x-z-plane for some distance in the y axis.
This allows to add planar outer sections to any SDF which fits inside of the ring.
"""
mutable struct RingSDF{T} <: AbstractSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    inner_radius::T
    hwidth::T
    hthickness::T
end

"""
    RingSDF(inner_radius, width, thickness)

Constructs a ring with `inner_radius` with a `width` and some thickness.
"""
function RingSDF(inner_radius::R, width::W, thickness::T) where {R, W, T}
    TT = promote_type(R, W, T)
    return RingSDF{TT}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        zeros(T, 3),
        inner_radius + width / 2,
        width / 2,
        thickness / 2)
end

function sdf_box(p, b)
    d = abs.(p) - b
    return norm(max.(d, zero(eltype(p)))) + min(max(d[1], d[2]), zero(eltype(p)))
end

function sdf(ring::RingSDF, point)
    p = _world_to_sdf(ring, point)

    return sdf_box(Point2(norm(Point2(p[1], p[3]))- ring.inner_radius, p[2]), Point2(ring.hwidth, ring.hthickness))
end

"""
    UnionSDF{T, TT <: Tuple} <: AbstractSDF{T}

This SDF represents the merging of two or more SDFs. If the constituent SDFs do not overlap
(they can and should touch) the resulting SDF should be still exact if the constituent SDFs
are exact.

The intended way to construct these is not explicitely but by just adding two `AbstractSDFs`
using the regular `+` operator.

```@example
s1 = SphereSDF(1.0)
translate3d!(s1, Point3(0, 1.0, 0.0))

s2 = SphereSDF(1.0)

# will result in a SDF with two spheres touching each other.
s_merged = s1 + s2
```

"""
mutable struct UnionSDF{T, TT <: Tuple} <: AbstractSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    sdfs::TT
end

function UnionSDF{T}(sdfs::Vararg{AbstractSDF{T}, N}) where {T, N}
    UnionSDF{T, typeof(sdfs)}(
        SMatrix{3,3}(one(T)*I),
        SMatrix{3,3}(one(T)*I),
        Point3{T}(zero(T)),
        sdfs
        )
end

function sdf(s::UnionSDF, pos)
    # sdf to world transform handled by sub-SDFs
    return minimum(sdf(_sdf, pos) for _sdf in s.sdfs)
end

Base.:+(s1::AbstractSDF{T}, s2::AbstractSDF{T}) where T = UnionSDF{T}(s1, s2)
Base.:+(union::UnionSDF{T}, sdf::AbstractSDF{T}) where T = UnionSDF{T}(union.sdfs..., sdf)
Base.:+(sdf::AbstractSDF{T}, union::UnionSDF{T}) where T = UnionSDF{T}(sdf, union.sdfs...)
Base.:+(u1::UnionSDF{T}, u2::UnionSDF{T}) where T = UnionSDF{T}(u1.sdfs..., u2.sdfs...)

function translate3d!(u::UnionSDF, offset)
    position!(u, position(u) .+ offset)
    translate3d!.(u.sdfs, Ref(offset))
    return nothing
end

function rotate3d!(u::UnionSDF, axis, θ)
    R = rotate3d(axis, θ)
    # Update group orientation
    orientation!(u, orientation(u) * R)
    # Rotate all sub-SDFs around union center
    for s in u.sdfs
        rotate3d!(s, axis, θ)
        v = position(s) - position(u)
        # Translate group around pivot point
        v = (R * v) - v
        translate3d!(s, v)
    end
    return nothing
end

abstract type RotationallySymmetricLensSDF{T} <: AbstractSDF{T} end

"""
    AbstractSphericalLensSDF{T} <: AbstractSDF{T}

An abstract type for SDFs which represent spherical lenses, i.e. biconvex or plano-concave.
"""
abstract type AbstractSphericalLensSDF{T} <: RotationallySymmetricLensSDF{T} end

"""
    sag(r::Real, l::Real)

Calculates the sag of a cut circle with radius `r` and chord length `l`
"""
sag(r::Real, l::Real) = r - sqrt(r^2 - 0.25 * l^2)

check_sag(r, d) = 2 * abs(r) < d ? error("r=$(r) must be ≥ d/2, d=$(d)!") : nothing

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
    return (front + back)
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
    return (front + mid + back)
end

"""
    BiConcaveLensSDF <: AbstractSphericalLensSDF

Implements a cylindrical lens SDF with two concave surfaces and a cylindrical mid section.
"""
mutable struct BiConcaveLensSDF{T} <: AbstractSphericalLensSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
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
- `md`: mechanical lens diameter, defaults to be identical to the lens diameter, Otherwise
        an outer ring section will be added to the lens, if `md` > `d`.

The spherical surfaces are constructed flush with the cylinder surface.
"""
function BiConcaveLensSDF(r1::L, r2::M, l::N, d::O = 1inch, md::MD = d) where {L, M, N, O, MD}
    T = promote_type(L, M, N, O, MD)
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

    lens = BiConcaveLensSDF{T}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        zeros(T, 3),
        front,
        back,
        mid)

    if md > d
        # add an outer ring
        ring = RingSDF(d/2, (md - d) / 2, l)
        lens += ring
    end

    return lens
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
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
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
    return (front + back)
end

"""
    PlanoConcaveLensSDF <: AbstractSphericalLensSDF

Implements a cylindrical lens SDF with one concave and one planar surface.
"""
mutable struct PlanoConcaveLensSDF{T} <: AbstractSphericalLensSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
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
- `md`: mechanical lens diameter, defaults to be identical to the lens diameter, Otherwise
        an outer ring section will be added to the lens, if `md` > `d`.

The spherical surface is constructed flush with the cylinder surface.
"""
function PlanoConcaveLensSDF(r::R, l::L, d::D = 1inch, md::MD = d) where {R, L, D, MD}
    T = promote_type(R, L, D, MD)
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
    lens = PlanoConcaveLensSDF{T}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        zeros(T, 3),
        front,
        back)

    if md > d
        # add an outer ring
        ring = RingSDF(d/2, (md - d) / 2, l)
        lens += ring
    end

    return lens
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
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
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
- `md`: mechanical lens diameter, defaults to be identical to the lens diameter, Otherwise
        an outer ring section will be added to the lens, if `md` > `d`.

The spherical surface is constructed flush with the cylinder surface.
"""
function ConvexConcaveLensSDF(r1::L, r2::M, l::N, d::O = 1inch, md::MD = d) where {L, M, N, O, MD}
    T = promote_type(L, M, N, O, MD)
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
    lens = ConvexConcaveLensSDF{T}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        zeros(T, 3),
        front,
        back,
        mid)

    if md > d
        # add an outer ring
        ring = RingSDF(d/2, (md - d) / 2, l)
        lens += ring
    end

    return lens
end

function sdf(ccl::ConvexConcaveLensSDF, pos)
    p = _world_to_sdf(ccl, pos)
    return max(min(sdf(ccl.front, p), sdf(ccl.mid, p)), -sdf(ccl.back, p))
end

abstract type AbstractAsphericalSurfaceSDF{T} <: AbstractSDF{T} end

"""
    ConvexAsphericalSurfaceSDF(r1, r2, l, d=1inch)

Constructs an aspheric lens with a convex-like surface according to ISO10110.

!!! note
    Currently, it is assumed that the aspheric surface is convex if the radius is positive.
    There might be unexpected effects for complex shapes which do not show a generally convex
    behavior.

- `coefficients` : (even) coefficients of the asphere.
- `radius` : radius of the lens
- `conic_constant` : conic constant of the lens surface
- `diameter`: lens diameter
"""
mutable struct ConvexAsphericalSurfaceSDF{T} <: AbstractAsphericalSurfaceSDF{T}
    coefficients::Vector{T}
    radius::T
    conic_constant::T
    diameter::T
    pos::Point3{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
end

# Constructor for ConvexAsphericalSurfaceSDF
function ConvexAsphericalSurfaceSDF(coefficients::Vector{T}, radius::T, conic_constant::T, diameter::T) where {T}
    return ConvexAsphericalSurfaceSDF{T}(coefficients, radius, conic_constant, diameter, Point3{T}(0), Matrix{T}(I, 3, 3), Matrix{T}(I, 3, 3))
end

"""
    ConcaveAsphericalSurfaceSDF(r1, r2, l, d=1inch)

Constructs an aspheric lens with a concave-like surface according to ISO10110.

!!! note
    Currently, it is assumed that the aspheric surface is concave if the radius is negative.
    There might be unexpected effects for complex shapes which do not show a generally concave
    behavior.

- `coefficients` : (even) coefficients of the asphere.
- `radius` : radius of the lens (negative!)
- `conic_constant` : conic constant of the lens surface
- `diameter`: lens diameter
- `mechanical_diameter`: mechanical lens diameter, defaults to be identical to the lens diameter, Otherwise
        an outer ring section will be added to the lens, if `mechanical_diameter` > `diameter`.
"""
mutable struct ConcaveAsphericalSurfaceSDF{T} <: AbstractAsphericalSurfaceSDF{T}
    coefficients::Vector{T}
    radius::T
    conic_constant::T
    diameter::T
    mechanical_diameter::T
    pos::Point3{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
end

# Constructor for ConcaveAsphericalSurfaceSDF
function ConcaveAsphericalSurfaceSDF(coefficients::V, radius::T, conic_constant::T, diameter::T, mechanical_diameter::T = diameter) where {T, V <: AbstractVector{T}}
    return ConcaveAsphericalSurfaceSDF(coefficients, radius, conic_constant, diameter, mechanical_diameter, Point3{T}(0), SMatrix{3,3}(one(T)*I), SMatrix{3,3}(one(T)*I))
end

"""
aspheric_equation(r, c, k, α_coeffs)

The aspheric surface equation. The asphere is defined by:
- `c` : The curvature (1/radius) of the surface
- `k` : The conic constant of the surface
- `α_coeffs` : The (even) aspheric coefficients, starting with A4.

This function returns NaN if the square root argument becomes negative.

!!! note
Only even aspheres are implemented at the moment. This will change soon.

"""
function aspheric_equation(r, c, k, α_coeffs)
    r2 = r^2
    sqrt_arg = 1 - (1 + k) * c^2 * r2
    if sqrt_arg < 0
        return NaN  # Handle as desired
    end
    sum_α = sum(i->α_coeffs[i] * r2^(i), eachindex(α_coeffs))
    return c * r2 / (1 + sqrt(sqrt_arg)) + sum_α
end

aspheric_equation(r::Real, a::AbstractAsphericalSurfaceSDF) = aspheric_equation(r, 1/a.radius, a.conic_constant, a.coefficients)

function gradient_aspheric_equation(r, c, k, α_coeffs)
    Ri = 1 / c
    sqrt_arg = 1-r^2*(1+k)/Ri^2
    sqrt_arg < 0 && return NaN
    gr = 2*r / (Ri*(√(sqrt_arg) + 1)) + r^3*(1+k)/(Ri^3*√(sqrt_arg)*(√(sqrt_arg)+1)^2)
    sum_r = sum(m->2*(m)*α_coeffs[m]*r^(2(m-1)+1), eachindex(α_coeffs))

    return Point2(-sum_r - gr, 1)
end

"""
    op_revolve(p, sdf2d::Function, offset)

Calculates the SDF at point `p` for the given 2D-SDF function with `offset` by revolving
the 2D shape around the z-axis.

"""
function op_revolve(p::Point3{T}, sdf2d::Function, offset=zero(T)) where {T <: Real}
    q = Point2(norm(Point2(p[1], p[2])) - offset, p[3])
    return sdf2d(q)
end

"""
    sd_line_segment(p, a, b)

Returns the signed distance from point `p` to the line segment described by the points `a`
and `b`.
"""
function sd_line_segment(p, a, b)
    pa, ba = p - a, b - a
    h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0)

    return norm(pa - h * ba)
end

"""
    convex_aspheric_surface_distance(r, z, c, k, d, α_coeffs)

Calculates the 2D distance field for an aspheric surface at radius `r` away from the optical axis
position `z`. The asphere is defined by:
- `c` : The curvature (1/radius) of the surface
- `k` : The conic constant of the surface
- `d` : The diameter of the asphere
- `α_coeffs` : The (even) aspheric coefficients, starting with A2.

Note that this is not just an infinite aspheric surface and also not a surface segment but
a closed 2D perimeter.

It is intended to pair the SDF derived from this distance field with a cylinder SDF to build a real lens.
"""
function convex_aspheric_surface_distance(r, z, c, k, d, α_coeffs)
    r2 = r^2
    r2_bound = (d/2)^2
    z_aspheric_value = aspheric_equation(r, c, k, α_coeffs)
    grad_z = gradient_aspheric_equation(r, c, k, α_coeffs)
    is_concave = c < 0

    z_aspheric_boundary = aspheric_equation(d/2, c, k, α_coeffs)
    grad_z_boundary = gradient_aspheric_equation(d/2, c, k, α_coeffs)

    # Handling NaN for points outside the aspheric surface
    if isnan(z_aspheric_value) || isnan(grad_z) || r2 > r2_bound
        distance_to_boundary = sqrt((r - sign(r)*d/2)^2 + (z - z_aspheric_boundary)^2)
        return distance_to_boundary / norm(grad_z_boundary)
    end

    # Calculate distance to the aspheric surface
    distance_to_aspheric = abs(z - z_aspheric_value) / norm(grad_z)
    # If we are inside the aperture, let's close the perimeter with a line segment
    a, b = Point2(d/2, z_aspheric_boundary), Point2(-d/2, z_aspheric_boundary)
    sdl = sd_line_segment(Point2(r, z), a, b) / norm(grad_z_boundary)
    if sign(c)*z_aspheric_value < sign(c)*z < sign(c)*z_aspheric_boundary
        return -min(sdl, distance_to_aspheric)
    else
        return min(sdl, distance_to_aspheric)
    end
end

function concave_aspheric_surface_distance(r, z, c, k, d, α_coeffs)
    r2 = r^2
    r2_bound = (d/2)^2
    z_aspheric_value = aspheric_equation(r, c, k, α_coeffs)
    grad_z = gradient_aspheric_equation(r, c, k, α_coeffs)

    z_aspheric_boundary = aspheric_equation(d/2, c, k, α_coeffs)
    grad_z_boundary = gradient_aspheric_equation(d/2, c, k, α_coeffs)

    # Handling NaN for points outside the aspheric surface
    if isnan(z_aspheric_value) || isnan(grad_z)
        distance_to_boundary = sqrt((r - sign(r)*d/2)^2 + (z - z_aspheric_boundary)^2)
        return distance_to_boundary / norm(grad_z_boundary)
    end

    distance_to_aspheric = abs(z - z_aspheric_value) / norm(grad_z)
    p = Point2(r, z)

    a1, b1 = Point2(d/2, z_aspheric_boundary), Point2(d/2, 0.0)
    sdl1 = sd_line_segment(p, a1, b1) / norm(grad_z_boundary)
    a2, b2 = Point2(d/2, 0.0), Point2(-d/2, 0.0)
    sdl2 = sd_line_segment(p, a2, b2) / norm(grad_z_boundary)
    a3, b3 = Point2(-d/2, 0.0), Point2(-d/2, z_aspheric_boundary)
    sdl3 = sd_line_segment(p, a3, b3) / norm(grad_z_boundary)

    if r2 > r2_bound
        # we are outside of the aperture, so the shortest distance can be determined
        # by the perimeter line segments which include the boundary points of the asphere
        return min(sdl1, sdl2, sdl3)
    else
        # return a negative sign within the asphere
        if z_aspheric_value < z < 0.0
            return -min(distance_to_aspheric, sdl1, sdl2, sdl3)
        else
            return min(distance_to_aspheric, sdl1, sdl2, sdl3)
        end
    end
end

function sdf(surface::ConvexAsphericalSurfaceSDF{T}, point) where {T}
    p_local = _world_to_sdf(surface, point)
    # standard logic is to have the z-axis as optical axis, so the aspheric code is written
    # with that convention in mind so everything matches ISO10110/textbook definitions.
    # We reinterpret the coordinates here, so xz is the transversal direction and y is the
    # optical axis.
    _pp = Point3{T}(p_local[1], p_local[3], p_local[2]) # xzy
    # rotate 2D sdf around the optical axis
    sdf_v = op_revolve(_pp,
        x->convex_aspheric_surface_distance(
            x[1],
            x[2],
            1/surface.radius,
            surface.conic_constant,
            surface.diameter,
            surface.coefficients
            ), zero(T))
    return sdf_v
end

function sdf(surface::ConcaveAsphericalSurfaceSDF{T}, point) where {T}
    p_local = _world_to_sdf(surface, point)
    # standard logic is to have the z-axis as optical axis, so the aspheric code is written
    # with that convention in mind so everything matches ISO10110/textbook definitions.
    # We reinterpret the coordinates here, so xz is the transversal direction and y is the
    # optical axis.
    _pp = Point3{T}(p_local[1], p_local[3], p_local[2]) # xzy
    # rotate 2D sdf around the optical axis
    sdf_v = op_revolve(_pp,
        x->concave_aspheric_surface_distance(
            x[1],
            x[2],
            1/surface.radius,
            surface.conic_constant,
            surface.diameter,
            surface.coefficients
            ), zero(T))
    return sdf_v
end

function render_object!(axis, asp::AbstractAsphericalSurfaceSDF; color=:red)
    radius = asp.diameter/2
    v = LinRange(0, 2π, 100)
    r = LinRange(1e-12, radius, 50)
    # Calculate beam surface at origin along y-axis, swap w and u
    y = aspheric_equation.(r, Ref(asp))
    u = y
    w = collect(r)
    if isa(asp, ConvexAsphericalSurfaceSDF)
        push!(u, u[end])
        push!(w, 1e-12)
    elseif isa(asp, ConcaveAsphericalSurfaceSDF)
        push!(u, 0, 0)
        push!(w, radius, 1e-12)
    else
        @warn "No suitable render fct. for $(typeof(asp))"
        return nothing
    end
    X = [w[i] * cos(v) for (i, u) in enumerate(u), v in v]
    Y = [u for u in u, v in v]
    Z = [w[i] * sin(v) for (i, u) in enumerate(u), v in v]
    # Transform into global coords
    R = asp.dir
    P = asp.pos
    Xt = R[1, 1] * X + R[1, 2] * Y + R[1, 3] * Z .+ P[1]
    Yt = R[2, 1] * X + R[2, 2] * Y + R[2, 3] * Z .+ P[2]
    Zt = R[3, 1] * X + R[3, 2] * Y + R[3, 3] * Z .+ P[3]
    render_surface!(axis, Xt, Yt, Zt; transparency = true, colormap = [color, color])
    return nothing
end

"""
    PlanoConvexAsphericalLensSDF(r, l, d=1inch)

Constructs a plano-convex aspheric lens SDF with:

- `r` > 0: front radius
- `l`: lens thickness
- `d`: lens diameter
- `k` : The conic constant of the surface
- `α_coeffs` : The (even) aspheric coefficients, starting with A4.

The spherical surface is constructed flush with the cylinder surface.
"""
function PlanoConvexAsphericalLensSDF(r::R, l::L, d::D, k::K, α_coeffs::AbstractVector{A}) where {R, L, D, K, A}
    T = promote_type(R, L, D, K, A)
    s = aspheric_equation(d/2, 1/r, k, α_coeffs)
    # Calculate length of cylindrical section
    l < 0 && error("Specified thickness is shorter than the lens sagitta at the edge")
    front = ConvexAsphericalSurfaceSDF(α_coeffs, r, k, d)
    back = CylinderSDF(d / 2, (l - abs(s)) / 2)
    # Shift and rotate cut spheres into position
    translate3d!(front, [0, -sign(r)*(l / 2 + abs(s) / 2), 0])
    return (front + back)
end

"""
    PlanoConcaveAsphericalLensSDF(r, l, d=1inch)

Constructs a plano-concave aspheric lens SDF with:

- `r` > 0: front radius
- `l`: lens thickness
- `d`: lens diameter
- `cz`: aspheric surface chip zone
- `k` : The conic constant of the surface
- `α_coeffs` : The (even) aspheric coefficients, starting with A4.
- `md`: lens mechanical diameter (default: md = d)

The spherical surface is constructed flush with the cylinder surface.
"""
function PlanoConcaveAsphericalLensSDF(r::R, l::L, d::D, k::K, α_coeffs::AbstractVector{A}, md::MD = d) where {R, L, D, K, A, MD}
    s = aspheric_equation(d/2, 1/r, k, α_coeffs)
    # Calculate length of cylindrical section
    l < 0 && error("Specified thickness is shorter than the lens sagitta at the edge")
    front = ConcaveAsphericalSurfaceSDF(α_coeffs, r, k, d)
    # Shift and rotate asphere into position
    translate3d!(front, [0, -(l - s)/2, 0])
    # add outer planar ring, if required
    if md > d
        # add an outer ring
        ring = RingSDF(d/2, (md - d) / 2, abs(s))
        translate3d!(ring, [0, -(l/2 - s), 0])
        front += ring
    elseif md < d
        md = d
        @warn "The lens mechanical diameter is less than the clear optical diameter. Parameter has been ignored."
    end
    back = CylinderSDF(md / 2, (l - s) / 2)

    return (front + back)
end

function render_object!(axis, s::UnionSDF)
    for sdf in s.sdfs
        render_object!(axis, sdf)
    end
end
