abstract type AbstractRotationallySymmetricLensSDF{T} <: AbstractSDF{T} end

"""
    sag(r::Real, l::Real)

Calculates the sag of a cut circle with radius `r` and chord length `l`
"""
sag(r::Real, l::Real) = r - sqrt(r^2 - 0.25 * l^2)

function check_sag(r, d)
    if abs(2r) < d 
        throw(ArgumentError("Radius of curvature (r = $(r)) must be ≥ than half the diameter (d = $(d)) or an illegal shape results!"))
    else
        return nothing
    end
end

"""
    AbstractSphericalLensSDF{T} <: AbstractSDF{T}

An abstract type for SDFs which represent spherical lenses, i.e. biconvex or plano-concave.
"""
abstract type AbstractSphericalLensSDF{T} <: AbstractRotationallySymmetricLensSDF{T} end
abstract type AbstractSphericalSurfaceSDF{T} <: AbstractRotationallySymmetricLensSDF{T} end

mutable struct ConcaveSphericalSurfaceSDF{T} <: AbstractSphericalSurfaceSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    radius::T
    diameter::T
    sag::T
end

function ConcaveSphericalSurfaceSDF(r::R, d::D) where {R, D}
    T = promote_type(R, D)
    check_sag(r, d)
    s = sag(r, d)
    return ConcaveSphericalSurfaceSDF{T}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        Point3{T}(0),
        r, d, s
    )
end

function sdf(css::ConcaveSphericalSurfaceSDF, point)
    p = _world_to_sdf(css, point)
    # cylinder sdf
    ps = p + Point3(0, css.sag/2, 0)
    d = abs.(Point2(norm(Point2(ps[1], ps[3])), ps[2])) -
        Point2(css.diameter/2, css.sag/2)
    sdf1 = min(maximum(d), 0) + norm(max.(d, 0))
    # sphere sdf
    ps = p + Point3(0, css.radius, 0)
    sdf2 = norm(ps) - css.radius
    return max(sdf1, -sdf2)
end

mutable struct ConvexSphericalSurfaceSDF{T} <: AbstractSphericalSurfaceSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    radius::T
    height::T
    w::T
end

function ConvexSphericalSurfaceSDF(radius::R, diameter::D) where {R, D}
    T = promote_type(R, D)
    check_sag(radius, diameter)
    _sag = sag(radius, diameter)
    cutoff_height = radius - _sag 
    w = sqrt(radius^2 - cutoff_height^2)
    return ConvexSphericalSurfaceSDF{T}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        zeros(T, 3),
        radius,
        cutoff_height,
        w)
end

function sdf(css::ConvexSphericalSurfaceSDF, point)
    p = _world_to_sdf(css, point)
    # -p[2] to align surface with neg. y-axis
    q = Point2(norm(Point2(p[1], p[3])), -p[2] + css.height)
    s = max((css.height - css.radius) * q[1]^2 + css.w^2 * (css.height + css.radius - 2 * q[2]),
        css.height * q[1] - css.w * q[2])
    if s < 0
        return norm(q) - css.radius
    elseif q[1] < css.w
        return css.height - q[2]
    else
        return norm(q - Point2(css.w, css.height))
    end
end

"""
    ThinLensSDF(r1, r2, d=1inch)

Constructs a bi-convex thin lens SDF-based shape with:

- `r1 > 0`: radius of convex front
- `r2 > 0`: radius of convex back
- `d`: lens diameter, default value is one inch

The spherical surfaces are constructed flush.
"""
function ThinLensSDF(r1::L, r2::M, d::O = 1inch) where {L, M, O}
    check_sag(r1, d)
    check_sag(r2, d)
    # Create lens halves
    front = ConvexSphericalSurfaceSDF(r1, d)
    back = ConvexSphericalSurfaceSDF(r2, d)
    # Rotate and move back cut sphere
    zrotate3d!(back, π)
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
    check_sag(r1, d)
    check_sag(r2, d)
    s1 = sag(r1, d)
    s2 = sag(r2, d)
    # Calculate length of cylindrical section
    l = l - (s1 + s2)
    front = ConvexSphericalSurfaceSDF(r1, d)
    back = ConvexSphericalSurfaceSDF(r2, d)
    mid = CylinderSDF(d / 2, l / 2)
    # Shift and rotate elements into position
    translate3d!(front, [0, -l / 2, 0])
    zrotate3d!(back, π)
    translate3d!(back, [0, l / 2, 0])
    return (front + mid + back)
end

function BiConcaveLensSDF(r1::L, r2::M, l::N, d::O = 1inch) where {L, M, N, O}
    # create segments
    front = ConcaveSphericalSurfaceSDF(r1, d)
    back = ConcaveSphericalSurfaceSDF(r2, d)
    mid = CylinderSDF(d / 2, l / 2)
    # Shift and rotate subtraction spheres into position
    translate3d!(front, [0, -l / 2, 0])
    zrotate3d!(back, π)
    translate3d!(back, [0, l / 2, 0])
    return (front + mid + back)
end

"""
    BiConcaveLensSDF(r1, r2, l, d=1inch)

Constructs a bi-concave lens SDF with:

- `r1` > 0: radius of concave front
- `r2` > 0: radius of convex back
- `l`: lens thickness
- `d`: lens diameter, default value is one inch
- `md`: mechanical lens diameter, adds an outer ring section to the lens, if `md` > `d`.

The spherical surfaces are constructed flush with the cylinder surface.
"""
function BiConcaveLensSDF(r1::L, r2::M, l::N, d::O, md::MD) where {L, M, N, O, MD}
    if md ≤ d
        throw(ArgumentError("Mech. diameter must be larger than lens diameter!"))
    end
    # generate ring-less shape
    shape = BiConcaveLensSDF(r1, r2, l, d)
    # add an outer ring
    l += shape.sdfs[1].sag + shape.sdfs[3].sag
    ring = RingSDF(d/2, (md - d) / 2, l)
    shape += ring
    return shape
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
    _sag = sag(r, d)
    # Calculate length of cylindrical section
    l = l - _sag
    front = ConvexSphericalSurfaceSDF(r, d)
    back = CylinderSDF(d / 2, l / 2)
    # Shift and rotate cut spheres into position
    translate3d!(front, [0, -l / 2, 0])
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