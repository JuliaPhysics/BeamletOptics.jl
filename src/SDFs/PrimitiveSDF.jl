"""
    BoxSDF <: AbstractSDF

Implements the box SDF with edge lengths `x`, `y`, and `z`.
Note that these values are stored in the `dimensions` field as:

- `dimensions`::Point3 = (
    len_in_x/2,
    len_in_y/2,
    len_in_z/2,
)
"""
mutable struct BoxSDF{T} <: AbstractSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    dimensions::Point3{T}
end

"""
    BoxSDF(x, y, z)

Creates a [`BoxSDF`](@ref) with:

- `x`: x-dir. edge length in [m]
- `y`: y-dir. edge length in [m]
- `z`: z-dir. edge length in [m]
"""
function BoxSDF(x::X, y::Y, z::Z) where {X<:Real, Y<:Real, Z<:Real}
    T = float(promote_type(X, Y, Z))
    return BoxSDF{T}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        Point3{T}(0),
        Point3{T}(x/2, y/2, z/2)
    )
end

thickness(s::BoxSDF) = 2*s.dimensions[2]

function sdf(box::BoxSDF{T}, point) where T
    p = _world_to_sdf(box, point)
    q = abs.(p) - box.dimensions
    l = norm(max.(q, zero(T))) + min(max(q[1], max(q[2], q[3])), zero(T))
    return l
end

function bounding_sphere(s::BoxSDF{T}) where T
    return (Point3{T}(0), norm(s.dimensions))
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

function CylinderSDF(r::R, h::H) where {R<:Real, H<:Real}
    T = float(promote_type(R, H))
    return CylinderSDF{T}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        Point3{T}(0),
        r,
        h)
end

function sdf(cylinder::CylinderSDF{T}, point) where T
    p = _world_to_sdf(cylinder, point)
    d = abs.(Point2(norm(Point2(p[1], p[3])), p[2])) -
        Point2(cylinder.radius, cylinder.height)
    return min(maximum(d), zero(T)) + norm(max.(d, zero(T)))
end

function bounding_sphere(s::CylinderSDF{T}) where T
    return (Point3{T}(0), sqrt(s.radius^2 + (s.height / 2)^2))
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

function bounding_sphere(s::CutSphereSDF{T}) where T
    return (Point3{T}(0), s.radius)
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

function bounding_sphere(s::RingSDF{T}) where T
    r = sqrt((s.inner_radius + s.hwidth)^2 + s.hthickness^2)
    return (Point3{T}(0), r)
end

"""
    RightAnglePrismSDF <: AbstractSDF

Implements the `SDF` of a right angle prism with symmetric leg length `l` and height `h`.
Note that these values are stored in the `dimensions` field as:

dimensions::Point3 = (
    leg_length,     # dim in x
    leg_length,     # dim in y
    height,         # dim in z
)

!!! info "Alignment"
    Note that the prism is not aligned with the positive y-axis!
"""
mutable struct RightAnglePrismSDF{T} <: AbstractSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    dimensions::Point3{T}
end

"""
    RightAnglePrismSDF(leg_length, height)

Constructs a symmetric right angle prism with `leg_length` in x and y and `height` z in [m].
"""
function RightAnglePrismSDF(leg_length::L, height::H) where {L, H}
    T = promote_type(L, H)
    return RightAnglePrismSDF{T}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        Point3{T}(0),
        Point3{T}(leg_length/2, leg_length/2, height/2))
end

function sdf(prism:: RightAnglePrismSDF{T}, point) where T
    p = _world_to_sdf(prism, point)
    q = abs.(p) - prism.dimensions
    box_dist = norm(max.(q, zero(T))) + min(max(q[1], max(q[2], q[3])), zero(T))
    pln_dist = (p[1] + p[2]) / sqrt(2)
    return max(box_dist, pln_dist)
end

function bounding_sphere(s::RightAnglePrismSDF{T}) where T
    return (Point3{T}(0), norm(s.dimensions))
end

# Analytic normal3d definitions for Primitive SDFs

# BoxSDF
function normal3d(box::BoxSDF{T}, point) where {T}
    p = _world_to_sdf(box, point)
    q = abs.(p) ./ box.dimensions
    max_idx = argmax(q)
    n_local = Point3{T}(
        max_idx == 1 ? sign(p[1]) : 0,
        max_idx == 2 ? sign(p[2]) : 0,
        max_idx == 3 ? sign(p[3]) : 0
    )
    return box.dir * n_local
end

# CylinderSDF
function _cylinder_local_normal(cylinder::CylinderSDF{T}, p::Point3{T}) where {T}
    half_h = cylinder.height / 2
    d_top = abs(p[2] - half_h)
    d_bottom = abs(p[2] + half_h)
    r_xz = norm(Point2(p[1], p[3]))
    d_side = abs(r_xz - cylinder.radius)
    
    if d_top <= d_bottom && d_top <= d_side
        return Point3{T}(0, 1, 0)
    elseif d_bottom <= d_top && d_bottom <= d_side
        return Point3{T}(0, -1, 0)
    else
        r_inv = r_xz > 0 ? inv(r_xz) : zero(T)
        return Point3{T}(p[1] * r_inv, 0, p[3] * r_inv)
    end
end

function normal3d(cylinder::CylinderSDF{T}, point) where {T}
    p = _world_to_sdf(cylinder, point)
    n_local = _cylinder_local_normal(cylinder, p)
    return cylinder.dir * n_local
end

# RightAnglePrismSDF
function _prism_face_index(prism::RightAnglePrismSDF{T}, p::Point3{T}) where {T}
    d_hyp = abs((p[1] + p[2]) / sqrt(T(2)))
    d_leg1 = abs(p[1] + prism.dimensions[1])
    d_leg2 = abs(p[2] + prism.dimensions[2])
    d_side = abs(abs(p[3]) - prism.dimensions[3])

    min_d = min(d_hyp, min(d_leg1, min(d_leg2, d_side)))
    if min_d == d_hyp
        return 1 # hypotenuse
    elseif min_d == d_leg1
        return 2 # leg1
    elseif min_d == d_leg2
        return 3 # leg2
    else
        return 4 # side
    end
end

function normal3d(prism::RightAnglePrismSDF{T}, point) where {T}
    p = _world_to_sdf(prism, point)
    idx = _prism_face_index(prism, p)
    n_local = if idx == 1
        Point3{T}(1/sqrt(T(2)), 1/sqrt(T(2)), 0)
    elseif idx == 2
        Point3{T}(-1, 0, 0)
    elseif idx == 3
        Point3{T}(0, -1, 0)
    else
        Point3{T}(0, 0, sign(p[3]))
    end
    return prism.dir * n_local
end
