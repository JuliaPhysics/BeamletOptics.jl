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
    T = promote_type(X, Y, Z)
    return BoxSDF{T}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        Point3{T}(0),
        Point3{T}(x/2, y/2, z/2)
    )
end

function sdf(box::BoxSDF{T}, point) where T
    p = _world_to_sdf(box, point)
    q = abs.(p) - box.dimensions
    l = norm(max.(q, zero(T))) + min(max(q[1], max(q[2], q[3])), zero(T))
    return l
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

function sdf(cylinder::CylinderSDF{T}, point) where T
    p = _world_to_sdf(cylinder, point)
    d = abs.(Point2(norm(Point2(p[1], p[3])), p[2])) -
        Point2(cylinder.radius, cylinder.height)
    return min(maximum(d), zero(T)) + norm(max.(d, zero(T)))
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
