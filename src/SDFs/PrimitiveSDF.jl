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

"""
    CutCylinderSDF <: AbstractSDF

Implements the `SDF` of a cut cylinder with radius `r`, diameter `d` and height `h`.
"""
mutable struct CutCylinderSDF{T} <: AbstractSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    radius::T
    diameter::T
    height::T
end

"""
    CutCylinderSDF(radius, diameter, height)

Constructs a cut cylinder with radius `r`, diameter `d` and height `h` in [m].
"""
function CutCylinderSDF(radius::R, diameter::D, height::H) where {R, D, H}
    T = promote_type(R, D, H)
    s = CutCylinderSDF{T}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        Point3{T}(0),
        radius,
        diameter,
        height
    )
    # rotate/shift in position
    xrotate3d!(s, -π / 2)
    translate3d!(s, [0, -(radius - √(radius^2 - (diameter/2)^2)), 0])

    return s
end

function sdf(s::CutCylinderSDF{T}, point) where T
    p = _world_to_sdf(s, point)

    return op_extrude_x(
        p,
        _p -> sdf_cut_disk(
            _p,
            s.radius,
            √(s.radius^2 - (s.diameter/2)^2)
        ),
        s.height / 2
    )
end

function sdf_cut_disk(point::Point2, r, h)
    w = √(r^2-h^2)
    p = Point2(abs(point[1]), point[2])

    s = max((h-r)*p[1]^2+w^2*(h+r-2*p[2]), h*p[1] -w*p[2])

    return (s < 0) ?        norm(p) - r  :
            (p[1] < w) ?    h - p[2] :
                            norm(p-Point2(w, h))

end

abstract type AbstractCylindricSurface{T} <: AbstractSurface{T} end

"""
    CylindricSurface{T} <: AbstractCylindricSurface{T}

A type representing a cylindric optical surface defined by its radius of curvature, diameter,
height and mechanical diameter

# Fields
- `radius::T`: The radius of curvature of the curved surface.
- `diameter::T`: The clear (optical) aperture of the surface.
- `height::T` : The height/length of the uncurved surface direction
- `mechanical_diameter::T`: The overall mechanical diameter of the surface. In many cases, this is equal
  to the optical diameter, but it can be set independently if the mechanical mount requires a larger dimension.

"""
struct CylindricSurface{T} <: AbstractCylindricSurface{T}
    radius::T
    diameter::T
    height::T
    mechanical_diameter::T
end

"""
    CylindricSurface(radius::T, diameter::T, height::T) where T

Construct a `CylindricSurface` given the radius of curvature, optical diameter and height.
This constructor automatically sets the mechanical diameter equal to the optical diameter.

# Arguments
- `radius::T`: The radius of curvature of the curved surface.
- `diameter::T`: The clear (optical) diameter of the surface.
- `height::T`: The height/length of the uncurved surface direction.
"""
CylindricSurface(radius::T, diameter::T, height::T) where T = CylindricSurface{T}(radius, diameter, height, diameter)

mechanical_diameter(s::CylindricSurface) = s.mechanical_diameter

edge_sag(::CylindricSurface, sd::CutCylinderSDF) = abs(sd.radius - √(sd.radius^2 - (sd.diameter/2)^2))

function sdf(s::CylindricSurface, ::ForwardOrientation)
    front = if radius(s) > 0
        CutCylinderSDF(radius(s), diameter(s), s.height)
    else
        throw(ArgumentError("CylindricSurface with negative radius as front surface not supported"))
    end

    return front
end

function sdf(s::CylindricSurface, ::BackwardOrientation)
    back = if radius(s) > 0
        throw(ArgumentError("CylindricSurface with positive radius as back surface not supported"))
    else
        CutCylinderSDF(radius(s), diameter(s), s.height)
    end

    return front
end
