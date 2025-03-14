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
    xrotate3d!(s, π / 2)
    translate3d!(s, [0, radius, 0])

    return s
end

thickness(s::CutCylinderSDF) = abs(sag(s.radius, s.diameter))

function sdf(s::CutCylinderSDF{T}, point) where {T}
    p = _world_to_sdf(s, point)

    return op_extrude_x(
        p,
        _p -> sdf_cut_disk(
            _p,
            s.radius,
            √(s.radius^2 - (s.diameter / 2)^2)
        ),
        s.height / 2
    )
end

function sdf_cut_disk(point::Point2, r, h)
    w = √(r^2 - h^2)
    p = Point2(abs(point[1]), point[2])

    s = max((h - r) * p[1]^2 + w^2 * (h + r - 2 * p[2]), h * p[1] - w * p[2])

    return (s < 0) ? norm(p) - r :
           (p[1] < w) ? h - p[2] :
           norm(p - Point2(w, h))
end


"""
    ConcaveCylinderSDF <: AbstractSDF

Implements the `SDF` of a concave cut cylinder with radius `r`, diameter `d` and height `h`.
"""
mutable struct ConcaveCylinderSDF{T} <: AbstractSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    radius::T
    diameter::T
    height::T
end

"""
    ConcaveCylinderSDF(radius, diameter, height)

Constructs a concave cut cylinder with radius `r`, diameter `d` and height `h` in [m].
"""
function ConcaveCylinderSDF(radius::R, diameter::D, height::H) where {R, D, H}
    T = promote_type(R, D, H)
    s = ConcaveCylinderSDF{T}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        Point3{T}(0),
        radius,
        diameter,
        height
    )

    return s
end

thickness(::ConcaveCylinderSDF{T}) where T = zero(T)

function sdf(s::ConcaveCylinderSDF{T}, point) where {T}
    p = _world_to_sdf(s, point)

    # cylinder
    _sag = sag(abs(s.radius), s.diameter)
    ps = p + Point3(zero(T), -s.radius, zero(T))
    d = abs.(Point2(norm(Point2(p[3], ps[2])), ps[1])) -
        Point2(abs(s.radius), s.height/2)
    c = min(maximum(d), zero(T)) + norm(max.(d, zero(T)))

    # box
    pp = p + Point3(0, -_sag/2*sign(s.radius), 0)
    q = abs.(pp) - Point3(s.height/2, _sag/2, s.diameter/2)
    l = norm(max.(q, zero(T))) + min(max(q[1], max(q[2], q[3])), zero(T))

    return max(l, -c)
end

# Surface API implementation

abstract type AbstractCylindricSurface{T} <: AbstractSurface{T} end

mechanical_diameter(s::AbstractCylindricSurface) = s.mechanical_diameter
radius(s::AbstractCylindricSurface) = s.radius
diameter(s::AbstractCylindricSurface) = s.diameter
height(s::AbstractCylindricSurface) = s.height

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
CylindricSurface(radius::T, diameter::T, height::T) where {T} = CylindricSurface{T}(
    radius, diameter, height, diameter)

edge_sag(::CylindricSurface, sd::CutCylinderSDF) = thickness(sd)
edge_sag(::CylindricSurface, sd::ConcaveCylinderSDF) = thickness(sd.cut_cylinder_sdf)

function sdf(s::CylindricSurface, ot::AbstractOrientationType)
    isinf(radius(s)) && return nothing

    return _sdf(s, ot)
end

function _sdf(s::CylindricSurface, ::ForwardOrientation)
    front = if radius(s) > 0
        CutCylinderSDF(radius(s), diameter(s), height(s))
    else
        ConcaveCylinderSDF(radius(s), diameter(s), height(s))
    end

    return front
end

function _sdf(s::CylindricSurface, ::BackwardOrientation)
    back = if radius(s) > 0
        ConcaveCylinderSDF(radius(s), diameter(s), height(s))
    else
        CutCylinderSDF(-radius(s), diameter(s), s.height)
    end

    return back
end
