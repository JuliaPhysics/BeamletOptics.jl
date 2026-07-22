"""
    OffAxisParaboloidSDF{T} <: AbstractSDF{T}

Signed distance function representation of an Off-Axis Parabolic (OAP) mirror substrate.
The front concave surface represents a paraboloid of revolution with parent focal length `f`
offset by distance `x_off` along the parent x-axis.

# Fields

- `f`: Parent paraboloid focal length [m]
- `x_off`: Off-axis distance from parent paraboloid vertex to aperture center [m]
- `diameter`: Mirror aperture diameter [m]
- `thickness`: Substrate thickness in y-direction [m]
- `pos`: Position point in world coordinates [m]
- `dir`: Rotation matrix in world coordinates
- `transposed_dir`: Transposed rotation matrix
"""
mutable struct OffAxisParaboloidSDF{T} <: AbstractSDF{T}
    f::T
    x_off::T
    diameter::T
    thickness::T
    pos::Point3{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
end

function OffAxisParaboloidSDF(f::F, x_off::X, diameter::D, thickness::Th) where {F<:Real, X<:Real, D<:Real, Th<:Real}
    T = float(promote_type(F, X, D, Th))
    return OffAxisParaboloidSDF{T}(
        T(f), T(x_off), T(diameter), T(thickness),
        Point3{T}(0),
        SMatrix{3, 3, T, 9}(I),
        SMatrix{3, 3, T, 9}(I)
    )
end

function bounding_sphere(s::OffAxisParaboloidSDF{T}) where {T}
    r_max = s.diameter / 2
    sag_max = abs(-(((r_max + s.x_off)^2) / (4 * s.f) - (s.x_off^2) / (4 * s.f)))
    y_center = (s.thickness - sag_max) / 2
    r_bound = sqrt(r_max^2 + ((s.thickness + sag_max) / 2)^2) + T(0.05)
    return Point3{T}(0, y_center, 0), r_bound
end

function sdf(s::OffAxisParaboloidSDF{T}, point) where {T}
    p = _world_to_sdf(s, point)
    x, y, z = p[1], p[2], p[3]

    y_surf = -(((x + s.x_off)^2 + z^2) / (4 * s.f) - (s.x_off^2) / (4 * s.f))
    grad_norm = sqrt(one(T) + ((x + s.x_off) / (2 * s.f))^2 + (z / (2 * s.f))^2)

    d_front = (y_surf - y) / grad_norm
    d_cyl = sqrt(x^2 + z^2) - s.diameter / 2
    d_back = y - s.thickness

    return max(d_front, d_cyl, d_back)
end
