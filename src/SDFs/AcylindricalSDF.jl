abstract type AbstractAcylindricSurfaceSDF{T} <: AbstractCylindricalSurfaceSDF{T} end

function aspheric_equation(r::Real, a::AbstractAcylindricSurfaceSDF)
    aspheric_equation(r, 1 / a.radius, a.conic_constant, a.coefficients)
end

"""
    AconvexCylinderSDF{T} <: AbstractAcylindricSurfaceSDF{T}

Implements the `SDF` of a cut cylinder with radius `r`, diameter `d` and height `h`.
"""
mutable struct AconvexCylinderSDF{T} <: AbstractAcylindricSurfaceSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    radius::T
    diameter::T
    height::T
    conic_constant::T
    coefficients::Vector{T}
    max_sag::Point2{T}
end

"""
    AconvexCylinderSDF(radius, diameter, height)

Constructs an aconvex cylinder with radius `r`, diameter `d` and height `h` in [m].
The acylindric shape is defined by its `conic_constant` and the `coefficients` for the
even aspheric polynomoial.
"""
function AconvexCylinderSDF(radius::R, diameter::D, height::H, conic_constant::CC, coefficients::AbstractVector{TT}) where {R, D, H, CC, TT}
    T = promote_type(R, D, H, CC, TT)
    s = AconvexCylinderSDF{T}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        Point3{T}(0),
        radius,
        diameter,
        height,
        conic_constant,
        coefficients,
        Point2(max_aspheric_value(1/radius, conic_constant, coefficients, diameter))
    )

    return s
end

function thickness(s::AconvexCylinderSDF)
    abs(aspheric_equation(s.diameter / 2, 1 / s.radius, s.conic_constant, s.coefficients))
end

function sdf(s::AconvexCylinderSDF{T}, point) where {T}
    p = _world_to_sdf(s, point)

    _pp = Point3{T}(p[1], p[3], p[2]) # xzy
    # rotate 2D sdf around the optical axis
    sdf_v = op_extrude_x(_pp,
        x -> convex_aspheric_surface_distance(
            x[1],
            x[2],
            1 / s.radius,
            s.conic_constant,
            s.diameter,
            s.coefficients,
            s.max_sag
        ),
        s.height/2
        )

    return sdf_v
end

# """
#     AconcaveCylinderSDF{T} <: AbstractAcylindricSurfaceSDF{T}

# Implements the `SDF` of a concave cylinder with radius `r`, diameter `d` and height `h`.
# """
# mutable struct AconcaveCylinderSDF{T} <: AbstractAcylindricSurfaceSDF{T}
#     dir::SMatrix{3, 3, T, 9}
#     transposed_dir::SMatrix{3, 3, T, 9}
#     pos::Point3{T}
#     radius::T
#     diameter::T
#     height::T
#     coefficients::Vector{T}
#     conic_constant::T
#     max_sag::Point2{T}
# end

# """
#     AconcaveCylinderSDF(radius, diameter, height)

# Constructs a aconcave cylinder with radius `r`, diameter `d` and height `h` in [m].
# The acylindric shape is defined by its `conic_constant` and the `coefficients` for the
# even aspheric polynomoial.
# """
# function AconcaveCylinderSDF(radius::R, diameter::D, height::H) where {R, D, H}
#     T = promote_type(R, D, H)
#     s = AconcaveCylinderSDF{T}(
#         Matrix{T}(I, 3, 3),
#         Matrix{T}(I, 3, 3),
#         Point3{T}(0),
#         radius,
#         diameter,
#         height
#     )

#     return s
# end

# thickness(::AconcaveCylinderSDF{T}) where T = zero(T)

# function sdf(s::AconcaveCylinderSDF{T}, point) where {T}
#     p = _world_to_sdf(s, point)

#     # cylinder
#     _sag = sag(abs(s.radius), s.diameter)
#     ps = p + Point3(zero(T), -s.radius, zero(T))
#     d = abs.(Point2(norm(Point2(p[3], ps[2])), ps[1])) -
#         Point2(abs(s.radius), s.height/2)
#     c = min(maximum(d), zero(T)) + norm(max.(d, zero(T)))

#     # box
#     pp = p + Point3(0, -_sag/2*sign(s.radius), 0)
#     q = abs.(pp) - Point3(s.height/2, _sag/2, s.diameter/2)
#     l = norm(max.(q, zero(T))) + min(max(q[1], max(q[2], q[3])), zero(T))

#     return max(l, -c)
# end

# Surface API implementation

abstract type AbstractAcylindricSurface{T} <: AbstractCylindricSurface{T} end

mechanical_diameter(s::AbstractCylindricSurface) = s.mechanical_diameter
radius(s::AbstractCylindricSurface) = s.radius
diameter(s::AbstractCylindricSurface) = s.diameter
height(s::AbstractCylindricSurface) = s.height
conic_constant(s::AbstractAcylindricSurface) = s.conic_constant
coefficients(s::AbstractAcylindricSurface) = s.coefficients

"""
    AcylindricSurface{T} <: AbstractAcylindricSurface{T}

A type representing an acylindric optical surface defined by its radius of curvature, diameter,
height, mechanical diameter, conic constant and even aspheric coefficients.
It is therefore a cylindric surface with a deviation from the perfect cylindric shape.

# Fields
- `radius::T`: The radius of curvature of the curved surface.
- `diameter::T`: The clear (optical) aperture of the surface.
- `height::T` : The height/length of the uncurved surface direction
- `conic_constant::T` : The conic_constant of the curved surface
- `coefficients::Vector{T}` : The coefficients of the even aspherical equation for the curved surface.
- `mechanical_diameter::T`: The overall mechanical diameter of the surface. In many cases, this is equal
  to the optical diameter, but it can be set independently if the mechanical mount requires a larger dimension.

"""
struct AcylindricSurface{T} <: AbstractAcylindricSurface{T}
    radius::T
    diameter::T
    height::T
    conic_constant::T
    coefficients::Vector{T}
    mechanical_diameter::T
end

"""
    AcylindricSurface(radius, diameter, height, conic_constant, coefficients)

Construct a `CylindricSurface` given the radius of curvature, optical diameter and height.
This constructor automatically sets the mechanical diameter equal to the optical diameter.

# Arguments
- `radius`: The radius of curvature of the curved surface.
- `diameter`: The clear (optical) diameter of the surface.
- `height`: The height/length of the uncurved surface direction.
- `conic_constant::T` : The conic_constant of the curved surface
- `coefficients::Vector{T}` : The coefficients of the even aspherical equation for the curved surface.
"""
function AcylindricSurface(radius::T1, diameter::T2, height::T3, conic_constant::T4, coefficients::AbstractVector{T5}) where {T1, T2, T3, T4, T5}
    T = promote_type(T1, T2, T3, T4, T5)

    return AcylindricSurface{T}(
    radius, diameter, height, conic_constant, coefficients, diameter)
end

#edge_sag(::CylindricSurface, sd::AconvexCylinderSDF) = thickness(sd)
#edge_sag(::CylindricSurface, sd::AconcaveCylinderSDF) = thickness(sd.cut_cylinder_sdf)

function sdf(s::AcylindricSurface, ot::AbstractOrientationType)
    isinf(radius(s)) && return nothing

    return _sdf(s, ot)
end

function _sdf(s::AcylindricSurface, ::ForwardOrientation)
    front = if radius(s) > 0
        AconvexCylinderSDF(radius(s), diameter(s), height(s), conic_constant(s), coefficients(s))
    else
        AconcaveCylinderSDF(radius(s), diameter(s), height(s), conic_constant(s), coefficients(s))
    end

    return front
end

function _sdf(s::AcylindricSurface, ::BackwardOrientation)
    back = if radius(s) > 0
        AconcaveCylinderSDF(radius(s), diameter(s), height(s), conic_constant(s), coefficients(s))
    else
        AconvexCylinderSDF(-radius(s), diameter(s), height(s), conic_constant(s), coefficients(s))
    end

    return back
end

# custom render code
function render_object!(axis, acyl::AbstractAcylindricSurfaceSDF; color = :red)
    r = diameter(acyl) / 2
    w_vals = LinRange(-r, r, 100)           # aperture coordinate
    z_vals = LinRange(-height(acyl)/2, height(acyl)/2, 20)  # extrusion coordinate

    # Parameterize lateral surface as (w, aspheric_equation(w), z)
    # but swap the w and z directions to fix the 90° rotation.
    # That is, let X_local come from z_vals and Z_local from w_vals.
    X_local = [z for w in w_vals, z in z_vals]
    Y_local = [aspheric_equation(w, acyl) for w in w_vals, z in z_vals]
    Z_local = [w for w in w_vals, z in z_vals]

    # Transform to global coordinates:
    R = orientation(acyl)
    P = position(acyl)
    Xt = R[1,1] .* X_local .+ R[1,2] .* Y_local .+ R[1,3] .* Z_local .+ P[1]
    Yt = R[2,1] .* X_local .+ R[2,2] .* Y_local .+ R[2,3] .* Z_local .+ P[2]
    Zt = R[3,1] .* X_local .+ R[3,2] .* Y_local .+ R[3,3] .* Z_local .+ P[3]

    render_surface!(axis, Xt, Yt, Zt; transparency = true, colormap = [color, color])
    render_acylindric_caps!(axis, acyl; color)
    return nothing
end

function render_acylindric_cap!(axis, acyl::AbstractAcylindricSurfaceSDF, top::Bool; color=:red)
    r = diameter(acyl) / 2
    h = height(acyl)
    Xval = top ? (h/2) : -(h/2)   # local X coordinate for top or bottom cap
    Nw   = 60                     # resolution along w
    Ny   = 60                     # resolution along y

    w_vals = LinRange(-r, r, Nw)
    X_local = fill(Xval, Nw, Ny)
    Y_local = similar(X_local)
    Z_local = similar(X_local)
    # aspheric value at the edge
    awe = aspheric_equation(r, acyl)
    for i in 1:Nw
        w = w_vals[i]
        aw = aspheric_equation(w, acyl)
        # If aw>0, sweep from awe to aw; if aw<0, sweep from aw to awe
        y_range = aw >= 0 ? LinRange(awe, aw, Ny) : LinRange(aw, awe, Ny)
        for j in 1:Ny
            X_local[i,j] = Xval
            Y_local[i,j] = y_range[j]
            Z_local[i,j] = w
        end
    end

    # Transform local -> global
    R = orientation(acyl)
    P = position(acyl)
    Xt = R[1,1].*X_local .+ R[1,2].*Y_local .+ R[1,3].*Z_local .+ P[1]
    Yt = R[2,1].*X_local .+ R[2,2].*Y_local .+ R[2,3].*Z_local .+ P[2]
    Zt = R[3,1].*X_local .+ R[3,2].*Y_local .+ R[3,3].*Z_local .+ P[3]

    render_surface!(axis, Xt, Yt, Zt; transparency=true, colormap=[color, color])
    return nothing
end

function render_acylindric_caps!(axis, acyl::AbstractAcylindricSurfaceSDF; color=:red)
    render_acylindric_cap!(axis, acyl, true;  color=color)  # top
    render_acylindric_cap!(axis, acyl, false; color=color)  # bottom
    return nothing
end
