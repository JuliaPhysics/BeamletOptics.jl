abstract type AbstractAsphericalSurfaceSDF{T} <: AbstractLensSDF{T} end

"""
    ConvexAsphericalSurfaceSDF

Constructs an aspheric lens with a convex-like surface according to ISO10110.

!!! note
    Currently, it is assumed that the aspheric surface is convex if the radius is positive.
    There might be unexpected effects for complex shapes which do not show a generally convex
    behavior.

- `coefficients` : (even) coefficients of the asphere, starting with A2 (r²).
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
    max_sag::Point2{T}
    edge::Point2{T}
end

function thickness(s::ConvexAsphericalSurfaceSDF)
    sag = aspheric_equation(s.diameter / 2, 1 / s.radius, s.conic_constant, s.coefficients)
    return s.max_sag[1] > 0 && sag < 0 ? s.max_sag[1] : abs(sag)
end

# Constructor for ConvexAsphericalSurfaceSDF
function ConvexAsphericalSurfaceSDF(
        coefficients::Vector{T}, radius::T, conic_constant::T, diameter::T) where {T}
    return ConvexAsphericalSurfaceSDF{T}(
        coefficients,
        radius,
        conic_constant,
        diameter,
        Point3{T}(0),
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        Point2(max_aspheric_value(1/radius, conic_constant, coefficients, diameter)),
        aspheric_edge(1/radius, conic_constant, coefficients, diameter)
    )
end

function max_aspheric_value(c, k, α_coeffs, d)
    f(r) = aspheric_equation(r, c, k, α_coeffs)
    # Use the first component of the gradient as the derivative with respect to r.
    fprime(r) = gradient_aspheric_equation(r, c, k, α_coeffs)[1]
    # Avoid r=0 (which might be problematic) by starting at a small positive value.
    a = 1e-8
    b = d/2
    # If there is no sign change, fallback to the endpoint
    if sign(fprime(a)) == sign(fprime(b))
        r_max = (abs(f(a)) > abs(f(b)) ? a : b)
    else
        r_max = find_zero_bisection(fprime, a, b)
    end
    return f(r_max), r_max
end

"""
    aspheric_edge(c, k, α_coeffs, d)

Returns the sag of the asphere at its edge `r = d/2` and the norm of the gradient of its 2D
distance field there, as a `Point2`. Both are constant for a surface, hence the SDFs store them
instead of evaluating the asphere at the edge for every distance.
"""
function aspheric_edge(c, k, α_coeffs, d)
    z = aspheric_equation(d / 2, c, k, α_coeffs)
    return Point2(z, norm(gradient_aspheric_equation(d / 2, c, k, α_coeffs)))
end

"""
    ConcaveAsphericalSurfaceSDF

Constructs an aspheric lens with a concave-like surface according to ISO10110.

!!! note
    Currently, it is assumed that the aspheric surface is concave if the radius is negative.
    There might be unexpected effects for complex shapes which do not show a generally concave
    behavior.

- `coefficients` : (even) coefficients of the asphere, starting with A2 (r²).
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
    max_sag::Point2{T}
    edge::Point2{T}
end

function thickness(s::ConcaveAsphericalSurfaceSDF{T}) where {T}
    sag = aspheric_equation(s.diameter / 2, 1 / s.radius, s.conic_constant, s.coefficients)
    return s.max_sag[1] > 0 && sag < 0 ? abs(sag) : zero(T)
end

# Constructor for ConcaveAsphericalSurfaceSDF
function ConcaveAsphericalSurfaceSDF(
        coefficients::V, radius::T, conic_constant::T, diameter::T,
        mechanical_diameter::T = diameter) where {T, V <: AbstractVector{T}}
    return ConcaveAsphericalSurfaceSDF(
        coefficients,
        radius,
        conic_constant,
        diameter,
        mechanical_diameter,
        Point3{T}(0),
        SMatrix{3, 3}(one(T) * I),
        SMatrix{3, 3}(one(T) * I),
        Point2(max_aspheric_value(1/radius, conic_constant, coefficients, diameter)),
        aspheric_edge(1/radius, conic_constant, coefficients, diameter)
    )
end

"""
aspheric_equation(r, c, k, α_coeffs)

The aspheric surface equation. The asphere is defined by:
- `c` : The curvature (1/radius) of the surface
- `k` : The conic constant of the surface
- `α_coeffs` : The (even) aspheric coefficients of the terms `Σᵢ αᵢ r²ⁱ`, starting with A2 (r²).

This function returns NaN if the square root argument becomes negative.

!!! note
Only even aspheres are implemented at the moment. This will change soon.

"""
function aspheric_equation(r, c, k, α_coeffs)
    r2 = r^2
    _aspheric_sqrt_arg(r2, c, k) < 0 && return NaN
    return _aspheric_sag_r2(r2, c, k, α_coeffs)
end

function aspheric_equation(r::Real, a::AbstractAsphericalSurfaceSDF)
    aspheric_equation(r, 1 / a.radius, a.conic_constant, a.coefficients)
end

function gradient_aspheric_equation(r, c, k, α_coeffs)
    r2 = r^2
    _aspheric_sqrt_arg(r2, c, k) < 0 && return NaN
    return Point2(-r * aspheric_slope_r2(r2, c, k, α_coeffs), 1)
end

# The argument of the square root in the sag of the asphere; the asphere is not defined where it is
# negative.
_aspheric_sqrt_arg(r2, c, k) = 1 - (1 + k) * c^2 * r2

# The sag, the slope divided by the radius, G = (dz/dr) / r, and the squared slope (dz/dr)² = r² G²
# of the asphere as functions of r². They are smooth on the axis, also for automatic
# differentiation, in contrast to functions of the radius r.
function _aspheric_sag_r2(r2, c, k, α_coeffs)
    # Even aspheric terms Σᵢ αᵢ r²ⁱ via Horner in r²
    acc = zero(r2 * zero(eltype(α_coeffs)))
    for i in reverse(eachindex(α_coeffs))
        acc = muladd(acc, r2, α_coeffs[i])
    end
    return c * r2 / (1 + sqrt(_aspheric_sqrt_arg(r2, c, k))) + acc * r2
end

function aspheric_slope_r2(r2, c, k, α_coeffs)
    sq = sqrt(_aspheric_sqrt_arg(r2, c, k))
    # Even aspheric terms Σᵢ 2i αᵢ r²⁽ⁱ⁻¹⁾ via Horner in r²
    acc = zero(r2 * zero(eltype(α_coeffs)))
    for i in reverse(eachindex(α_coeffs))
        acc = muladd(acc, r2, 2i * α_coeffs[i])
    end
    return 2c / (sq + 1) + r2 * (1 + k) * c^3 / (sq * (sq + 1)^2) + acc
end

aspheric_slope2_r2(r2, c, k, α_coeffs) = r2 * aspheric_slope_r2(r2, c, k, α_coeffs)^2

# The smaller of the distance to the aspheric surface, `d_surface`, and the distance to the segments
# that close the 2D perimeter, `d_closure`, divided, i.e. scaled down, by the factor `g` (≥ 1) as the
# distance of the surface is. Within the surface threshold, where a closing segment may touch the
# surface, e.g. at the vertex of a concave asphere, the surface is taken if it is not farther than the
# unscaled segment: both give the same normal, but only the surface gives its change, e.g. its
# derivative in automatic differentiation. The value changes by less than the threshold.
function _min_surface_first(d_surface, d_closure, g)
    near = d_surface ≤ Config.get_sdf_surface_threshold()
    near && d_surface ≤ d_closure + Config.get_sdf_raymarch_eps() / 100 && return d_surface
    return min(d_surface, d_closure / g)
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
    convex_aspheric_surface_distance(r, z, c, k, d, α_coeffs, max_sag, edge)

Calculates the 2D distance field for an aspheric surface at radius `r` away from the optical axis
position `z`. The asphere is defined by:
- `c` : The curvature (1/radius) of the surface
- `k` : The conic constant of the surface
- `d` : The diameter of the asphere
- `α_coeffs` : The (even) aspheric coefficients, starting with A2.
- `max_sag` : the extremal sag and its radius, see `max_aspheric_value`
- `edge` : the sag and the gradient norm at the edge, see `aspheric_edge`

Note that this is not just an infinite aspheric surface and also not a surface segment but
a closed 2D perimeter.

It is intended to pair the SDF derived from this distance field with a cylinder SDF to build a real lens.
"""
function convex_aspheric_surface_distance(r, z, c, k, d, α_coeffs, max_sag, edge, r2 = r^2)
    r2_bound = (d / 2)^2

    # Get edge sag and gradient norm
    z_aspheric_boundary, n_grad_z_boundary = edge

    # Handle points outside the asphere
    if _aspheric_sqrt_arg(r2, c, k) < 0 || r2 > r2_bound
        if z < z_aspheric_boundary
            distance_to_boundary = sqrt((r - sign(r) * d / 2)^2 + (z - z_aspheric_boundary)^2)
        elseif z_aspheric_boundary < z < 0
            distance_to_boundary = sqrt((r - sign(r) * d / 2)^2)
        elseif z > 0 && (sign(c) == 1 && z_aspheric_boundary < 0)
            distance_to_boundary = sqrt((r - sign(r) * d / 2)^2 + z^2)
        else
            distance_to_boundary = sqrt((r - sign(r) * d / 2)^2 + (z - z_aspheric_boundary)^2)
        end
        return distance_to_boundary / n_grad_z_boundary
    end

    # Distance to asphere via r²
    z_aspheric_value = _aspheric_sag_r2(r2, c, k, α_coeffs)
    distance_to_aspheric = abs(z - z_aspheric_value) / sqrt(1 + aspheric_slope2_r2(r2, c, k, α_coeffs))
    if sign(c) == 1 && z_aspheric_boundary < 0
        # Close perimeter at max sag
        p = Point2(r, z)

        a1, b1 = Point2(d / 2, z_aspheric_boundary), Point2(d / 2, max_sag[1])
        a2, b2 = Point2(d / 2, max_sag[1]), Point2(-d / 2, max_sag[1])
        a3, b3 = Point2(-d / 2, max_sag[1]), Point2(-d / 2, z_aspheric_boundary)
        # Distance to closing segments
        d_closure = min(sd_line_segment(p, a1, b1), sd_line_segment(p, a2, b2), sd_line_segment(p, a3, b3))

        if z_aspheric_value < z < max_sag[1]
            return -_min_surface_first(distance_to_aspheric, d_closure, n_grad_z_boundary)
        else
            return _min_surface_first(distance_to_aspheric, d_closure, n_grad_z_boundary)
        end
    else
        # Close perimeter with line segment
        a, b = Point2(d / 2, z_aspheric_boundary), Point2(-d / 2, z_aspheric_boundary)
        d_closure = sd_line_segment(Point2(r, z), a, b)
        if sign(c) * z_aspheric_value < sign(c) * z < sign(c) * z_aspheric_boundary
            return -_min_surface_first(distance_to_aspheric, d_closure, n_grad_z_boundary)
        else
            return _min_surface_first(distance_to_aspheric, d_closure, n_grad_z_boundary)
        end
    end
end

function concave_aspheric_surface_distance(r, z, c, k, d, α_coeffs, max_sag, edge, r2 = r^2)
    r2_bound = (d / 2)^2

    # Get edge sag and gradient norm
    z_aspheric_boundary, n_grad_z_boundary = edge
    # Handle points outside the asphere
    if _aspheric_sqrt_arg(r2, c, k) < 0
        if z < 0
            distance_to_boundary = sqrt((r - sign(r) * d / 2)^2 + z^2)
        elseif 0 < z < z_aspheric_boundary
            distance_to_boundary = sqrt((r - sign(r) * d / 2)^2)
        else
            distance_to_boundary = sqrt((r - sign(r) * d / 2)^2 + (z - z_aspheric_boundary)^2)
        end
        return distance_to_boundary / n_grad_z_boundary
    end

    # Distance to asphere via r²
    z_aspheric_value = _aspheric_sag_r2(r2, c, k, α_coeffs)
    distance_to_aspheric = abs(z - z_aspheric_value) / sqrt(1 + aspheric_slope2_r2(r2, c, k, α_coeffs))
    p = Point2(r, z)

    if max_sag[1] > 0 && z_aspheric_boundary < 0
        a, b = Point2(d / 2, z_aspheric_boundary), Point2(-d / 2, z_aspheric_boundary)
        d_closure = sd_line_segment(p, a, b)
        sdl = d_closure / n_grad_z_boundary

        if r2 > r2_bound
            # Outside aperture, use perimeter segments
            return sdl
        else
            # Negate inside asphere
            if z_aspheric_boundary < z < z_aspheric_value
                return -_min_surface_first(distance_to_aspheric, d_closure, n_grad_z_boundary)
            elseif z_aspheric_boundary > 0 && (0.0 < z < z_aspheric_value)
                return -_min_surface_first(distance_to_aspheric, d_closure, n_grad_z_boundary)
            else
                return _min_surface_first(distance_to_aspheric, d_closure, n_grad_z_boundary)
            end
        end
    else
        a1, b1 = Point2(d / 2, z_aspheric_boundary), Point2(d / 2, 0.0)
        a2, b2 = Point2(d / 2, 0.0), Point2(-d / 2, 0.0)
        a3, b3 = Point2(-d / 2, 0.0), Point2(-d / 2, z_aspheric_boundary)
        # Distance to closing segments
        d_closure = min(sd_line_segment(p, a1, b1), sd_line_segment(p, a2, b2), sd_line_segment(p, a3, b3))

        if r2 > r2_bound
            # Outside aperture, use perimeter segments
            return d_closure / n_grad_z_boundary
        else
            # Negate inside asphere
            if z_aspheric_boundary < 0 && (z_aspheric_value < z < 0.0)
                return -_min_surface_first(distance_to_aspheric, d_closure, n_grad_z_boundary)
            elseif z_aspheric_boundary > 0 && (0.0 < z < z_aspheric_value)
                return -_min_surface_first(distance_to_aspheric, d_closure, n_grad_z_boundary)
            else
                return _min_surface_first(distance_to_aspheric, d_closure, n_grad_z_boundary)
            end
        end
    end
end

"""
    bounding_sphere_of(s::AbstractAsphericalSurfaceSDF)

Returns the sphere around the cylinder of the lens diameter that spans the axial extent of the
surface volume, centered on the axis.
"""
function bounding_sphere_of(s::AbstractAsphericalSurfaceSDF{T}) where T
    lo, hi = axial_extent(s)
    return SingleBoundingSphere(s, Point3{T}(0, (lo + hi) / 2, 0), sqrt((s.diameter / 2)^2 + ((hi - lo) / 2)^2))
end

function sdf(surface::ConvexAsphericalSurfaceSDF{T}, point) where {T}
    p_local = _world_to_sdf(surface, point)
    # Revolve 2D distance around y
    x, z = p_local[1], p_local[3]
    return convex_aspheric_surface_distance(
        _radial(x, z),
        p_local[2],
        1 / surface.radius,
        surface.conic_constant,
        surface.diameter,
        surface.coefficients,
        surface.max_sag,
        surface.edge,
        x^2 + z^2)
end

function sdf(surface::ConcaveAsphericalSurfaceSDF{T}, point) where {T}
    p_local = _world_to_sdf(surface, point)
    # Revolve 2D distance around y
    x, z = p_local[1], p_local[3]
    return concave_aspheric_surface_distance(
        _radial(x, z),
        p_local[2],
        1 / surface.radius,
        surface.conic_constant,
        surface.diameter,
        surface.coefficients,
        surface.max_sag,
        surface.edge,
        x^2 + z^2)
end

"""
    PlanoConvexAsphericalLensSDF(r, l, d=1inch)

Constructs a plano-convex aspheric lens SDF with:

- `r` > 0: front radius
- `l`: lens thickness
- `d`: lens diameter
- `k` : The conic constant of the surface
- `α_coeffs` : The (even) aspheric coefficients, starting with A2 (r²).

The spherical surface is constructed flush with the cylinder surface.
"""
function PlanoConvexAsphericalLensSDF(
        r::R, l::L, d::D, k::K, α_coeffs::AbstractVector{A}) where {R, L, D, K, A}
    T = promote_type(R, L, D, K, A)
    s = aspheric_equation(d / 2, 1 / r, k, α_coeffs)
    # Calculate length of cylindrical section
    l < 0 && error("Specified thickness is shorter than the lens sagitta at the edge")
    front = ConvexAsphericalSurfaceSDF(α_coeffs, r, k, d)
    back = CylinderSDF(d / 2, (l - abs(s)) / 2)
    # Shift and rotate cut spheres into position
    translate3d!(front, [0, -sign(r) * (l / 2 + abs(s) / 2), 0])
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
- `α_coeffs` : The (even) aspheric coefficients, starting with A2 (r²).
- `md`: lens mechanical diameter (default: md = d)

The spherical surface is constructed flush with the cylinder surface.
"""
function PlanoConcaveAsphericalLensSDF(r::R, l::L, d::D, k::K, α_coeffs::AbstractVector{A},
        md::MD = d) where {R, L, D, K, A, MD}
    s = aspheric_equation(d / 2, 1 / r, k, α_coeffs)
    # Calculate length of cylindrical section
    l < 0 && error("Specified thickness is shorter than the lens sagitta at the edge")
    front = ConcaveAsphericalSurfaceSDF(α_coeffs, r, k, d)
    # Shift and rotate asphere into position
    translate3d!(front, [0, -(l - s) / 2, 0])
    # add outer planar ring, if required
    if md > d
        # add an outer ring
        ring = RingSDF(d / 2, (md - d) / 2, abs(s))
        translate3d!(ring, [0, -(l / 2 - s), 0])
        front += ring
    elseif md < d
        md = d
        @warn "The lens mechanical diameter is less than the clear optical diameter. Parameter has been ignored."
    end
    back = CylinderSDF(md / 2, (l - s) / 2)

    return (front + back)
end

"""
    EvenAsphericalSurface{T} <: AbstractRotationallySymmetricSurface{T}

A type representing an aspherical optical surface defined by its radius of curvature,
clear (optical) diameter, conic constant, aspheric coefficients and mechanical diameter.
This surface is rotationally symmetric about its optical axis.

# Fields
- `spherical::SphericalSurface{T}`: The base spherical surface portion of the asphere.
- `conic_constant::T`: The conic constant defining the deviation from a spherical shape.
- `coefficients::AbstractVector{T}`: A vector of even aspheric coefficients for higher-order corrections,
  starting with A2 (r²).

"""
struct EvenAsphericalSurface{T} <: AbstractRotationallySymmetricSurface{T}
    spherical::SphericalSurface{T}
    conic_constant::T
    coefficients::Vector{T}
end

"""
    EvenAsphericalSurface(radius, diameter, conic_constant, coefficients::AbstractVector, mechanical_diameter = diameter)

Construct a `EvenAsphericalSurface` given the radius of curvature, the optical diameter, conic
constant, the aspheric coefficients and optionally the mechanical diameter.
This constructor automatically sets the mechanical diameter equal to the optical diameter.

# Arguments
- `radius`: The radius of curvature of the base spherical surface.
- `diameter`: The clear (optical) diameter of the surface.
- `conic_constant`: The conic constant defining the deviation from a spherical shape.
- `coefficients::AbstractVector`: The even aspheric coefficients of the terms `Σᵢ αᵢ r²ⁱ` in SI units, starting
  with A2 (r²), e.g. `[0, A4 * 1e3^3, A6 * 1e3^5]` for mm-based catalog values without an r² term.
- `mechanical_diameter`: The mechanical diameter of the surface; if not provided, it defaults to `diameter`.

"""
function EvenAsphericalSurface(radius::T1, diameter::T2, conic_constant::T3, coefficients::AbstractVector{T4}, mechanical_diameter::T5=diameter) where {T1, T2, T3, T4, T5}
    T = promote_type(T1, T2, T3, T4, T5)
    st = SphericalSurface{T}(radius, diameter, mechanical_diameter)

    return EvenAsphericalSurface{T}(
        st,
        conic_constant,
        coefficients
    )
end
radius(s::EvenAsphericalSurface) = radius(s.spherical)
diameter(s::EvenAsphericalSurface) = diameter(s.spherical)
mechanical_diameter(s::EvenAsphericalSurface) = mechanical_diameter(s.spherical)

function edge_sag(s::EvenAsphericalSurface{T}, ::AbstractAsphericalSurfaceSDF) where T
    sag = aspheric_equation(
        diameter(s) / 2,
        1 / radius(s),
        s.conic_constant,
        s.coefficients
    )

    return sag
end

function sdf(s::EvenAsphericalSurface, ot::AbstractOrientationType)
    isinf(radius(s)) && return nothing

    return _sdf(s, ot)
end

function _sdf(s::EvenAsphericalSurface, ::ForwardOrientation)
    front = if radius(s) > 0
        ConvexAsphericalSurfaceSDF(s.coefficients, radius(s), s.conic_constant, diameter(s))
    else
        ConcaveAsphericalSurfaceSDF(s.coefficients, radius(s), s.conic_constant, diameter(s))
    end

    return front
end

function _sdf(s::EvenAsphericalSurface, ::BackwardOrientation)
    back = if radius(s) > 0
        ConcaveAsphericalSurfaceSDF(s.coefficients, radius(s), s.conic_constant, diameter(s))
    else
        ConvexAsphericalSurfaceSDF(s.coefficients, radius(s), s.conic_constant, diameter(s))
    end

    return back
end

# The remaining orientation types are only requested for the MeniscusLensSDF
function _sdf(::EvenAsphericalSurface, ::AbstractOrientationType)
    throw(ArgumentError("A meniscus lens, whose center thickness does not exceed the sagitta of its convex surface, \
        must consist of spherical surfaces. Aspherical surfaces are not supported in this case."))
end