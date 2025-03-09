abstract type AbstractAsphericalSurfaceSDF{T} <: AbstractLensSDF{T} end

"""
    ConvexAsphericalSurfaceSDF

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
    max_sag::Point2{T}
end

function thickness(s::ConvexAsphericalSurfaceSDF)
    abs(aspheric_equation(s.diameter / 2, 1 / s.radius, s.conic_constant, s.coefficients))
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
        Point2(max_aspheric_value(1/radius, conic_constant, coefficients, diameter))
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
    ConcaveAsphericalSurfaceSDF

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
    max_sag::Point2{T}
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
        Point2(max_aspheric_value(1/radius, conic_constant, coefficients, diameter))
    )
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
    sum_α = sum(i -> α_coeffs[i] * r2^(i), eachindex(α_coeffs))
    return c * r2 / (1 + sqrt(sqrt_arg)) + sum_α
end

function aspheric_equation(r::Real, a::AbstractAsphericalSurfaceSDF)
    aspheric_equation(r, 1 / a.radius, a.conic_constant, a.coefficients)
end

function gradient_aspheric_equation(r, c, k, α_coeffs)
    Ri = 1 / c
    sqrt_arg = 1 - r^2 * (1 + k) / Ri^2
    sqrt_arg < 0 && return NaN
    gr = 2 * r / (Ri * (√(sqrt_arg) + 1)) +
         r^3 * (1 + k) / (Ri^3 * √(sqrt_arg) * (√(sqrt_arg) + 1)^2)
    sum_r = sum(m -> 2 * (m) * α_coeffs[m] * r^(2(m - 1) + 1), eachindex(α_coeffs))

    return Point2(-sum_r - gr, 1)
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
function convex_aspheric_surface_distance(r, z, c, k, d, α_coeffs, max_sag)
    r2 = r^2
    r2_bound = (d / 2)^2
    z_aspheric_value = aspheric_equation(r, c, k, α_coeffs)
    grad_z = gradient_aspheric_equation(r, c, k, α_coeffs)

    z_aspheric_boundary = aspheric_equation(d / 2, c, k, α_coeffs)
    grad_z_boundary = gradient_aspheric_equation(d / 2, c, k, α_coeffs)

    # Handling NaN for points outside the aspheric surface
    if isnan(z_aspheric_value) || isnan(grad_z) || r2 > r2_bound
        if z < z_aspheric_boundary
            distance_to_boundary = sqrt((r - sign(r) * d / 2)^2 + (z - z_aspheric_boundary)^2)
        elseif z_aspheric_boundary < z < 0
            distance_to_boundary = sqrt((r - sign(r) * d / 2)^2)
        elseif z > 0 && (sign(c) == 1 && z_aspheric_boundary < 0)
            distance_to_boundary = sqrt((r - sign(r) * d / 2)^2 + z^2)
        else
            distance_to_boundary = sqrt((r - sign(r) * d / 2)^2 + (z - z_aspheric_boundary)^2)
        end
        return distance_to_boundary / norm(grad_z_boundary)
    end

    # Calculate distance to the aspheric surface
    distance_to_aspheric = abs(z - z_aspheric_value) / norm(grad_z)
    if sign(c) == 1 && z_aspheric_boundary < 0
        # the asphere curves towards negative sag so we need to close the perimeter
        # at its max sag value instead of the boundary
        p = Point2(r, z)
        n_gzb = norm(grad_z_boundary)

        a1, b1 = Point2(d / 2, z_aspheric_boundary), Point2(d / 2, max_sag[1])
        sdl1 = sd_line_segment(p, a1, b1) / n_gzb
        a2, b2 = Point2(d / 2, max_sag[1]), Point2(-d / 2, max_sag[1])
        sdl2 = sd_line_segment(p, a2, b2) / n_gzb
        a3, b3 = Point2(-d / 2, max_sag[1]), Point2(-d / 2, z_aspheric_boundary)
        sdl3 = sd_line_segment(p, a3, b3) / n_gzb

        if z_aspheric_value < z < max_sag[1]
            return -min(distance_to_aspheric, sdl1, sdl2, sdl3)
        else
            return min(distance_to_aspheric, sdl1, sdl2, sdl3)
        end
    else
        # If we are inside the aperture, let's close the perimeter with a line segment
        a, b = Point2(d / 2, z_aspheric_boundary), Point2(-d / 2, z_aspheric_boundary)
        sdl = sd_line_segment(Point2(r, z), a, b) / norm(grad_z_boundary)
        if sign(c) * z_aspheric_value < sign(c) * z < sign(c) * z_aspheric_boundary
            return -min(sdl, distance_to_aspheric)
        else
            return min(sdl, distance_to_aspheric)
        end
    end
end

function concave_aspheric_surface_distance(r, z, c, k, d, α_coeffs, max_sag)
    r2 = r^2
    r2_bound = (d / 2)^2
    z_aspheric_value = aspheric_equation(r, c, k, α_coeffs)
    grad_z = gradient_aspheric_equation(r, c, k, α_coeffs)

    z_aspheric_boundary = aspheric_equation(d / 2, c, k, α_coeffs)
    grad_z_boundary = gradient_aspheric_equation(d / 2, c, k, α_coeffs)
    # Handling NaN for points outside the aspheric surface
    if isnan(z_aspheric_value) || isnan(grad_z)
        if z < 0
            distance_to_boundary = sqrt((r - sign(r) * d / 2)^2 + z^2)
        elseif 0 < z < z_aspheric_boundary
            distance_to_boundary = sqrt((r - sign(r) * d / 2)^2)
        else
            distance_to_boundary = sqrt((r - sign(r) * d / 2)^2 + (z - z_aspheric_boundary)^2)
        end
        return distance_to_boundary / norm(grad_z_boundary)
        # distance_to_boundary = sqrt((r - sign(r) * d / 2)^2 + z^2)
        # return distance_to_boundary / norm(grad_z_boundary)
    end

    distance_to_aspheric = abs(z - z_aspheric_value) / norm(grad_z)
    p = Point2(r, z)

    if max_sag[1] > 0 && z_aspheric_boundary < 0
        a, b = Point2(d / 2, z_aspheric_boundary), Point2(-d / 2, z_aspheric_boundary)
        sdl = sd_line_segment(p, a, b) / norm(grad_z_boundary)

        if r2 > r2_bound
            # we are outside of the aperture, so the shortest distance can be determined
            # by the perimeter line segments which include the boundary points of the asphere
            return sdl
        else
            # return a negative sign within the asphere
            if z_aspheric_boundary < z < z_aspheric_value
                return -min(distance_to_aspheric, sdl)
            elseif z_aspheric_boundary > 0 && (0.0 < z < z_aspheric_value)
                return -min(distance_to_aspheric, sdl)
            else
                return min(distance_to_aspheric, sdl)
            end
        end
    else
        a1, b1 = Point2(d / 2, z_aspheric_boundary), Point2(d / 2, 0.0)
        sdl1 = sd_line_segment(p, a1, b1) / norm(grad_z_boundary)
        a2, b2 = Point2(d / 2, 0.0), Point2(-d / 2, 0.0)
        sdl2 = sd_line_segment(p, a2, b2) / norm(grad_z_boundary)
        a3, b3 = Point2(-d / 2, 0.0), Point2(-d / 2, z_aspheric_boundary)
        sdl3 = sd_line_segment(p, a3, b3) / norm(grad_z_boundary)

        if r2 > r2_bound
            # we are outside of the aperture, so the shortest distance can be determined
            # by the perimeter line segments which include the boundary points of the asphere
            return min(sdl1, sdl2, sdl3)
        else
            # return a negative sign within the asphere
            if z_aspheric_boundary < 0 && (z_aspheric_value < z < 0.0)
                return -min(distance_to_aspheric, sdl1, sdl2, sdl3)
            elseif z_aspheric_boundary > 0 && (0.0 < z < z_aspheric_value)
                return -min(distance_to_aspheric, sdl1, sdl2, sdl3)
            else
                return min(distance_to_aspheric, sdl1, sdl2, sdl3)
            end
        end
    end
end

function sdf(surface::ConvexAsphericalSurfaceSDF{T}, point) where {T}
    p_local = _world_to_sdf(surface, point)
    # spherical logic is to have the z-axis as optical axis, so the aspheric code is written
    # with that convention in mind so everything matches ISO10110/textbook definitions.
    # We reinterpret the coordinates here, so xz is the transversal direction and y is the
    # optical axis.
    _pp = Point3{T}(p_local[1], p_local[3], p_local[2]) # xzy
    # rotate 2D sdf around the optical axis
    sdf_v = op_revolve_z(_pp,
        x -> convex_aspheric_surface_distance(
            x[1],
            x[2],
            1 / surface.radius,
            surface.conic_constant,
            surface.diameter,
            surface.coefficients,
            surface.max_sag
        ), zero(T))
    return sdf_v
end

function sdf(surface::ConcaveAsphericalSurfaceSDF{T}, point) where {T}
    p_local = _world_to_sdf(surface, point)
    # spherical logic is to have the z-axis as optical axis, so the aspheric code is written
    # with that convention in mind so everything matches ISO10110/textbook definitions.
    # We reinterpret the coordinates here, so xz is the transversal direction and y is the
    # optical axis.
    _pp = Point3{T}(p_local[1], p_local[3], p_local[2]) # xzy
    # rotate 2D sdf around the optical axis
    sdf_v = op_revolve_z(_pp,
        x -> concave_aspheric_surface_distance(
            x[1],
            x[2],
            1 / surface.radius,
            surface.conic_constant,
            surface.diameter,
            surface.coefficients,
            surface.max_sag
        ), zero(T))
    return sdf_v
end

function render_object!(axis, asp::AbstractAsphericalSurfaceSDF; color = :red)
    radius = asp.diameter / 2
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
- `α_coeffs` : The (even) aspheric coefficients, starting with A4.
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
- `coefficients::AbstractVector{T}`: A vector of even aspheric coefficients for higher-order corrections.

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
- `coefficients::AbstractVector`: A vector of even aspheric coefficients for higher-order corrections.
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

struct EvenAsphericSurface{T} <: AbstractRotationallySymmetricSurface{T}
    standard::StandardSurface{T}
    conic_constant::T
    coefficients::Vector{T}
end
function EvenAsphericSurface(radius::T1, diameter::T2, conic_constant::T3, coefficients::AbstractVector{T4}, mechanical_diameter::T5=diameter) where {T1, T2, T3, T4, T5}
    T = promote_type(T1, T2, T3, T4, T5)
    st = StandardSurface{T}(radius, diameter, mechanical_diameter)

    return EvenAsphericSurface{T}(
        st,
        conic_constant,
        coefficients
    )
end
radius(s::EvenAsphericSurface) = radius(s.standard)
diameter(s::EvenAsphericSurface) = diameter(s.standard)
mechanical_diameter(s::EvenAsphericSurface) = mechanical_diameter(s.standard)

function edge_thickness(s::EvenAsphericSurface, ::AbstractAsphericalSurfaceSDF)
    return abs(aspheric_equation(
        diameter(s) / 2,
        1 / radius(s),
        s.conic_constant,
        s.coefficients
    ))
end

function sdf(s::EvenAsphericSurface, ot::AbstractOrientationType)
    isinf(radius(s)) && return nothing

    return _sdf(s, ot)
end

function _sdf(s::EvenAsphericSurface, ::ForwardOrientation)
    front = if radius(s) > 0
        ConvexAsphericalSurfaceSDF(s.coefficients, radius(s), s.conic_constant, diameter(s))
    else
        ConcaveAsphericalSurfaceSDF(s.coefficients, radius(s), s.conic_constant, diameter(s))
    end

    return front
end

function _sdf(s::EvenAsphericSurface, ::BackwardOrientation)
    back = if radius(s) > 0
        ConcaveAsphericalSurfaceSDF(s.coefficients, radius(s), s.conic_constant, diameter(s))
    else
        ConvexAsphericalSurfaceSDF(s.coefficients, radius(s), s.conic_constant, diameter(s))
    end

    return back
end
