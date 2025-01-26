# New abstract type for extended aspherical surfaces
abstract type AbstractExtendedAsphericalSurfaceSDF{T} <: AbstractAsphericalSurfaceSDF{T} end

"""
    ExtendedConvexAsphericalSurfaceSDF

Represents a convex-like aspheric surface that can include a more general polynomial expansion
with both even and odd terms.

- `coefficients`: Vector of polynomial coefficients for the aspheric departure.
  The terms are assumed to start from A2, A3, ... so that the aspheric departure is:
  z(r) = c*r² / [1 + sqrt(1 - (1+k)*c²*r²)] + sum(A_n * r^n, n≥2)
- `radius`: radius of curvature (positive for convex surfaces)
- `conic_constant`: conic constant (k)
- `diameter`: lens diameter
"""
mutable struct ExtendedConvexAsphericalSurfaceSDF{T} <: AbstractExtendedAsphericalSurfaceSDF{T}
    coefficients::Vector{T}
    radius::T
    conic_constant::T
    diameter::T
    pos::Point3{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
end

function ExtendedConvexAsphericalSurfaceSDF(coefficients::Vector{T}, radius::T, conic_constant::T, diameter::T) where {T}
    return ExtendedConvexAsphericalSurfaceSDF{T}(
        coefficients, radius, conic_constant, diameter,
        Point3{T}(0),
        SMatrix{3,3}(one(T)*I),
        SMatrix{3,3}(one(T)*I)
    )
end

"""
    ExtendedConcaveAsphericalSurfaceSDF

Represents a concave-like aspheric surface that can include a more general polynomial expansion
with both even and odd terms.

- `coefficients`: Vector of polynomial coefficients for the aspheric departure.
  The terms are assumed to start from A2, A3, ... so that the aspheric departure is:
  z(r) = c*r² / [1 + sqrt(1 - (1+k)*c²*r²)] + sum(A_n * r^n, n≥2)
- `radius`: radius of curvature (negative for concave surfaces)
- `conic_constant`: conic constant (k)
- `diameter`: clear aperture diameter of the lens
- `mechanical_diameter`: mechanical lens diameter (can be larger than the clear aperture)
"""
mutable struct ExtendedConcaveAsphericalSurfaceSDF{T} <: AbstractExtendedAsphericalSurfaceSDF{T}
    coefficients::Vector{T}
    radius::T
    conic_constant::T
    diameter::T
    mechanical_diameter::T
    pos::Point3{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
end

function ExtendedConcaveAsphericalSurfaceSDF(coefficients::V, radius::T, conic_constant::T, diameter::T, mechanical_diameter::T = diameter) where {T, V <: AbstractVector{T}}
    return ExtendedConcaveAsphericalSurfaceSDF{T}(
        coefficients, radius, conic_constant, diameter, mechanical_diameter,
        Point3{T}(0),
        SMatrix{3,3}(one(T)*I),
        SMatrix{3,3}(one(T)*I)
    )
end

"""
    extended_aspheric_equation(r, c, k, coeffs)

Generalized aspheric surface equation that can handle both even and odd polynomial terms.

The base form of the asphere (ISO 10110 style) is:
z(r) = c*r² / [1 + sqrt(1 - (1+k)*c²*r²)] + Σ A_n * r^n, for n≥2.

- `r`: radial coordinate
- `c`: curvature (1/radius)
- `k`: conic constant
- `coeffs`: vector of polynomial coefficients A_n, starting from n=2. That is,
            coeffs[1] = A₂, coeffs[2] = A₃, ..., so that the exponent corresponds to index+1.

Returns NaN if the square root argument is negative.
"""
function extended_aspheric_equation(r, c, k, coeffs)
    r2 = r^2
    sqrt_arg = 1 - (1 + k)*c^2*r2
    if sqrt_arg < 0
        return NaN
    end

    # polynomial expansion:
    # If coeffs = [A₂, A₃, A₄, ...],
    # term index i in coeffs corresponds to A_(i+1+1), i.e. A₂ for i=1, A_3 for i=2, etc.
    # The exponent for A_n term is n, which equals i+1 (since i starts at 1 for A₂).
    poly_term = zero(eltype(coeffs))
    for (i, a) in enumerate(coeffs)
        n = i + 1  # i=1 -> n=2, i=2 -> n=3, ...
        poly_term += a * r^n
    end

    return (c*r2)/(1 + sqrt(sqrt_arg)) + poly_term
end

"""
    gradient_extended_aspheric_equation(r, c, k, coeffs)

Compute the gradient with respect to `z` and `r` for the extended aspheric surface.
This is used in the SDF computation to determine local surface orientation.

Follows similar logic as `gradient_aspheric_equation` but now includes terms for all polynomial orders.

Returns NaN if the square root argument is negative.
"""
function gradient_extended_aspheric_equation(r, c, k, coeffs)
    Ri = 1/c
    sqrt_arg = 1 - r^2*(1+k)/Ri^2
    sqrt_arg < 0 && return NaN

    # Derivative of the base conic part wrt r:
    # z(r) base = c*r² / [1 + sqrt(1-(1+k)*c²*r²)]
    # Let’s define an intermediate variable for clarity:
    denom = (Ri*(√(sqrt_arg)+1))
    base_deriv = 2*r/denom + (r^3*(1+k)/(Ri^3*√(sqrt_arg)*(√(√(sqrt_arg))+1)^2))
    # Note: Care must be taken with the chain rule on these terms. If necessary,
    # verify the exact derivative steps as in standard asphere math references
    # (see e.g. Malacara, "Handbook of Optical Engineering" [CRC Press, 2001] for standard forms).

    # Derivative of the polynomial terms:
    # For coeffs[i] = A_(i+1+1), the exponent is n = i+1+1 = i+1 starting from i=1:
    # d/d(r) Σ A_n * r^n = Σ n*A_n * r^(n-1)
    poly_deriv = zero(eltype(coeffs))
    for (i, a) in enumerate(coeffs)
        n = i + 1
        poly_deriv += n*a*r^(n-1)
    end

    # The gradient direction in 2D (r,z) context is typically given as dZ/dR, but we return a point for compatibility.
    # We treat (-sum of partial derivatives wrt surface coordinate, 1) similar to original code,
    # so the gradient vector in the R-Z plane is perpendicular to the surface.
    # The logic is consistent with the original gradient code.

    total_deriv = base_deriv + poly_deriv
    return Point2(-total_deriv, 1)
end

extended_aspheric_equation(r::Real, a::AbstractExtendedAsphericalSurfaceSDF) = extended_aspheric_equation(r, 1/a.radius, a.conic_constant, a.coefficients)

# 2D distance functions for extended aspheres, following similar structure as original code.
function extended_convex_aspheric_surface_distance(r, z, c, k, d, coeffs)
    r2 = r^2
    r2_bound = (d/2)^2
    z_aspheric_value = extended_aspheric_equation(r, c, k, coeffs)
    grad_z = gradient_extended_aspheric_equation(r, c, k, coeffs)

    z_aspheric_boundary = extended_aspheric_equation(d/2, c, k, coeffs)
    grad_z_boundary = gradient_extended_aspheric_equation(d/2, c, k, coeffs)

    if isnan(z_aspheric_value) || isnan(grad_z) || r2 > r2_bound
        distance_to_boundary = sqrt((r - sign(r)*d/2)^2 + (z - z_aspheric_boundary)^2)
        return distance_to_boundary / norm(grad_z_boundary)
    end

    distance_to_aspheric = abs(z - z_aspheric_value) / norm(grad_z)
    # perimeter line segment
    a, b = Point2(d/2, z_aspheric_boundary), Point2(-d/2, z_aspheric_boundary)
    sdl = sd_line_segment(Point2(r,z), a, b) / norm(grad_z_boundary)

    if sign(c)*z_aspheric_value < sign(c)*z < sign(c)*z_aspheric_boundary
        return -min(sdl, distance_to_aspheric)
    else
        return min(sdl, distance_to_aspheric)
    end
end

function extended_concave_aspheric_surface_distance(r, z, c, k, d, coeffs)
    r2 = r^2
    r2_bound = (d/2)^2
    z_aspheric_value = extended_aspheric_equation(r, c, k, coeffs)
    grad_z = gradient_extended_aspheric_equation(r, c, k, coeffs)

    z_aspheric_boundary = extended_aspheric_equation(d/2, c, k, coeffs)
    grad_z_boundary = gradient_extended_aspheric_equation(d/2, c, k, coeffs)

    if isnan(z_aspheric_value) || isnan(grad_z)
        distance_to_boundary = sqrt((r - sign(r)*d/2)^2 + (z - z_aspheric_boundary)^2)
        return distance_to_boundary / norm(grad_z_boundary)
    end

    distance_to_aspheric = abs(z - z_aspheric_value) / norm(grad_z)
    p = Point2(r,z)

    a1, b1 = Point2(d/2, z_aspheric_boundary), Point2(d/2, 0.0)
    sdl1 = sd_line_segment(p, a1, b1)/norm(grad_z_boundary)
    a2, b2 = Point2(d/2, 0.0), Point2(-d/2, 0.0)
    sdl2 = sd_line_segment(p, a2, b2)/norm(grad_z_boundary)
    a3, b3 = Point2(-d/2, 0.0), Point2(-d/2, z_aspheric_boundary)
    sdl3 = sd_line_segment(p, a3, b3)/norm(grad_z_boundary)

    if r2 > r2_bound
        return min(sdl1, sdl2, sdl3)
    else
        if z_aspheric_value < z < 0.0
            return -min(distance_to_aspheric, sdl1, sdl2, sdl3)
        else
            return min(distance_to_aspheric, sdl1, sdl2, sdl3)
        end
    end
end

function sdf(surface::ExtendedConvexAsphericalSurfaceSDF{T}, point) where {T}
    p_local = _world_to_sdf(surface, point)
    _pp = Point3{T}(p_local[1], p_local[3], p_local[2]) # xzy
    return op_revolve(_pp,
        x->extended_convex_aspheric_surface_distance(
            x[1],
            x[2],
            1/surface.radius,
            surface.conic_constant,
            surface.diameter,
            surface.coefficients
        ), zero(T))
end

function sdf(surface::ExtendedConcaveAsphericalSurfaceSDF{T}, point) where {T}
    p_local = _world_to_sdf(surface, point)
    _pp = Point3{T}(p_local[1], p_local[3], p_local[2]) # xzy
    return op_revolve(_pp,
        x->extended_concave_aspheric_surface_distance(
            x[1],
            x[2],
            1/surface.radius,
            surface.conic_constant,
            surface.diameter,
            surface.coefficients
        ), zero(T))
end

# Rendering function for extended surfaces similar to previous ones.
function render_object!(axis, asp::AbstractExtendedAsphericalSurfaceSDF; color=:blue)
    radius = asp.diameter/2
    v = LinRange(0, 2π, 100)
    r = LinRange(1e-12, radius, 50)
    y = extended_aspheric_equation.(r, Ref(asp))
    u = y
    w = collect(r)
    if isa(asp, ExtendedConvexAsphericalSurfaceSDF)
        push!(u, u[end])
        push!(w, 1e-12)
    elseif isa(asp, ExtendedConcaveAsphericalSurfaceSDF)
        push!(u, 0, 0)
        push!(w, radius, 1e-12)
    else
        @warn "No suitable render function for $(typeof(asp))"
        return nothing
    end

    X = [w[i] * cos(v) for (i, u) in enumerate(u), v in v]
    Y = [u for u in u, v in v]
    Z = [w[i] * sin(v) for (i, u) in enumerate(u), v in v]
    R = asp.dir
    P = asp.pos
    Xt = R[1,1]*X + R[1,2]*Y + R[1,3]*Z .+ P[1]
    Yt = R[2,1]*X + R[2,2]*Y + R[2,3]*Z .+ P[2]
    Zt = R[3,1]*X + R[3,2]*Y + R[3,3]*Z .+ P[3]
    render_surface!(axis, Xt, Yt, Zt; transparency = true, colormap = [color, color])
    return nothing
end
