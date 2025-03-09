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
    # standard logic is to have the z-axis as optical axis, so the aspheric code is written
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
    # standard logic is to have the z-axis as optical axis, so the aspheric code is written
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
    generalized_lens_shape_constructor(r1, r2, l, d1, d2=d1; md1=d1, md2=d2,
                                     front_kind::Symbol = :spherical, front_k = 0, front_coeffs = nothing,
                                     back_kind::Symbol  = :spherical, back_k  = 0, back_coeffs  = nothing)

Constructs an `<: AbstractLensSDF` representing a rotationally symmetric lens whose
surfaces may be spherical or aspherical.

# Arguments
- `r1`: Front surface curvature (radius). A positive value indicates a convex-like front surface.
- `r2`: Back surface curvature (radius). A positive value indicates a concave-like back surface.
- `l`: Axial thickness of the lens.
- `d1`: Clear (optical) diameter for the front surface.
- `d2`: Clear (optical) diameter for the back surface (defaults to `d1`).
- `md1`: Mechanical diameter for the front; if greater than `d1`, an outer ring is added (defaults to `d1`).
- `md2`: Mechanical diameter for the back; if greater than `d2`, an outer ring is added (defaults to `d2`).

# Keyword Options (for surface type)
- `front_kind`: Either `:spherical` (default) or `:aspherical` for the front surface.
- `front_k`: Conic constant for the front surface (used if `front_kind == :aspherical`).
- `front_coeffs`: Vector of even aspheric coefficients for the front surface (used if `front_kind == :aspherical`).
- `back_kind`: Either `:spherical` (default) or `:aspherical` for the back surface.
- `back_k`: Conic constant for the back surface (used if `back_kind == :aspherical`).
- `back_coeffs`: Vector of even aspheric coefficients for the back surface (used if `back_kind == :aspherical`).

# Description
The function constructs the lens shape by first computing a central plano (cylindrical)
section with a diameter equal to the smaller of the two clear apertures (`min(d1, d2)`).
The curved front and back surfaces are constructed using the specified surface type
(spherical or aspherical) with their respective clear apertures. If aspheric surfaces
are used, the sag at the clear edge is computed via the aspheric equation.
Mechanical diameters (`md1`, `md2`) allow for adding an outer ring if they exceed
the corresponding optical diameters. If the remaining cylindrical section length
(after subtracting the sag contributions) is non-positive, the function falls back
to constructing a meniscus lens using `MeniscusLensSDF`.

Returns an object of type `<: AbstractLensSDF` representing the composite lens.
"""
function generalized_lens_shape_constructor(r1, r2, l, d1, d2 = d1;
        md1 = d1, md2 = d2,
        front_kind::Symbol = :spherical, front_k = 0, front_coeffs = nothing,
        back_kind::Symbol = :spherical, back_k = 0, back_coeffs = nothing)

    # Define effective (optical and mechanical) diameters:
    d_mid = min(d1, d2)
    md_mid = max(md1, md2)

    # Initialize remaining cylindrical section length.
    l0 = l

    # --- Front Surface ---
    if isinf(r1)
        front = nothing
    elseif front_kind === :spherical
        if r1 > 0
            front = BeamletOptics.ConvexSphericalSurfaceSDF(r1, d1)
            l0 -= BeamletOptics.sag(front)
        else
            front = BeamletOptics.ConcaveSphericalSurfaceSDF(abs(r1), d1)
        end
    elseif front_kind === :aspherical
        if front_coeffs === nothing
            throw(ArgumentError("front_coeffs must be provided for aspherical front surface"))
        end
        if r1 > 0
            front = ConvexAsphericalSurfaceSDF(front_coeffs, r1, front_k, d1)
            s_front = aspheric_equation(d1 / 2, 1 / r1, front_k, front_coeffs)
            l0 -= abs(s_front)
        else
            front = ConcaveAsphericalSurfaceSDF(front_coeffs, r1, front_k, d1)
        end
    else
        throw(ArgumentError("Unsupported front_kind: " * front_kind))
    end

    # --- Back Surface ---
    if isinf(r2)
        back = nothing
    elseif back_kind === :spherical
        if r2 > 0
            back = BeamletOptics.ConcaveSphericalSurfaceSDF(r2, d2)
            zrotate3d!(back, π)
        else
            back = BeamletOptics.ConvexSphericalSurfaceSDF(abs(r2), d2)
            zrotate3d!(back, π)
            l0 -= BeamletOptics.sag(back)
        end
    elseif back_kind === :aspherical
        if back_coeffs === nothing
            throw(ArgumentError("back_coeffs must be provided for aspherical back surface"))
        end
        if r2 > 0
            back = ConcaveAsphericalSurfaceSDF(back_coeffs, r2, back_k, d2)
            l0 -= thickness(back)
        else
            back = ConvexAsphericalSurfaceSDF(back_coeffs, r2, back_k, d2)
            s_back = aspheric_equation(d2 / 2, 1 / r2, back_k, back_coeffs)
            l0 -= abs(s_back)
        end
    else
        throw(ArgumentError("Unsupported back_kind: " * back_kind))
    end

    # --- Use MeniscusLensSDF if cylinder length is non-positive ---
    if l0 ≤ 0
        if sign(r1) == sign(r2)
            return generalized_meniscus_lens_sdf(r1, r2, l, d1, md1, d2, md2;
                front_kind = front_kind, front_k = front_k, front_coeffs = front_coeffs,
                back_kind = back_kind, back_k = back_k, back_coeffs = back_coeffs)
        end
        throw(ArgumentError("Lens parameters lead to cylinder section length of ≤ 0, use ThinLens instead."))
    end

    # --- Construct central plano surface and add spherical/aspherical parts ---
    mid = PlanoSurfaceSDF(l0, d_mid)
    if front !== nothing
        translate3d!(mid, [0, thickness(front), 0])
        mid += front
    end
    if back !== nothing
        translate3d!(back, [0, thickness(mid) + thickness(back), 0])
        mid += back
    end

    # --- Add mechanical ring if md_mid > d_mid ---
    if md_mid > d_mid
        ring_thickness = thickness(mid)
        ring_center = position(mid)[2] + ring_thickness / 2
        if front !== nothing
            if front_kind === :spherical
                s = abs(sag(front))
            elseif front_kind === :aspherical
                s = aspheric_equation(d1 / 2, 1 / r1, front_k, front_coeffs)
            else
                throw(ArgumentError("Unsupported front_kind: " * front_kind))
            end
            ring_thickness -= s
            ring_center += s / 2
        end
        if back !== nothing
            if back_kind === :spherical
                s = abs(sag(back))
            elseif back_kind === :aspherical
                s = aspheric_equation(d2 / 2, 1 / r2, back_k, back_coeffs)
            else
                throw(ArgumentError("Unsupported back_kind: " * back_kind))
            end
            ring_thickness += s
            ring_center += s / 2
        end
        ring = RingSDF(d_mid / 2, (md_mid - d_mid) / 2, ring_thickness)
        translate3d!(ring, [0, ring_center, 0])
        mid += ring
    elseif md_mid < d_mid
        @warn "Mechanical diameter is less than clear aperture; parameter md has been ignored."
    end

    return mid
end

"""
    generalized_meniscus_lens_sdf(r1, r2, l, d1, d2=d1; md1=d1, md2=d2,
                                  front_kind::Symbol = :spherical, front_k = 0, front_coeffs = nothing,
                                  back_kind::Symbol  = :spherical, back_k  = 0, back_coeffs  = nothing)

Constructs an [`AbstractLensSDF`] meniscus lens whose surfaces can be either spherical or aspherical.
This function supports independent specification of the clear (optical) diameters and mechanical diameters
for the front and back surfaces. The clear diameter for the front surface is given by `d1` and for the
back surface by `d2` (defaulting to `d1`). Similarly, the mechanical diameters are specified by `md1` and
`md2` (defaulting to `d1` and `d2`, respectively). If a mechanical diameter exceeds the corresponding clear
diameter, an outer ring is added on that side. The curvature signs follow ISO 10110 (a positive value means
the center lies to the right). If the computed cylindrical section thickness becomes non-positive, an error is thrown.

# Arguments
- `r1`: Front surface curvature (radius). A positive value produces a convex-like front surface.
- `r2`: Back surface curvature (radius). A positive value produces a concave-like back surface.
- `l`: Total axial thickness of the lens.
- `d1`: Clear (optical) diameter of the front surface.
- `d2`: Clear (optical) diameter of the back surface (default: `d1`).
- `md1`: Mechanical diameter for the front surface (default: `d1`).
- `md2`: Mechanical diameter for the back surface (default: `d2`).

# Keyword Options (for surface types)
- `front_kind`: Surface type for the front; either `:spherical` (default) or `:aspherical`.
- `front_k`: Conic constant for the front surface (used if `front_kind == :aspherical`).
- `front_coeffs`: Vector of even aspheric coefficients for the front surface (used if `front_kind == :aspherical`).
- `back_kind`: Surface type for the back; either `:spherical` (default) or `:aspherical`.
- `back_k`: Conic constant for the back surface (used if `back_kind == :aspherical`).
- `back_coeffs`: Vector of even aspheric coefficients for the back surface (used if `back_kind == :aspherical`).

# Description
This function constructs a composite meniscus lens by:
  1. Building a front surface using `r1` and clear diameter `d1` (with the specified surface model),
  2. Building a back surface using `r2` and clear diameter `d2` (with the specified surface model),
  3. Constructing a central plano (cylindrical) section whose diameter is the minimum of `d1` and `d2`,
     with an effective thickness adjusted by subtracting the sag from the front and adding the sag from the back.
Mechanical diameters (`md1` and `md2`) are used to determine if an outer ring should be added; the effective
mechanical diameter is taken as `max(md1, md2)`, and if it exceeds the clear aperture (i.e. `min(d1, d2)`),
an outer ring is appended. If the resulting cylindrical section length is non-positive, an error is thrown.

Returns a meniscus lens object of type [`MeniscusLensSDF{T, S1, S2}`] representing the composite lens.
"""
function generalized_meniscus_lens_sdf(
        r1::R1, r2::R2, l::L, d1::D, d2::D = d1, md1::MD = d1, md2::MD = d2;
        front_kind::Symbol = :spherical, front_k = 0, front_coeffs = nothing,
        back_kind::Symbol = :spherical, back_k = 0, back_coeffs = nothing) where {
        R1, R2, L, D, MD}
    T = promote_type(R1, R2, L, D)
    # Determine orientation: for left‐facing (r₁, r₂ > 0) front is convex and back is concave;
    # for right‐facing (r₁, r₂ < 0) front is concave and back is convex.
    if sign(r1) == sign(r2) > 0
        orientation = :left_facing
        convex_val = r1
        concave_val = r2
    elseif sign(r1) == sign(r2) < 0
        orientation = :right_facing
        convex_val = abs(r2)
        concave_val = abs(r1)
    else
        throw(ArgumentError("Invalid sign combination for r₁ and r₂"))
    end

    # Compute sag values at the clear aperture:
    # Use d1 for the front and d2 for the back.
    convex_sag = front_kind == :spherical ?
                 sag(convex_val, d1) :
                 abs(aspheric_equation(d1 / 2, 1 / convex_val, front_k, front_coeffs))
    concave_sag = back_kind == :spherical ?
                  sag(concave_val, d2) :
                  abs(aspheric_equation(d2 / 2, 1 / concave_val, back_k, back_coeffs))

    cylinder_l = l - convex_sag + concave_sag
    if cylinder_l ≤ 0
        throw(ErrorException("Lens parameters lead to zero lens edge thickness"))
    end

    # Spawn sub-shapes.
    if orientation == :left_facing
        front_surface = front_kind == :spherical ?
                        ConvexSphericalSurfaceSDF(convex_val, d1) :
                        ConvexAsphericalSurfaceSDF(front_coeffs, convex_val, front_k, d1)
        back_surface = back_kind == :spherical ?
                       SphereSDF(concave_val) :
                       ConcaveAsphericalSurfaceSDF(back_coeffs, concave_val, back_k, d2)
    else  # right-facing: front is concave, back is convex.
        front_surface = front_kind == :spherical ?
                        SphereSDF(concave_val) :
                        ConcaveAsphericalSurfaceSDF(front_coeffs, concave_val, front_k, d2)
        back_surface = back_kind == :spherical ?
                       ConvexSphericalSurfaceSDF(convex_val, d1) :
                       ConvexAsphericalSurfaceSDF(back_coeffs, convex_val, back_k, d1)
    end

    # Use effective optical diameter d_mid = min(d1, d2)
    d_mid = min(d1, d2)
    cylinder = PlanoSurfaceSDF(cylinder_l, d_mid)

    # Position the sub-shapes.
    if orientation == :left_facing
        translate3d!(cylinder, [0, thickness(front_surface), 0])
        translate3d!(back_surface, [0, concave_val + l, 0])
        convex_shape = front_surface
        concave_shape = back_surface
    else
        translate3d!(back_surface, [0, -concave_val, 0])
        translate3d!(cylinder, [0, -concave_sag, 0])
        zrotate3d!(front_surface, π)
        translate3d!(front_surface, [0, thickness(cylinder) - concave_sag + convex_sag, 0])
        convex_shape = back_surface
        concave_shape = front_surface
    end

    # Build the composite meniscus lens.
    lens = MeniscusLensSDF{T, typeof(convex_shape), typeof(concave_shape)}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        Point3{T}(0),
        convex_shape,
        cylinder,
        concave_shape,
        l
    )

    # Effective mechanical diameter: md_mid = max(md1, md2)
    md_mid = max(md1, md2)
    if md_mid > d_mid
        # For a meniscus lens the overall optical thickness is l.
        ring_thickness = l
        ring_center = l / 2
        # Adjust for the sag of the front surface (use d1).
        if front_surface !== nothing
            if front_kind == :spherical
                s = abs(sag(front_surface))
            elseif front_kind == :aspherical
                s = abs(aspheric_equation(d1 / 2, 1 / r1, front_k, front_coeffs))
            else
                throw(ArgumentError("Unsupported front_kind: $front_kind"))
            end
            ring_thickness += (s < 0 ? abs(s) : 0.0)
            ring_center += (s < 0 ? s / 2 : 0.0)
        end
        # Adjust for the sag of the back surface (use d2).
        if back_surface !== nothing
            if back_kind == :spherical
                s = abs(sag(back_surface))
            elseif back_kind == :aspherical
                s = abs(aspheric_equation(d2 / 2, 1 / r2, back_k, back_coeffs))
            else
                throw(ArgumentError("Unsupported back_kind: $back_kind"))
            end
            ring_thickness += (s < 0 ? abs(s) : 0.0)
            ring_center += (s < 0 ? s / 2 : 0.0)
        end
        ring = RingSDF(d_mid / 2, (md_mid - d_mid) / 2, ring_thickness)
        translate3d!(ring, [0, ring_center, 0])
        lens += ring
    elseif md_mid < d_mid
        @warn "Mechanical diameter is less than clear aperture; parameter md has been ignored."
    end

    return lens
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
