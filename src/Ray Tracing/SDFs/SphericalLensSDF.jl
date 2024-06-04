abstract type AbstractRotationallySymmetricLensSDF{T} <: AbstractSDF{T} end

"""
    sag(r::Real, l::Real)

Calculates the sag of a cut circle with radius `r` and chord length `l`
"""
sag(r::Real, l::Real) = r - sqrt(r^2 - 0.25 * l^2)

function check_sag(r, d)
    if abs(2r) < d 
        throw(ArgumentError("Radius of curvature (r = $(r)) must be ≥ than half the diameter (d = $(d)) or an illegal shape results!"))
    else
        return nothing
    end
end

"""
    AbstractSphericalSurfaceSDF{T} <: AbstractSDF{T}

An abstract type for SDF-based volumes which represent spherical lens surfaces, i.e. convex or concave.
"""
abstract type AbstractSphericalSurfaceSDF{T} <: AbstractRotationallySymmetricLensSDF{T} end

"""Returns the sagitta of the `AbstractSphericalSurfaceSDF`"""
sag(s::AbstractSphericalSurfaceSDF) = s.sag

"""Returns the radius of curvature of the `AbstractSphericalSurfaceSDF`"""
radius(s::AbstractSphericalSurfaceSDF) = s.radius

"""Returns the lens diameter of the `AbstractSphericalSurfaceSDF`"""
diameter(s::AbstractSphericalSurfaceSDF) = s.diameter

mutable struct ConcaveSphericalSurfaceSDF{T} <: AbstractSphericalSurfaceSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    radius::T
    diameter::T
    sag::T
end

function ConcaveSphericalSurfaceSDF(r::R, d::D) where {R, D}
    T = promote_type(R, D)
    check_sag(r, d)
    s = sag(r, d)
    return ConcaveSphericalSurfaceSDF{T}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        Point3{T}(0),
        r, d, s
    )
end

function sdf(css::ConcaveSphericalSurfaceSDF, point)
    p = _world_to_sdf(css, point)
    # cylinder sdf
    ps = p + Point3(0, sag(css)/2, 0)
    d = abs.(Point2(norm(Point2(ps[1], ps[3])), ps[2])) -
        Point2(diameter(css)/2, sag(css)/2)
    sdf1 = min(maximum(d), 0) + norm(max.(d, 0))
    # sphere sdf
    ps = p + Point3(0, radius(css), 0)
    sdf2 = norm(ps) - radius(css)
    return max(sdf1, -sdf2)
end

mutable struct ConvexSphericalSurfaceSDF{T} <: AbstractSphericalSurfaceSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    radius::T
    diameter::T
    sag::T
    height::T
end

function ConvexSphericalSurfaceSDF(radius::R, diameter::D) where {R, D}
    T = promote_type(R, D)
    check_sag(radius, diameter)
    _sag = sag(radius, diameter)
    cutoff_height = radius - _sag
    return ConvexSphericalSurfaceSDF{T}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        zeros(T, 3),
        radius,
        diameter,
        _sag,
        cutoff_height
    )
end

function sdf(css::ConvexSphericalSurfaceSDF, point)
    p = _world_to_sdf(css, point)
    # -p[2] to align surface with neg. y-axis
    q = Point2(norm(Point2(p[1], p[3])), -p[2] + css.height)
    s = max((css.height - radius(css)) * q[1]^2 + (diameter(css)/2)^2 * (css.height + radius(css) - 2 * q[2]),
        css.height * q[1] - diameter(css)/2 * q[2])
    if s < 0
        return norm(q) - radius(css)
    elseif q[1] < diameter(css)/2
        return css.height - q[2]
    else
        return norm(q - Point2(diameter(css)/2, css.height))
    end
end

"""
    ThinLensSDF(r1, r2, d=1inch)

Constructs a bi-convex thin lens SDF-based shape with:

- `r1 > 0`: radius of convex front
- `r2 > 0`: radius of convex back
- `d`: lens diameter, default value is one inch

The spherical surfaces are constructed flush.
"""
function ThinLensSDF(r1::L, r2::M, d::O = 1inch) where {L, M, O}
    check_sag(r1, d)
    check_sag(r2, d)
    # Create lens halves
    front = ConvexSphericalSurfaceSDF(r1, d)
    back = ConvexSphericalSurfaceSDF(r2, d)
    # Rotate and move back cut sphere
    zrotate3d!(back, π)
    return (front + back)
end

"""
    BiConvexLensSDF(r1, r2, l, d=1inch)

Constructs a cylindrical bi-convex lens SDF with:

- `r1` > 0: radius of convex front
- `r2` > 0: radius of convex back
- `l`: lens thickness
- `d`: lens diameter, default value is one inch

The spherical surfaces are constructed flush with the cylinder surface.
"""
function BiConvexLensSDF(r1::L, r2::M, l::N, d::O = 1inch) where {L, M, N, O}
    check_sag(r1, d)
    check_sag(r2, d)
    s1 = sag(r1, d)
    s2 = sag(r2, d)
    # Calculate length of cylindrical section
    l = l - (s1 + s2)
    front = ConvexSphericalSurfaceSDF(r1, d)
    back = ConvexSphericalSurfaceSDF(r2, d)
    mid = CylinderSDF(d / 2, l / 2)
    # Shift and rotate elements into position
    translate3d!(front, [0, -l / 2, 0])
    zrotate3d!(back, π)
    translate3d!(back, [0, l / 2, 0])
    return (front + mid + back)
end

function BiConcaveLensSDF(r1::L, r2::M, l::N, d::O = 1inch) where {L, M, N, O}
    # create segments
    front = ConcaveSphericalSurfaceSDF(r1, d)
    back = ConcaveSphericalSurfaceSDF(r2, d)
    mid = CylinderSDF(d / 2, l / 2)
    # Shift and rotate subtraction spheres into position
    translate3d!(front, [0, -l / 2, 0])
    zrotate3d!(back, π)
    translate3d!(back, [0, l / 2, 0])
    return (front + mid + back)
end

"""
    BiConcaveLensSDF(r1, r2, l, d=1inch)

Constructs a bi-concave lens SDF with:

- `r1` > 0: radius of concave front
- `r2` > 0: radius of convex back
- `l`: lens thickness
- `d`: lens diameter, default value is one inch
- `md`: mechanical lens diameter, adds an outer ring section to the lens, if `md` > `d`.

The spherical surfaces are constructed flush with the cylinder surface.
"""
function BiConcaveLensSDF(r1::L, r2::M, l::N, d::O, md::MD) where {L, M, N, O, MD}
    if md ≤ d
        throw(ArgumentError("Mech. diameter must be larger than lens diameter!"))
    end
    # generate ring-less shape
    shape = BiConcaveLensSDF(r1, r2, l, d)
    # add an outer ring
    l += sag(shape.sdfs[1]) + sag(shape.sdfs[3])
    ring = RingSDF(d/2, (md - d) / 2, l)
    shape += ring
    return shape
end

"""
    PlanoConvexLensSDF(r, l, d=1inch)

Constructs a plano-convex lens SDF with:

- `r` > 0: front radius
- `l`: lens thickness
- `d`: lens diameter, default value is one inch

The spherical surface is constructed flush with the cylinder surface.
"""
function PlanoConvexLensSDF(r::R, l::L, d::D = 1inch) where {R, L, D}
    _sag = sag(r, d)
    # Calculate length of cylindrical section
    l = l - _sag
    front = ConvexSphericalSurfaceSDF(r, d)
    back = CylinderSDF(d / 2, l / 2)
    # Shift and rotate cut spheres into position
    translate3d!(front, [0, -l / 2, 0])
    return (front + back)
end

function PlanoConcaveLensSDF(r::R, l::L, d::D = 1inch) where {R, L, D}
    front = CylinderSDF(d / 2, l / 2)
    back = ConcaveSphericalSurfaceSDF(r, d)
    # Shift and rotate elements into position
    zrotate3d!(back, π)
    translate3d!(back, [0, l / 2, 0])
    return (front + back)
end

"""
    PlanoConcaveLensSDF(r, l, d=1inch)

Constructs a plano-concave lens SDF with:

- `r` > 0: front radius
- `l`: lens thickness
- `d`: lens diameter, default value is one inch
- `md`: mechanical lens diameter, defaults to be identical to the lens diameter, Otherwise
        an outer ring section will be added to the lens, if `md` > `d`.

The spherical surface is constructed flush with the cylinder surface.
"""
function PlanoConcaveLensSDF(r::R, l::L, d::D, md::MD) where {R, L, D, MD}
    if md ≤ d
        throw(ArgumentError("Mech. diameter must be larger than lens diameter!"))
    end
    # generate ring-less shape
    shape = PlanoConcaveLensSDF(r, l, d)
    # add an outer ring
    _sag = sag(shape.sdfs[2])
    _l = l + _sag
    ring = RingSDF(d/2, (md - d) / 2, _l)
    translate3d!(ring, [0, l/2, 0])
    shape += ring
    return shape
end

function ConvexConcaveLensSDF(r1::R1, r2::R2, l::L, d::D = 1inch) where {R1, R2, L, D}
    # Create front and back surfaces
    front = ConvexSphericalSurfaceSDF(r1, d)
    back = ConcaveSphericalSurfaceSDF(r2, d)
    # Calculate length of cylindrical section
    _l = l - sag(front)
    mid = CylinderSDF(d / 2, _l / 2)
    # Move elements into position
    translate3d!(front, [0, -_l / 2, 0])
    translate3d!(back, [0, _l / 2, 0])
    zrotate3d!(back, π)
    return (front + mid + back)
end

"""
    ConvexConcaveLensSDF(r1, r2, l, d=1inch)

Constructs a positive/negative meniscus lens SDF with:

- `r1` > 0: radius of convex front
- `r2` > 0: radius of concave back
- `l`: lens thickness
- `d`: lens diameter, default value is one inch
- `md`: mechanical lens diameter, defaults to be identical to the lens diameter, Otherwise
        an outer ring section will be added to the lens, if `md` > `d`.

The spherical surface is constructed flush with the cylinder surface.
"""
function ConvexConcaveLensSDF(r1::R1, r2::R2, l::L, d::D, md::MD) where {R1, R2, L, D, MD}
    if md ≤ d
        throw(ArgumentError("Mech. diameter must be larger than lens diameter!"))
    end
    # generate ring-less shape
    shape = ConvexConcaveLensSDF(r1, r2, l, d)
    _l = 2*shape.sdfs[2].height + sag(shape.sdfs[3])
    ring = RingSDF(d/2, (md - d) / 2, _l)
    translate3d!(ring, [0, shape.sdfs[2].height, 0])
    return (shape + ring)
end

function render_object!(axis, css::ConcaveSphericalSurfaceSDF; color=:white)
    _radius = diameter(css)/2
    v = LinRange(0, 2π, 100)
    r = LinRange(1e-12, 1, 100) .^ (1/2) * _radius
    # Calculate beam surface at origin along y-axis, swap w and u
    y = sqrt.(radius(css)^2 .- r.^2) .- radius(css)
    u = y
    w = collect(r)
    # Close conture
    push!(u, 0) # push!(u, 0, 0)
    push!(w, _radius) # push!(w, radius, 1e-12)
    X = [w[i] * cos(v) for (i, u) in enumerate(u), v in v]
    Y = [u for u in u, v in v]
    Z = [w[i] * sin(v) for (i, u) in enumerate(u), v in v]
    # Transform into global coords
    R = css.dir
    P = css.pos
    Xt = R[1, 1] * X + R[1, 2] * Y + R[1, 3] * Z .+ P[1]
    Yt = R[2, 1] * X + R[2, 2] * Y + R[2, 3] * Z .+ P[2]
    Zt = R[3, 1] * X + R[3, 2] * Y + R[3, 3] * Z .+ P[3]
    render_surface!(axis, Xt, Yt, Zt; transparency = true, colormap = [color, color])
    return nothing
end

function render_object!(axis, css::ConvexSphericalSurfaceSDF; color=:white)
    v = LinRange(0, 2π, 100)
    r = LinRange(1e-12, 1, 100) .^ (1/2) * diameter(css)/2
    # Calculate beam surface at origin along y-axis, swap w and u
    y = -(sqrt.(radius(css)^2 .- r.^2) .- radius(css) .+ sag(css))
    u = y
    w = collect(r)
    # Close conture
    # push!(u, 0)
    # push!(w, 1e-12)
    X = [w[i] * cos(v) for (i, u) in enumerate(u), v in v]
    Y = [u for u in u, v in v]
    Z = [w[i] * sin(v) for (i, u) in enumerate(u), v in v]
    # Transform into global coords
    R = css.dir
    P = css.pos
    Xt = R[1, 1] * X + R[1, 2] * Y + R[1, 3] * Z .+ P[1]
    Yt = R[2, 1] * X + R[2, 2] * Y + R[2, 3] * Z .+ P[2]
    Zt = R[3, 1] * X + R[3, 2] * Y + R[3, 3] * Z .+ P[3]
    render_surface!(axis, Xt, Yt, Zt; transparency = true, colormap = [color, color])
    return nothing
end
