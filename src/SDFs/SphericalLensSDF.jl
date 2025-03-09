"""
    AbstractLensSDF

A class of [`AbstractSDF`](@ref)s that can be used to represent rotationally symmetric lens surfaces.
It is implicity assumed that all surfaces are represented by **closed volumes** for ray-tracing correctness.

# Implementation reqs.

Subtypes of `AbstractLensSDF` must implement the following:

## Functions:

- `thickness`: this function returns the material thickness of the element along its symmetry axis
- `diameter`: this function returns the outer diameter of the element

!!! note "Shape orientation"
    For easy compatibility between subtypes, the follwing requirements should be fulfilled:
    1. Symmetry axis aligned onto the y-axis
    2. Surface contour aligned towards negative y-values
    3. Surface point with `min(y)` should satisfy `min(y) = 0` on the symmetry axis
"""
abstract type AbstractLensSDF{T} <: AbstractSDF{T} end

"""Returns the outer bounding diameter of the `AbstractLensSDF`"""
diameter(s::AbstractLensSDF) = s.diameter

"""Returns the on-axis thickness of the `AbstractLensSDF`"""
thickness(s::AbstractLensSDF) = s.thickness

"""
    PlanoSurfaceSDF

[`AbstractSDF`](@ref)-based representation of two flat optical surfaces, i.e. equivalent to the [`CylinderSDF`](@ref).
When constructed, it is assumed that the first flat surface lies at the origin and the optical axis is aligned with the positive `y`-axis.

## Fields:

- `diameter`: the outer diameter of the circular flat lens surface
- `thickness`: the distance between the flat surfaces
"""
mutable struct PlanoSurfaceSDF{T} <: AbstractLensSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    thickness::T
    diameter::T
end

function PlanoSurfaceSDF(thickness::T, diameter::D) where {T, D}
    F = promote_type(T, D)
    return PlanoSurfaceSDF{F}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        Point3{T}(0),
        thickness,
        diameter
    )
end

function sdf(ps::PlanoSurfaceSDF{T}, point) where T
    p = _world_to_sdf(ps, point)
    d = abs.(Point2(norm(Point2(p[1], p[3])), p[2] - thickness(ps)/2)) -
        Point2(diameter(ps)/2, thickness(ps)/2)
    return min(maximum(d), zero(T)) + norm(max.(d, zero(T)))
end

"""
    SphereSDF

Implements the SDF of a perfect sphere. Orientation is fixed to unity matrix.
"""
mutable struct SphereSDF{T} <: AbstractLensSDF{T}
    pos::Point3{T}
    radius::T
end

diameter(s::SphereSDF) = 2s.radius
thickness(s::SphereSDF) = 2s.radius

SphereSDF(r::T) where {T} = SphereSDF{T}(zeros(T, 3), r)

orientation(::SphereSDF{T}) where {T} = SArray{Tuple{3,3}, T}(I)
transposed_orientation(::SphereSDF{T}) where {T} = SArray{Tuple{3,3}, T}(I)
orientation!(::SphereSDF, ::Any) = nothing

function sdf(sphere::SphereSDF, point)
    p = _world_to_sdf(sphere, point)
    return norm(p) - sphere.radius
end

"""
    AbstractSphericalSurfaceSDF{T} <: AbstractSDF{T}

An abstract type for SDF-based volumes which represent spherical lens surfaces, i.e. [`ConvexSphericalSurfaceSDF`](@ref) or [`ConcaveSphericalSurfaceSDF`](@ref).

# Implementation reqs.

Subtypes of `AbstractSphericalSurfaceSDF` should implement all supertype reqs. as well as the following:

## Fields:

- `radius`: the radius of curvature
- `diameter`: the lens outer diameter
- `sag`: the lens sagitta

## Lens construction

It is intended that practical lens shapes are constructed from `AbstractSphericalSurfaceSDF`s using the [`UnionSDF`](@ref) type.
"""
abstract type AbstractSphericalSurfaceSDF{T} <: AbstractLensSDF{T} end

"""Returns the sagitta of the `AbstractSphericalSurfaceSDF`"""
sag(s::AbstractSphericalSurfaceSDF) = s.sag

"""Returns the radius of curvature of the `AbstractSphericalSurfaceSDF`"""
radius(s::AbstractSphericalSurfaceSDF) = s.radius

"""
    ConcaveSphericalSurfaceSDF

[`AbstractSDF`](@ref)-based representation of a concave spherical lens surface.
When constructed, it is assumed that the plano-surface lies at the origin and the optical axis is aligned with the `y`-axis.
The concave surface is orientated towards negative y-values for `R > 0` and vice versa.

## Fields:

- `radius`: the radius of curvature of the convex spherical surface.
- `diameter`: the outer diameter of the lens surface
- `sag`: the sagitta of the opposing convex shape
"""
mutable struct ConcaveSphericalSurfaceSDF{T} <: AbstractSphericalSurfaceSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    radius::T
    diameter::T
    sag::T
end

thickness(::ConcaveSphericalSurfaceSDF{T}) where T = zero(T)

"""
    ConcaveSphericalSurfaceSDF(radius, diameter)

Constructs a [`ConcaveSphericalSurfaceSDF`](@ref) with a specific `radius` of curvature and lens outer `diameter`.
"""
function ConcaveSphericalSurfaceSDF(radius::R, diameter::D) where {R, D}
    T = promote_type(R, D)
    check_sag(radius, diameter)
    s = sag(radius, diameter)
    return ConcaveSphericalSurfaceSDF{T}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        Point3{T}(0),
        radius, diameter, s
    )
end

function sdf(css::ConcaveSphericalSurfaceSDF{T}, point) where T
    p = _world_to_sdf(css, point)
    # cylinder sdf
    ps = p + Point3(zero(T), sag(css)/2, zero(T))
    d = abs.(Point2(norm(Point2(ps[1], ps[3])), ps[2])) -
        Point2(diameter(css)/2, sag(css)/2)
    sdf1 = min(maximum(d), zero(T)) + norm(max.(d, zero(T)))
    # sphere sdf
    ps = p + Point3(zero(T), radius(css), zero(T))
    sdf2 = norm(ps) - radius(css)
    return max(sdf1, -sdf2)
end

"""
    ConvexSphericalSurfaceSDF

[`AbstractSDF`](@ref)-based representation of a convex spherical lens surface.
When constructed, it is assumed that the plano-surface lies at the origin and the optical axis is aligned with the `y`-axis.
The convex surface is orientated towards negative y-values for `R > 0` and vice versa.

## Fields:

- `radius`: the radius of curvature of the concave spherical surface.
- `diameter`: the outer diameter of the lens surface
- `sag`: the sagitta of the convex shape
- `height`: the sphere cutoff height, see also [`CutSphereSDF`](@ref)
"""
mutable struct ConvexSphericalSurfaceSDF{T} <: AbstractSphericalSurfaceSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    radius::T
    diameter::T
    sag::T
    height::T
end

thickness(css::ConvexSphericalSurfaceSDF) = sag(css)

"""
    ConvexSphericalSurfaceSDF(radius, diameter)

Constructs a [`ConvexSphericalSurfaceSDF`](@ref) with a specific `radius` of curvature and lens outer `diameter`.
"""
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
    q = Point2(norm(Point2(p[1], p[3])), -p[2] + radius(css))
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
    # Create lens halves
    front = ConvexSphericalSurfaceSDF(r1, d)
    back = ConvexSphericalSurfaceSDF(r2, d)
    # Rotate and move back cut sphere
    translate3d!(back, [0, thickness(front) + thickness(back), 0])
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
    mid = PlanoSurfaceSDF(l, d)
    # Shift and rotate elements into position
    translate3d!(mid, [0, thickness(front), 0])
    zrotate3d!(back, π)
    translate3d!(back, [0, thickness(front) + thickness(mid) + thickness(back), 0])
    return (front + mid + back)
end

function BiConcaveLensSDF(r1::L, r2::M, l::N, d::O = 1inch) where {L, M, N, O}
    # create segments
    front = ConcaveSphericalSurfaceSDF(r1, d)
    back = ConcaveSphericalSurfaceSDF(r2, d)
    mid = PlanoSurfaceSDF(l, d)
    # Shift and rotate subtraction spheres into position
    zrotate3d!(back, π)
    translate3d!(back, [0, thickness(mid), 0])
    return (front + mid + back)
end

function check_md(d, md)
    if md ≤ d
        throw(ArgumentError("Mech. diameter must be larger than lens diameter!"))
    end
    return nothing
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
    check_md(d, md)
    # generate ring-less shape
    shape = BiConcaveLensSDF(r1, r2, l, d)
    # add an outer ring
    l0 = l + sag(shape.sdfs[1]) + sag(shape.sdfs[3])
    ring = RingSDF(d/2, (md - d) / 2, l0)
    translate3d!(ring, [0, l/2, 0])
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
    back = PlanoSurfaceSDF(l, d)
    # Shift and rotate cut spheres into position
    translate3d!(back, [0, thickness(front), 0])
    return (front + back)
end

function PlanoConcaveLensSDF(r::R, l::L, d::D = 1inch) where {R, L, D}
    front = PlanoSurfaceSDF(l, d)
    back = ConcaveSphericalSurfaceSDF(r, d)
    # Shift and rotate elements into position
    zrotate3d!(back, π)
    translate3d!(back, [0, thickness(front), 0])
    return (front + back)
end

"""
    PlanoConcaveLensSDF(r, l, d=1inch)

Constructs a plano-concave lens SDF with:

- `r` > 0: front radius
- `l`: lens thickness
- `d`: lens diameter, default value is one inch
- `md`: mechanical lens diameter, must be > d

The spherical surface is constructed flush with the cylinder surface.
"""
function PlanoConcaveLensSDF(r::R, l::L, d::D, md::MD) where {R, L, D, MD}
    check_md(d, md)
    # generate ring-less shape
    shape = PlanoConcaveLensSDF(r, l, d)
    # add an outer ring
    _sag = sag(shape.sdfs[2])
    _l = l + _sag
    ring = RingSDF(d/2, (md - d) / 2, _l)
    translate3d!(ring, [0, _l/2, 0])
    shape += ring
    return shape
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
    y = -(sqrt.(radius(css)^2 .- r.^2) .- radius(css))
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

"""
    SphericalLensShapeConstructor(r1, r2, l, d)

Creates an [`AbstractLensSDF`](@ref)-based spherical [`Lens`](@ref) shape.

!!! info "Radius of curvature (ROC) sign"
    The ROC is defined to be positive if the center is to the right of the surface. Otherwise it is negative.
    Plano-surfaces are created through the use of `Inf`.

# Input

- `r1`: front radius
- `r2`: back radius
- `l`: lens thickness
- `d`: lens diameter, default is one inch
"""
function SphericalLensShapeConstructor(r1, r2, l, d)
    # length tracking variable
    l0 = l
    # determine front shape
    if isinf(r1)
        front = nothing
    elseif r1 > 0
        front = BeamletOptics.ConvexSphericalSurfaceSDF(r1, d)
        l0 -= BeamletOptics.sag(front)
    else
        front = BeamletOptics.ConcaveSphericalSurfaceSDF(abs(r1), d)
    end
    # determine back shape
    if isinf(r2)
        back = nothing
    elseif r2 > 0
        back = BeamletOptics.ConcaveSphericalSurfaceSDF(r2, d)
        zrotate3d!(back, π)
    else
        back = BeamletOptics.ConvexSphericalSurfaceSDF(abs(r2), d)
        zrotate3d!(back, π)
        l0 -= BeamletOptics.sag(back)
    end
    # test if cylinder is legal
    if l0 ≤ 0
        # Catch MeniscusLensSDF
        if sign(r1) == sign(r2)
            return MeniscusLensSDF(r1, r2, l, d)
        end
        throw(ArgumentError("Lens parameters lead to cylinder section length of ≤ 0, use ThinLens instead."))
    end
    # design plano-plano surface, add spherical parts
    mid = BeamletOptics.PlanoSurfaceSDF(l0, d)
    if !isnothing(front)
        translate3d!(mid, [0, BeamletOptics.thickness(front), 0])
        mid += front
    end
    if !isnothing(back)
        translate3d!(back, [0,  BeamletOptics.thickness(mid) + BeamletOptics.thickness(back), 0])
        mid += back
    end
    # return shape
    return mid
end

struct StandardSurface{T} <: AbstractRotationallySymmetricSurface{T}
    radius::T
    diameter::T
    mechanical_diameter::T
end
StandardSurface(radius::T, diameter::T) where T = StandardSurface{T}(radius, diameter, diameter)

mechanical_diameter(s::StandardSurface) = s.mechanical_diameter

edge_thickness(::StandardSurface, sd::AbstractSphericalSurfaceSDF) = abs(sag(sd))

function sdf(s::StandardSurface, ot::AbstractOrientationType)
    isinf(radius(s)) && return nothing

    return _sdf(s, ot)
end

function _sdf(s::StandardSurface, ::ForwardOrientation)
    front = if radius(s) > 0
        ConvexSphericalSurfaceSDF(radius(s), diameter(s))
    else
        ConcaveSphericalSurfaceSDF(abs(radius(s)), diameter(s))
    end

    return front
end

function _sdf(s::StandardSurface, ::BackwardOrientation)
    back = if radius(s) > 0
        ConcaveSphericalSurfaceSDF(radius(s), diameter(s))
    else
        ConvexSphericalSurfaceSDF(abs(radius(s)), diameter(s))
    end
    zrotate3d!(back, π)

    return back
end

sdf(s::StandardSurface, ::ForwardLeftMeniscusOrientation) = ConvexSphericalSurfaceSDF(radius(s), diameter(s))
sdf(s::StandardSurface, ::BackwardLeftMeniscusOrientation) = SphereSDF(radius(s))

sdf(s::StandardSurface, ::ForwardRightMeniscusOrientation) = SphereSDF(abs(radius(s)))
sdf(s::StandardSurface, ::BackwardRightMeniscusOrientation) = ConvexSphericalSurfaceSDF(abs(radius(s)), diameter(s))
