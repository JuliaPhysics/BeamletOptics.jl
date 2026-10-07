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

thickness(s::BoxSDF) = 2*s.dimensions[2]

"""
    bounding_sphere_of(box::BoxSDF)

Returns the sphere through the eight corners of the box, centered at the local origin.
"""
bounding_sphere_of(box::BoxSDF{T}) where T = SingleBoundingSphere(box, Point3{T}(0), norm(box.dimensions))

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
    bounding_sphere_of(cylinder::CylinderSDF)

Returns the sphere through both rims of the cylinder, centered at the local origin. Note that the
`height` field stores half of the cylinder height.
"""
function bounding_sphere_of(cylinder::CylinderSDF{T}) where T
    return SingleBoundingSphere(cylinder, Point3{T}(0), sqrt(cylinder.radius^2 + cylinder.height^2))
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
    bounding_sphere_of(cs::CutSphereSDF)

Returns the smallest sphere around the spherical cap `y ≥ height`: the sphere through the rim,
centered in the cut plane, for a cap of at most a hemisphere (`height ≥ 0`), and the full sphere
otherwise.
"""
function bounding_sphere_of(cs::CutSphereSDF{T}) where T
    cs.height > 0 || return SingleBoundingSphere(cs, Point3{T}(0), cs.radius)
    return SingleBoundingSphere(cs, Point3{T}(0, cs.height, 0), cs.w)
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
        Matrix{TT}(I, 3, 3),
        Matrix{TT}(I, 3, 3),
        zeros(TT, 3),
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
    bounding_sphere_of(ring::RingSDF)

Returns the sphere through both outer rims of the ring, centered at the local origin. Note that the
`inner_radius` field stores the mean radius of the ring.
"""
function bounding_sphere_of(ring::RingSDF{T}) where T
    return SingleBoundingSphere(ring, Point3{T}(0), sqrt((ring.inner_radius + ring.hwidth)^2 + ring.hthickness^2))
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

thickness(s::RightAnglePrismSDF) = 2 * s.dimensions[2]

"""
    bounding_sphere_of(prism::RightAnglePrismSDF)

Returns the sphere through the four corners of the hypotenuse face, centered at the local origin,
which is the center of that face.
"""
bounding_sphere_of(prism::RightAnglePrismSDF{T}) where T = SingleBoundingSphere(prism, Point3{T}(0), norm(prism.dimensions))

"""
    PolygonPrismSDF <: AbstractSDF

Implements the exact `SDF` of a right prism, i.e. the extrusion of a convex polygon along the local
z-axis. The local origin is the origin of the coordinates given to the constructor, it is not
shifted.

# Fields

- `vertices`: corners of the cross-section in the local x-y-plane in [m], in counter-clockwise order
- `normals`: outward unit normals of the edges, where edge `i` runs from vertex `i` to vertex `i + 1`
- `height`: extent along the local z-axis in [m], the prism spans `z ∈ [-height/2, height/2]`
"""
mutable struct PolygonPrismSDF{T} <: AbstractSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    vertices::Vector{Point2{T}}
    normals::Vector{Point2{T}}
    height::T
end

"""
    PolygonPrismSDF(vertices, height)

Constructs the prism over the convex polygon given by `vertices`, a vector of 2D points `(x, y)` in
[m], extruded symmetrically by `height` in [m] along the local z-axis. Both vertex orders are
accepted and normalized to counter-clockwise. Throws an `ArgumentError` for fewer than three
vertices, for points that are not 2D, for degenerate (zero area or collinear) corners and for
non-convex or self-intersecting polygons.
"""
function PolygonPrismSDF(vertices::AbstractVector, height::H) where {H <: Real}
    length(vertices) >= 3 || throw(ArgumentError("a polygon needs at least 3 vertices"))
    all(v -> length(v) == 2, vertices) ||
        throw(ArgumentError("the vertices must be 2D points (x, y) in the local x-y-plane"))
    height > 0 || throw(ArgumentError("the prism height must be positive"))
    T = float(mapreduce(v -> promote_type(typeof(v[1]), typeof(v[2])), promote_type, vertices; init = H))
    pts = [Point2{T}(v[1], v[2]) for v in vertices]
    n = length(pts)
    # signed area (shoelace) of the polygon relative to its first vertex, which makes the test
    # invariant under translation (no cancellation for polygons far from the origin); the
    # orientation is normalized to counter-clockwise
    rel = [p - pts[1] for p in pts]
    area = sum(i -> rel[i][1] * rel[mod1(i + 1, n)][2] - rel[mod1(i + 1, n)][1] * rel[i][2], 1:n) / 2
    scale = maximum(p -> maximum(abs, p), rel)
    abs(area) > eps(T) * scale^2 * n || throw(ArgumentError("the polygon has zero area"))
    area < 0 && reverse!(pts)
    # convex: every corner turns left, and the edges turn by 2π in total (no star polygons)
    turning = zero(T)
    for i in 1:n
        e1 = pts[mod1(i + 1, n)] - pts[i]
        e2 = pts[mod1(i + 2, n)] - pts[mod1(i + 1, n)]
        c = e1[1] * e2[2] - e1[2] * e2[1]
        c > sqrt(eps(T)) * norm(e1) * norm(e2) ||
            throw(ArgumentError("the polygon must be strictly convex (no reflex or collinear vertices)"))
        turning += atan(c, dot(e1, e2))
    end
    isapprox(turning, 2π; atol = sqrt(eps(T))) ||
        throw(ArgumentError("the polygon must not be self-intersecting"))
    normals = map(1:n) do i
        e = pts[mod1(i + 1, n)] - pts[i]
        return normalize(Point2{T}(e[2], -e[1]))
    end
    return PolygonPrismSDF{T}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        Point3{T}(0),
        pts,
        normals,
        T(height))
end

function thickness(s::PolygonPrismSDF)
    lo, hi = extrema(v -> v[2], s.vertices)
    return hi - lo
end

"""
    sdf_convex_polygon(q, vertices, normals)

Exact signed distance of the 2D point `q` to the convex polygon with the counter-clockwise
`vertices` and the outward unit `normals` of its edges, see [`PolygonPrismSDF`](@ref).
"""
function sdf_convex_polygon(q, vertices, normals)
    # the largest signed distance to an edge line is exact inside and wherever the closest point
    # lies within that edge, which keeps the gradient smooth on the faces
    d, k = dot(normals[1], q - vertices[1]), 1
    @inbounds for i in 2:length(vertices)
        di = dot(normals[i], q - vertices[i])
        if di > d
            d, k = di, i
        end
    end
    d > 0 || return d
    # outside and beside that edge, the closest point is the vertex at the end the point lies beyond
    a, b = vertices[k], vertices[k == length(vertices) ? 1 : k + 1]
    t = dot(q - a, b - a)
    t < 0 && return norm(q - a)
    t > sum(abs2, b - a) && return norm(q - b)
    return d
end

function sdf(prism::PolygonPrismSDF, point)
    p = Point3(_world_to_sdf(prism, point))
    return op_extrude_z(p, q -> sdf_convex_polygon(q, prism.vertices, prism.normals), prism.height / 2)
end

"""
    _enclosing_circle(vertices)

Returns the center and the radius of a circle that encloses the 2D points `vertices` and passes
through at least one of them. Of the circles around the mean of the points (exact for regular
polygons) and around the center of their bounding box (exact for rectangles and right triangles)
the smaller one is returned, which is not the smallest enclosing circle in general.
"""
function _enclosing_circle(vertices::AbstractVector{Point2{T}}) where T
    lo = hi = acc = first(vertices)
    @inbounds for i in 2:length(vertices)
        v = vertices[i]
        lo, hi, acc = min.(lo, v), max.(hi, v), acc + v
    end
    c1, c2 = acc / length(vertices), (lo + hi) / 2
    q1 = q2 = zero(T)
    for v in vertices
        q1, q2 = max(q1, sum(abs2, v - c1)), max(q2, sum(abs2, v - c2))
    end
    return q1 ≤ q2 ? (c1, sqrt(q1)) : (c2, sqrt(q2))
end

"""
    bounding_sphere_of(prism::PolygonPrismSDF)

Returns the sphere around the prism over a circle that encloses the cross-section (the smaller one
of the circles around the mean of the vertices and around the center of their bounding box),
centered in the local plane `z = 0`. It passes through at least two corners of the prism.
"""
function bounding_sphere_of(prism::PolygonPrismSDF{T}) where T
    c, r = _enclosing_circle(prism.vertices)
    return SingleBoundingSphere(prism, Point3{T}(c[1], c[2], 0), sqrt(r^2 + (prism.height / 2)^2))
end
