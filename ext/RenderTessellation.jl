#=
Analytic tessellation of SDF shapes into triangle meshes.

Every supported shape is tessellated in its local coordinates and transformed into world
coordinates via its `orientation` and `position`. Composite shapes are assembled from the meshes
of their parts, see `_union` and `_difference`. The resulting meshes are rendered as one
`Makie.Mesh` plot per object, see `_render_mesh!`.
=#

"""
    _TriMesh

Triangle mesh as produced by [`_tessellate`](@ref), usually in world coordinates.

# Fields

- `points`: vertex positions
- `normals`: unit vertex normals, duplicated vertices along sharp edges
- `faces`: vertex index triples, counter-clockwise seen from the outside
- `tags`: per vertex, `0` for the main color and `1` for the substrate color, see `_vertex_colors`
- `two_sided`: lit from both sides, used for open (zero thickness) meshes
"""
struct _TriMesh
    points::Vector{Point3d}
    normals::Vector{Vec3d}
    faces::Vector{NTuple{3, Int}}
    tags::Vector{UInt8}
    two_sided::Bool
end

_TriMesh(; two_sided::Bool = false) = _TriMesh(Point3d[], Vec3d[], NTuple{3, Int}[], UInt8[], two_sided)

"""Number of segments of a full revolution, i.e. an angular step of 3°."""
const _N_THETA = 120

"""Number of radial rings of curved lens surfaces."""
const _N_RADIAL = 16

"""Maximum angular step of curved surfaces in [rad]."""
const _MAX_STEP = deg2rad(3)

"""
Tolerance of the composition relative to the size of the shape, see `_clip`. Must resolve thin
parts like the 80 nm cap of a lens surface with a radius of 1000 m.
"""
const _TAU = 1e-9

"""
    _has_mesh(x)

Returns `true` if [`_tessellate`](@ref) supports `x`, i.e. `x` is rendered from an analytic mesh
instead of the marching cubes fallback.
"""
_has_mesh(::Any) = false

"""
    _tessellate(s)

Returns the [`_TriMesh`](@ref) of the shape `s` in world coordinates. See `_has_mesh`.
"""
function _tessellate end

"""
    _mesh(s)

Returns the `GeometryBasics.Mesh` of the shape `s` in world coordinates, with per-vertex normals.
"""
_mesh(s) = _to_geometry(_tessellate(s))

function _to_geometry(m::_TriMesh)
    faces = [GLTriangleFace(f...) for f in m.faces]
    return Mesh(m.points, faces; normal = [Vec3f(n...) for n in m.normals])
end

function _add_vertex!(m::_TriMesh, p, n, tag::UInt8 = 0x00)
    push!(m.points, Point3d(p...))
    push!(m.normals, normalize(Vec3d(n...)))
    push!(m.tags, tag)
    return length(m.points)
end

"""
    _add_triangle!(m, a, b, c)

Adds the triangle `a`, `b`, `c` to `m`, oriented such that its geometric normal agrees with the
vertex normals. Degenerate triangles are skipped.
"""
function _add_triangle!(m::_TriMesh, a::Int, b::Int, c::Int)
    pa, pb, pc = m.points[a], m.points[b], m.points[c]
    g = cross(pb - pa, pc - pa)
    l2 = max(sum(abs2, pb - pa), sum(abs2, pc - pa), sum(abs2, pc - pb))
    norm(g) ≤ 1e-12 * l2 && return nothing
    if dot(g, m.normals[a] + m.normals[b] + m.normals[c]) < 0
        b, c = c, b
    end
    push!(m.faces, (a, b, c))
    return nothing
end

"""
    _add_rows!(m, rows)

Triangulates consecutive `rows` of vertex indices around a closed loop. A row with a single index
is a pole (e.g. the center of a disc).
"""
function _add_rows!(m::_TriMesh, rows::Vector{Vector{Int}})
    for k in 1:(length(rows) - 1)
        A, B = rows[k], rows[k + 1]
        n = max(length(A), length(B))
        n == 1 && continue
        for j in 1:n
            j2 = mod1(j + 1, n)
            if length(A) == 1
                _add_triangle!(m, A[1], B[j], B[j2])
            elseif length(B) == 1
                _add_triangle!(m, A[j], A[j2], B[1])
            else
                _add_triangle!(m, A[j], A[j2], B[j2])
                _add_triangle!(m, A[j], B[j2], B[j])
            end
        end
    end
    return m
end

"""
    _revolve!(m, segment; tag = 0x00)

Revolves the profile `segment` around the local y-axis and adds it to `m`. The `segment` is a
vector of `(r, y, nr, ny)`, i.e. radius, height and the outward normal in the `r`-`y`-plane.
Points with `r = 0` become poles. Sharp edges require separate segments.
"""
function _revolve!(m::_TriMesh, segment; tag::UInt8 = 0x00)
    θs = [2π * (j - 1) / _N_THETA for j in 1:_N_THETA]
    rows = Vector{Int}[]
    for (r, y, nr, ny) in segment
        if r ≤ 0
            push!(rows, [_add_vertex!(m, (0, y, 0), (0, sign(ny), 0), tag)])
        else
            push!(rows, [_add_vertex!(m, (r * cos(θ), y, r * sin(θ)), (nr * cos(θ), ny, nr * sin(θ)), tag)
                         for θ in θs])
        end
    end
    return _add_rows!(m, rows)
end

# Flat profile segments of revolved shapes, `s = ±1` is the sign of the outward normal
_disc(y, r, s) = [(0.0, y, 0.0, s), (r, y, 0.0, s)]
_annulus(y, r1, r2, s) = [(r1, y, 0.0, s), (r2, y, 0.0, s)]
_wall(r, y1, y2, s) = [(r, y1, s, 0.0), (r, y2, s, 0.0)]

"""
    _extrude!(m, ws, lo, hi, nlo, nhi, h, axes)

Extrudes the 2D region `lo(w) ≤ v ≤ hi(w)`, sampled at `ws`, between `t = ±h` and adds it to `m`.
`nlo` and `nhi` are the outward 2D normals `(nw, nv)` of the lower and upper boundary. Walls are
added at the first and last sample if `hi > lo` there. `axes = (iw, iv, it)` are the local
coordinate indices of `w`, `v` and `t`.
"""
function _extrude!(m::_TriMesh, ws, lo, hi, nlo, nhi, h, axes)
    iw, iv, it = axes
    function vec3(w, v, t)
        x = zeros(3)
        x[iw], x[iv], x[it] = w, v, t
        return x
    end
    function band!(curve)
        rows = [[_add_vertex!(m, vec3(w, v, -h), vec3(nw, nv, 0)),
                    _add_vertex!(m, vec3(w, v, h), vec3(nw, nv, 0))] for (w, v, nw, nv) in curve]
        for k in 1:(length(rows) - 1)
            a, b = rows[k]
            c, d = rows[k + 1]
            _add_triangle!(m, a, c, d)
            _add_triangle!(m, a, d, b)
        end
    end
    n = length(ws)
    band!([(ws[i], lo[i], nlo[i]...) for i in 1:n])
    band!([(ws[i], hi[i], nhi[i]...) for i in 1:n])
    tol = 1e-12 * max(maximum(abs, ws), maximum(abs, hi), maximum(abs, lo))
    hi[1] - lo[1] > tol && band!([(ws[1], lo[1], -1, 0), (ws[1], hi[1], -1, 0)])
    hi[n] - lo[n] > tol && band!([(ws[n], lo[n], 1, 0), (ws[n], hi[n], 1, 0)])
    # caps
    for s in (-1, 1)
        L = [_add_vertex!(m, vec3(ws[i], lo[i], s * h), vec3(0, 0, s)) for i in 1:n]
        H = [_add_vertex!(m, vec3(ws[i], hi[i], s * h), vec3(0, 0, s)) for i in 1:n]
        for i in 1:(n - 1)
            _add_triangle!(m, L[i], L[i + 1], H[i + 1])
            _add_triangle!(m, L[i], H[i + 1], H[i])
        end
    end
    return m
end

"""Transforms the local mesh `m` into world coordinates, i.e. `p ↦ R * p + P`, `n ↦ R * n`."""
function _transform!(m::_TriMesh, R, P)
    Rm = SMatrix{3, 3, Float64}(R)
    Pv = Vec3d(P...)
    for i in eachindex(m.points)
        m.points[i] = Point3d(Rm * m.points[i] + Pv)
        m.normals[i] = normalize(Vec3d(Rm * m.normals[i]))
    end
    return m
end

_transform!(m::_TriMesh, s::BMO.AbstractSDF) = _transform!(m, BMO.orientation(s), BMO.position(s))

"""Concatenates the `meshes` into a single [`_TriMesh`](@ref)."""
function _merge(meshes)
    out = _TriMesh(; two_sided = any(m -> m.two_sided, meshes))
    nv = sum(m -> length(m.points), meshes; init = 0)
    sizehint!(out.points, nv)
    sizehint!(out.normals, nv)
    sizehint!(out.tags, nv)
    sizehint!(out.faces, sum(m -> length(m.faces), meshes; init = 0))
    for m in meshes
        offset = length(out.points)
        append!(out.points, m.points)
        append!(out.normals, m.normals)
        append!(out.tags, m.tags)
        append!(out.faces, (f .+ offset for f in m.faces))
    end
    return out
end

"""Returns the lower and upper corner of the axis-aligned bounding box of the `meshes`."""
function _bbox(meshes...)
    lo = Vec3d(Inf)
    hi = Vec3d(-Inf)
    for m in meshes, p in m.points
        lo = min.(lo, p)
        hi = max.(hi, p)
    end
    all(isfinite, lo) || return Vec3d(0), Vec3d(0)
    return lo, hi
end

"""Flips the orientation of `m`, i.e. reverses the faces and negates the normals."""
function _flip!(m::_TriMesh)
    m.faces .= [(f[1], f[3], f[2]) for f in m.faces]
    m.normals .= .-m.normals
    return m
end

#=
Primitive shapes
=#

function _tessellate(s::BMO.BoxSDF)
    hx, hy, hz = s.dimensions
    m = _TriMesh()
    _extrude!(m, [-hx, hx], [-hy, -hy], [hy, hy], [(0, -1), (0, -1)], [(0, 1), (0, 1)], hz, (1, 2, 3))
    return _transform!(m, s)
end

# The prism is the box cut by x + y ≤ 0 with symmetric legs, see `RightAnglePrismSDF`
function _tessellate(s::BMO.RightAnglePrismSDF)
    hx, hy, hz = s.dimensions
    m = _TriMesh()
    n = (1 / sqrt(2), 1 / sqrt(2))
    _extrude!(m, [-hx, hx], [-hy, -hy], [hy, -hx], [(0, -1), (0, -1)], [n, n], hz, (1, 2, 3))
    return _transform!(m, s)
end

function _tessellate(s::BMO.CylinderSDF)
    r, h = s.radius, s.height
    m = _TriMesh()
    _revolve!(m, _disc(-h, r, -1))
    _revolve!(m, _wall(r, -h, h, 1))
    _revolve!(m, _disc(h, r, 1))
    return _transform!(m, s)
end

function _tessellate(s::BMO.PlanoSurfaceSDF)
    r, t = BMO.diameter(s) / 2, BMO.thickness(s)
    m = _TriMesh()
    _revolve!(m, _disc(0.0, r, -1))
    _revolve!(m, _wall(r, 0.0, t, 1))
    _revolve!(m, _disc(t, r, 1))
    return _transform!(m, s)
end

function _tessellate(s::BMO.RingSDF)
    ri, ro = s.inner_radius - s.hwidth, s.inner_radius + s.hwidth
    h = s.hthickness
    m = _TriMesh()
    _revolve!(m, _annulus(-h, ri, ro, -1))
    _revolve!(m, _wall(ro, -h, h, 1))
    _revolve!(m, _annulus(h, ri, ro, 1))
    _revolve!(m, _wall(ri, -h, h, -1))
    return _transform!(m, s)
end

function _tessellate(s::BMO.SphereSDF)
    R = s.radius
    n = ceil(Int, π / _MAX_STEP)
    αs = LinRange(0, π, n + 1)
    arc = [(i == 1 || i == n + 1 ? 0.0 : R * sin(α), -R * cos(α), sin(α), -cos(α))
           for (i, α) in enumerate(αs)]
    m = _TriMesh()
    _revolve!(m, arc)
    return _transform!(m, s)
end

# Sphere of `radius` cut by the plane y = `height`, the solid is y ≥ height, see `CutSphereSDF`
function _tessellate(s::BMO.CutSphereSDF)
    R, h, w = s.radius, s.height, s.w
    α_max = acos(h / R)
    n = max(16, ceil(Int, α_max / _MAX_STEP))
    # polar angles from the cut to the pole, including the equator (the extent in x and z)
    αs = sort(unique([collect(LinRange(α_max, 0, n + 1)); α_max > π / 2 ? [π / 2] : Float64[]]); rev = true)
    arc = [(α == 0 ? 0.0 : R * sin(α), R * cos(α), sin(α), cos(α)) for α in αs]
    arc[1] = (w, h, sin(α_max), cos(α_max))
    m = _TriMesh()
    _revolve!(m, _disc(h, w, -1))
    _revolve!(m, arc)
    return _transform!(m, s)
end

const _PrimitiveSDF = Union{BMO.BoxSDF, BMO.RightAnglePrismSDF, BMO.CylinderSDF,
    BMO.PlanoSurfaceSDF, BMO.RingSDF, BMO.SphereSDF, BMO.CutSphereSDF}

_has_mesh(::_PrimitiveSDF) = true

#=
Meshes of BeamletOptics.Mesh shapes
=#

_has_mesh(::BMO.AbstractMesh) = true

"""
    _tessellate(mesh::AbstractMesh)

Flat shaded copy of the `mesh`, i.e. one vertex per face corner with the face normal. The mesh is
lit from both sides, since meshes like the `Retroreflector` or the `Detector` are open.
"""
function _tessellate(mesh::BMO.AbstractMesh)
    V, F = BMO.vertices(mesh), BMO.faces(mesh)
    m = _TriMesh(; two_sided = true)
    vertex(i) = Point3d(V[i, 1], V[i, 2], V[i, 3])
    for f in axes(F, 1)
        p1, p2, p3 = vertex(F[f, 1]), vertex(F[f, 2]), vertex(F[f, 3])
        g = cross(p2 - p1, p3 - p1)
        norm(g) > 0 || continue
        n = normalize(g)
        push!(m.faces, (_add_vertex!(m, p1, n), _add_vertex!(m, p2, n), _add_vertex!(m, p3, n)))
    end
    return m
end

#=
Composition
=#

"""
    _crossing(g, pa, pb, fa, fb)

Returns the point on the segment `pa`-`pb` where `g` changes its sign, with `fa = g(pa)` and
`fb = g(pb)` of opposite sign (Illinois variant of the regula falsi).
"""
function _crossing(g, pa, pb, fa, fb)
    a, b = 0.0, 1.0
    t = 0.5
    for _ in 1:8
        t = (a * fb - b * fa) / (fb - fa)
        ft = g(pa + t * (pb - pa))
        abs(ft) < 1e-15 && break
        if sign(ft) == sign(fa)
            a, fa = t, ft
            fb /= 2
        else
            b, fb = t, ft
            fa /= 2
        end
    end
    return pa + t * (pb - pa)
end

"""
    _clip(m, g, own, τ; on = nothing, interface = nothing)

Returns the part of `m` where `g > τ`, e.g. outside of the other parts of a union. Triangles with
all vertices on the boundary (`|g| ≤ τ`, i.e. coincident faces) are kept if `g > τ` holds at their
centroid, which is projected onto the surface of `own` (the SDF of `m`, or `nothing`) first, since
the centroids of curved triangles lie inside the true surface. Triangles with vertices on both
sides are cut along `g = 0`.

Removed coincident triangles whose projected centroid also lies on the surface of `on`
(`|on| ≤ τ`) are added to the [`_TriMesh`](@ref) `interface`, e.g. the cemented faces of a lens.
"""
function _clip(m::_TriMesh, g, own, τ; on = nothing, interface::Union{Nothing, _TriMesh} = nothing)
    n = length(m.points)
    f = [g(p) for p in m.points]
    cls = [v > τ ? 1 : (v < -τ ? -1 : 0) for v in f]
    # vertices on the cuts get the indices n + 1, n + 2, …
    cut_points, cut_normals, cut_tags = Point3d[], Vec3d[], UInt8[]
    cuts = Dict{Tuple{Int, Int}, Int}()
    function cut(i, j)
        key = i < j ? (i, j) : (j, i)
        get!(cuts, key) do
            a, b = key
            p = _crossing(g, m.points[a], m.points[b], f[a], f[b])
            t = norm(p - m.points[a]) / norm(m.points[b] - m.points[a])
            push!(cut_points, p)
            push!(cut_normals, normalize((1 - t) * m.normals[a] + t * m.normals[b]))
            push!(cut_tags, m.tags[a])
            n + length(cut_points)
        end
    end
    kept = NTuple{3, Int}[]
    coincident = NTuple{3, Int}[]
    poly = Int[]
    for tri in m.faces
        c = map(i -> cls[i], tri)
        if all(≥(0), c) && any(>(0), c)
            push!(kept, tri)
        elseif all(≤(0), c) && any(<(0), c)
            continue
        elseif all(==(0), c)
            ctr = Point3d(sum(i -> Vec3d(m.points[i]), tri) / 3)
            if !isnothing(own)
                # each step reduces the distance by about α²/2, α the angle between n and the surface normal
                nc = normalize(sum(i -> m.normals[i], tri))
                for _ in 1:3
                    d = _sdf(own, ctr)
                    abs(d) ≤ τ && break
                    ctr = ctr - d * nc
                end
            end
            if g(ctr) > τ
                push!(kept, tri)
            elseif !isnothing(on) && abs(on(ctr)) ≤ τ
                push!(coincident, tri)
            end
        else
            # Sutherland-Hodgman clipping of the triangle against g > 0
            empty!(poly)
            for k in 1:3
                i, j = tri[k], tri[mod1(k + 1, 3)]
                f[i] > 0 && push!(poly, i)
                (f[i] > 0) != (f[j] > 0) && push!(poly, cut(i, j))
            end
            for k in 2:(length(poly) - 1)
                push!(kept, (poly[1], poly[k], poly[k + 1]))
            end
        end
    end
    # copy the referenced vertices only
    out = _TriMesh(; two_sided = m.two_sided)
    ids = zeros(Int, n + length(cut_points))
    function vertex(k)
        if ids[k] == 0
            push!(out.points, k ≤ n ? m.points[k] : cut_points[k - n])
            push!(out.normals, k ≤ n ? m.normals[k] : cut_normals[k - n])
            push!(out.tags, k ≤ n ? m.tags[k] : cut_tags[k - n])
            ids[k] = length(out.points)
        end
        return ids[k]
    end
    sizehint!(out.faces, length(kept))
    for tri in kept
        push!(out.faces, map(vertex, tri))
    end
    if !isnothing(interface) && !isempty(coincident)
        # the coincident triangles are not cut, i.e. consist of vertices of m only
        used = unique(Iterators.flatten(coincident))
        offset = length(interface.points)
        new_ids = Dict(k => offset + i for (i, k) in enumerate(used))
        append!(interface.points, m.points[used])
        append!(interface.normals, m.normals[used])
        append!(interface.tags, m.tags[used])
        append!(interface.faces, (map(k -> new_ids[k], tri) for tri in coincident))
    end
    return out
end

"""
    _refine(m, L)

Splits all edges of `m` longer than `L` at their midpoint, until no such edge remains. Shared
edges are split once, hence closed meshes stay closed.
"""
function _refine(m::_TriMesh, L)
    m = _TriMesh(copy(m.points), copy(m.normals), copy(m.faces), copy(m.tags), m.two_sided)
    while true
        mids = Dict{Tuple{Int, Int}, Int}()
        for tri in m.faces, k in 1:3
            i, j = tri[k], tri[mod1(k + 1, 3)]
            key = i < j ? (i, j) : (j, i)
            haskey(mids, key) && continue
            norm(m.points[i] - m.points[j]) > L || continue
            mids[key] = _add_vertex!(m, (m.points[i] + m.points[j]) / 2, m.normals[i] + m.normals[j], m.tags[i])
        end
        isempty(mids) && return m
        faces = NTuple{3, Int}[]
        for tri in m.faces
            mid(k) = get(mids, minmax(tri[k], tri[mod1(k + 1, 3)]), 0)
            split = [mid(k) for k in 1:3]
            ns = count(>(0), split)
            if ns == 0
                push!(faces, tri)
            elseif ns == 3
                a, b, c = tri
                ab, bc, ca = split
                append!(faces, [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)])
            else
                # rotate such that edge 1 (a-b) is split and, for two splits, edge 3 (c-a) is not
                r = ns == 1 ? findfirst(>(0), split) : findfirst(==(0), split) % 3 + 1
                a, b, c = tri[r], tri[mod1(r + 1, 3)], tri[mod1(r + 2, 3)]
                ab, bc = split[r], split[mod1(r + 1, 3)]
                if ns == 1
                    append!(faces, [(a, ab, c), (ab, b, c)])
                else
                    append!(faces, [(ab, b, bc), (a, ab, bc), (a, bc, c)])
                end
            end
        end
        empty!(m.faces)
        append!(m.faces, faces)
    end
end

"""
    _BoundedSDF(s, lo, hi, pad)

Callable `p -> sdf(s, p)` that skips the SDF evaluation outside of the box `lo`, `hi` padded by
`pad` and returns the (smaller, positive) distance to that box instead.
"""
struct _BoundedSDF{S}
    s::S
    lo::Vec3d
    hi::Vec3d
    pad::Float64
end

function (b::_BoundedSDF)(p)
    d = max.(b.lo .- b.pad .- p, p .- b.hi .- b.pad, 0.0)
    any(>(0), d) && return norm(d) + b.pad
    return _sdf(b.s, p)
end

# Type stable evaluation of `sdf(::UnionSDF)`, whose generator over the parts allocates
_sdf(s, p) = BMO.sdf(s, p)
_sdf(u::BMO.UnionSDF, p) = _min_sdf(u.sdfs, p)

# unrolled over the heterogeneous tuple of parts
_min_sdf(t::Tuple, p) = min(_sdf(first(t), p), _min_sdf(Base.tail(t), p))
_min_sdf(::Tuple{}, p) = Inf

"""Callable `p -> minimum(f(p) for f in fs)` for a tuple `fs`, type stable for the composition."""
struct _MinOf{F <: Tuple}
    fs::F
end

(m::_MinOf)(p) = _min_call(m.fs, p)

_min_call(t::Tuple, p) = min(first(t)(p), _min_call(Base.tail(t), p))
_min_call(::Tuple{}, p) = Inf

"""Callable `p -> -sdf(s, p)`, i.e. the SDF of the complement of `s`."""
struct _Outside{S}
    s::S
end

(o::_Outside)(p) = -BMO.sdf(o.s, p)

"""
    _union(parts; interface = nothing, cemented = nothing)

Returns the merged mesh of the union of `parts`, a vector of `(mesh, sdf)`. Triangles of a part
that lie inside or on the surface of another part are removed, see `_clip`. Parts with
`sdf = nothing` are neither filtered nor used for filtering.

If a [`_TriMesh`](@ref) `interface` is given, the faces where two `cemented` parts (one `Bool` per
part) touch are added to it once, i.e. the copy of the first of both parts.
"""
function _union(parts; interface::Union{Nothing, _TriMesh} = nothing, cemented = nothing)
    length(parts) == 1 && return first(first(parts))
    lo, hi = _bbox((m for (m, _) in parts)...)
    τ = _TAU * maximum(hi - lo)
    # the padding exceeds the chord error of the curved faces
    pad = 1e-3 * maximum(hi - lo)
    bounded = [isnothing(s) ? nothing : _BoundedSDF(s, _bbox(m)..., pad) for (m, s) in parts]
    meshes = _TriMesh[]
    for (i, (m, s)) in enumerate(parts)
        others = Tuple(bounded[j] for j in eachindex(parts) if j != i && !isnothing(bounded[j]))
        if isnothing(s) || isempty(others)
            push!(meshes, m)
            continue
        end
        on = nothing
        if !isnothing(interface) && cemented[i]
            later = Tuple(bounded[j] for j in (i + 1):length(parts) if cemented[j] && !isnothing(bounded[j]))
            isempty(later) || (on = _MinOf(later))
        end
        push!(meshes, _clip(m, _MinOf(others), s, τ; on, interface))
    end
    return _merge(meshes)
end

"""
    _difference(base, tools)

Returns the merged mesh of `base \\ (tool_1 ∪ …)`, where `base` and the `tools` are `(mesh, sdf)`.
The meshes are refined first, such that the cut along the tool surfaces follows the geometry.
The kept faces of the tools are flipped.
"""
function _difference(base, tools)
    bm, bs = base
    lo, hi = _bbox(bm)
    size = maximum(hi - lo)
    τ = _TAU * size
    L = size / 24
    sdfs = Tuple(p -> BMO.sdf(ts, p) for (_, ts) in tools)
    meshes = [_clip(_refine(bm, L), _MinOf(sdfs), bs, τ)]
    for (k, (tm, ts)) in enumerate(tools)
        g = _MinOf((_Outside(bs), (sdfs[j] for j in eachindex(sdfs) if j != k)...))
        push!(meshes, _flip!(_clip(_refine(tm, L), g, ts, τ)))
    end
    return _merge(meshes)
end

_has_mesh(u::BMO.UnionSDF) = all(_has_mesh, u.sdfs)

# the parts of a UnionSDF are in world coordinates
_tessellate(u::BMO.UnionSDF) = _union([(_tessellate(s), s) for s in u.sdfs])

_has_mesh(d::BMO.DifferenceSDF) = _has_mesh(d.base) && all(_has_mesh, d.tools)

# dispatch on the base, since a method for `DifferenceSDF{T, <:ConicSDF}` is not more specific
# than one for `DifferenceSDF` (the parameters of the latter are constrained by T)
_tessellate(d::BMO.DifferenceSDF) = _tessellate_difference(d.base, d)

function _tessellate_difference(::Any, d::BMO.DifferenceSDF)
    return _difference((_tessellate(d.base), d.base), [(_tessellate(t), t) for t in d.tools])
end

_has_mesh(ml::BMO.MeniscusLensSDF) = _has_mesh(ml.convex) && _has_mesh(ml.concave)

"""
    _tessellate(ml::MeniscusLensSDF)

The parts of `ml` are given in its local frame. The concave face is the spherical cap of the
`SphereSDF` within the lens diameter, such that its rim matches the rim of the cylinder.
"""
function _tessellate(ml::BMO.MeniscusLensSDF)
    convex, cylinder, concave = ml.convex, ml.cylinder, ml.concave
    if concave isa BMO.SphereSDF
        C = BMO.position(concave)
        Rc = concave.radius
        r = BMO.diameter(cylinder) / 2
        # the cap faces the lens, i.e. the cylinder
        s = sign(BMO.position(cylinder)[2] + BMO.thickness(cylinder) / 2 - C[2])
        cap = [(rk, C[2] + s * sqrt(Rc^2 - rk^2), -rk / Rc, -s * sqrt(Rc^2 - rk^2) / Rc)
               for rk in sqrt.(LinRange(0, 1, _N_RADIAL + 1)) .* r]
        cap_mesh = _revolve!(_TriMesh(), cap)
        g_convex(p) = min(BMO.sdf(cylinder, p), BMO.sdf(concave, p))
        g_cylinder(p) = min(BMO.sdf(convex, p), BMO.sdf(concave, p))
        mc, my = _tessellate(convex), _tessellate(cylinder)
        lo, hi = _bbox(mc, my)
        τ = _TAU * maximum(hi - lo)
        m = _merge([_clip(mc, g_convex, convex, τ), _clip(my, g_cylinder, cylinder, τ), cap_mesh])
    else
        m = _difference((_union([(_tessellate(convex), convex), (_tessellate(cylinder), cylinder)]),
                convex + cylinder), [(_tessellate(concave), concave)])
    end
    return _transform!(m, ml)
end

#=
Rendering
=#

"""
    _MESH_COLLECTOR

Collects the meshes of the parts of a `MultiShape` object while it is rendered, such that they can
be merged into one plot per color, see `render!(ax, ::MultiShape, obj)`. `nothing` outside.
"""
const _MESH_COLLECTOR = ScopedValue{Union{Nothing, Vector{Any}}}(nothing)

"""Default color of the analytic mesh of `s`, `nothing` for the Makie default."""
_default_color(::Any) = nothing

"""
    _vertex_colors(m, color)

Returns `color`, or per-vertex colors if `m` has substrate vertices (tag `1`), which are grey.
"""
function _vertex_colors(m::_TriMesh, color)
    all(iszero, m.tags) && return color
    c0 = Makie.to_color(isnothing(color) ? :silver : color)
    c1 = Makie.to_color(:grey)
    return [t == 0 ? c0 : c1 for t in m.tags]
end

"""
    _plot_mesh!(ax, m; color, kwargs...)

Plots the [`_TriMesh`](@ref) `m` as a single `mesh!`. Open meshes are lit from both sides.
"""
function _plot_mesh!(ax::_RenderEnv, m::_TriMesh; color = nothing, kwargs...)
    isempty(m.faces) && return nothing
    c = _vertex_colors(m, color)
    kw = isnothing(c) ? kwargs : (; color = c, kwargs...)
    if m.two_sided && !haskey(kw, :backlight)
        kw = (; backlight = 1.0f0, kw...)
    end
    mesh!(ax, _to_geometry(m); kw...)
    return nothing
end

"""
    _render_mesh!(ax, s; color = _default_color(s), edges = false, cemented = false, kwargs...)

Renders the analytic mesh of the shape `s` and, if `edges`, its feature edges, see `_plot_edges!`.
Inside of a `MultiShape` object the mesh is collected instead, see `_MESH_COLLECTOR`, together with
`edges`, i.e. whether the part contributes feature edges. `cemented` marks the parts whose contact
faces are cemented interfaces, see `_plot_collected!`.
"""
function _render_mesh!(ax::_RenderEnv, s; color = _default_color(s), edges::Bool = false,
        cemented::Bool = false, kwargs...)
    m = _tessellate(s)
    parts = _MESH_COLLECTOR[]
    if isnothing(parts)
        _plot_mesh!(ax, m; color, kwargs...)
        edges && _plot_edges!(ax, (m,); color, kwargs...)
    else
        push!(parts, (m, s isa BMO.AbstractSDF ? s : nothing, (; color, kwargs...), cemented, edges))
    end
    return nothing
end

"""
    _plot_collected!(ax, parts)

Plots the collected `parts` (see `_MESH_COLLECTOR`), one mesh per set of plot attributes. The parts
of a mesh are merged via `_union`, i.e. coincident interior faces are removed. The faces where two
`cemented` parts touch (e.g. the lenses of a doublet) are plotted once as separate interface meshes
with the `:interface` material, which overrides the look attributes of the parts. Finally, the
feature edges of all meshes of the parts with `edges` are plotted as a single plot, see
`_plot_edges!`.
"""
function _plot_collected!(ax::_RenderEnv, parts)
    keys = Any[]
    groups = Vector{Any}[]
    for (m, s, kw, cemented, edges) in parts
        # Parts with and without edges are merged separately
        key = (kw, m.two_sided, edges)
        i = findfirst(k -> isequal(k, key), keys)
        if isnothing(i)
            push!(keys, key)
            push!(groups, Any[])
            i = length(keys)
        end
        push!(groups[i], (m, s, cemented))
    end
    meshes = _TriMesh[]
    interfaces = Any[]
    edge_kw = nothing
    for (key, group) in zip(keys, groups)
        cemented = Bool[c for (_, _, c) in group]
        interface = count(cemented) > 1 ? _TriMesh(; two_sided = true) : nothing
        m = _union([(m, s) for (m, s, _) in group]; interface, cemented)
        _plot_mesh!(ax, m; key[1]...)
        if key[3]
            push!(meshes, m)
            isnothing(edge_kw) && (edge_kw = key[1])
        end
        isnothing(interface) || isempty(interface.faces) || push!(interfaces, (interface, key[1]))
    end
    for (interface, kw) in interfaces
        _plot_mesh!(ax, interface; kw..., _materials()[:interface]...)
    end
    isempty(meshes) || _plot_edges!(ax, meshes; edge_kw...)
    return nothing
end
