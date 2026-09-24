#=
CAD-like look of the rendered objects: material presets per component class, feature edge lines
and the studio lighting rig. The component classes are assigned in `RenderPresets.jl`.
=#

"""Plot attributes of a material preset, see `_MATERIALS`."""
const _Material = @NamedTuple{color::RGBf, alpha::Float32, transparency::Bool, diffuse::Float32,
    specular::Float32, shininess::Float32}

_preset(color, alpha, transparency, diffuse, specular, shininess) =
    _Material((color, alpha, transparency, diffuse, specular, shininess))

"""
    _MATERIALS

Material presets per component class, see `_material_class`. The attributes are passed to the mesh
plot of an object, explicit kwargs of `render!` (e.g. `color`) override them.

| class         | components                                |
|:--------------|:------------------------------------------|
| `:refractive` | lenses, prisms, plates, windows           |
| `:reflective` | mirrors, retroreflector                   |
| `:coating`    | beamsplitter coatings                     |
| `:polarizer`  | polarization filters                      |
| `:detector`   | detectors                                 |
| `:mechanics`  | mechanics, dummies and other objects      |
| `:interface`  | cemented interfaces of doublets, triplets |
"""
const _MATERIALS = Dict{Symbol, _Material}(
    :refractive => _preset(RGBf(0.72, 0.85, 0.92), 0.35, true, 0.6, 0.8, 96),
    :reflective => _preset(RGBf(0.82, 0.83, 0.85), 1, false, 0.5, 1.0, 128),
    :coating => _preset(RGBf(0.85, 0.55, 0.85), 0.5, true, 0.6, 0.6, 64),
    :polarizer => _preset(RGBf(0.15, 0.30, 0.35), 0.8, true, 0.7, 0.4, 32),
    :detector => _preset(RGBf(0.12, 0.22, 0.35), 1, false, 0.8, 0.2, 16),
    :mechanics => _preset(RGBf(0.55, 0.56, 0.58), 1, false, 0.9, 0.15, 16),
    :interface => _preset(RGBf(0.95, 0.80, 0.45), 0.25, true, 0.6, 0.3, 32),
)

"""
    _material_class(obj)

Returns the material class of the object `obj`, i.e. a key of `_MATERIALS`. Defaults to
`:mechanics`, the component classes are assigned in `RenderPresets.jl`.
"""
_material_class(::Any) = :mechanics

"""
    _material(obj, material = nothing)

Returns the plot attributes of the `material` (a key of `_MATERIALS`), or of the material class of
`obj` if `material` is `nothing`.
"""
function _material(obj, material = nothing)
    class = isnothing(material) ? _material_class(obj) : material
    if !(class isa Symbol && haskey(_MATERIALS, class))
        throw(ArgumentError("unknown material $(repr(class)), must be one of $(sort!(collect(keys(_MATERIALS))))"))
    end
    return _MATERIALS[class]
end

#=
Feature edges
=#

"""Minimum angle between the normals of adjacent faces of a feature edge, see `_feature_edges`."""
const _EDGE_ANGLE = deg2rad(30)

"""Color of the feature edge lines."""
const _EDGE_COLOR = RGBAf(0.1, 0.1, 0.12, 0.8)

"""Plot attributes of the object that are passed on to its feature edge lines."""
const _EDGE_KWARGS = (:visible, :clip_planes)

"""
    _weld(points, tol)

Welds the `points` by position within about `tol`, returns the representative index of each point.
"""
function _weld(points, tol)
    cell(p) = (floor(Int, p[1] / tol), floor(Int, p[2] / tol), floor(Int, p[3] / tol))
    grid = Dict{NTuple{3, Int}, Int}()
    sizehint!(grid, length(points))
    ids = zeros(Int, length(points))
    for (i, p) in enumerate(points)
        c = cell(p)
        rep = get(grid, c, 0)
        if rep == 0
            # points within tol may lie in a neighboring cell
            for d in Iterators.product(-1:1, -1:1, -1:1)
                j = get(grid, c .+ d, 0)
                if j != 0 && norm(points[j] - p) ≤ tol
                    rep = j
                    break
                end
            end
            rep == 0 && (rep = i)
            grid[c] = rep
        end
        ids[i] = rep
    end
    return ids
end

"""
    _feature_edges(m::_TriMesh; angle = 30°)

Returns the feature edges of `m` as pairs of points, see `_polylines`. The vertices are welded by
position first (tolerance `1e-9` of the mesh size). An edge is a feature edge if the normals of its
two adjacent faces differ by more than `angle`, or if it is a boundary edge (one adjacent face).
Non-manifold edges (more than two faces) are feature edges if not all faces are parallel, i.e.
the edges of coincident flat faces are skipped. The orientation of the faces of open
(`two_sided`) meshes is ignored.
"""
function _feature_edges(m::_TriMesh; angle = _EDGE_ANGLE)
    pts = Point3f[]
    isempty(m.faces) && return pts
    lo, hi = _bbox(m)
    tol = 1e-9 * max(maximum(hi - lo), floatmin(Float64))
    ids = _weld(m.points, tol)
    normals = Vector{Vec3d}(undef, length(m.faces))
    cmin = cos(angle)
    # edge => (number of faces, first face, second face, any face not parallel to the first face)
    edges = Dict{Tuple{Int, Int}, Tuple{Int, Int, Int, Bool}}()
    sizehint!(edges, 2 * length(m.faces))
    for (k, f) in enumerate(m.faces)
        a, b, c = ids[f[1]], ids[f[2]], ids[f[3]]
        (a == b || b == c || a == c) && continue
        g = cross(m.points[b] - m.points[a], m.points[c] - m.points[a])
        normals[k] = norm(g) > 0 ? g / norm(g) : Vec3d(0)
        for (i, j) in ((a, b), (b, c), (c, a))
            key = minmax(i, j)
            n, f1, f2, bent = get(edges, key, (0, 0, 0, false))
            if n == 0
                edges[key] = (1, k, 0, false)
            else
                bent |= abs(dot(normals[f1], normals[k])) < cmin
                edges[key] = (n + 1, f1, n == 1 ? k : f2, bent)
            end
        end
    end
    for ((i, j), (n, f1, f2, bent)) in edges
        if n == 2
            d = dot(normals[f1], normals[f2])
            (m.two_sided ? abs(d) : d) ≥ cmin && continue
        elseif n > 2
            # e.g. the coincident faces of touching parts, which are not removed within `_TAU`
            bent || continue
        end
        push!(pts, Point3f(m.points[i]), Point3f(m.points[j]))
    end
    return pts
end

"""
    _polylines(segments)

Chains the `segments` (pairs of points, see `_feature_edges`) into polylines, which are separated
by `NaN` points. Segments are joined at points that are shared by exactly two segments, e.g. the
segments of a circle are joined into a closed loop.
"""
function _polylines(segments::Vector{Point3f})
    ids = Dict{Point3f, Int}()
    vertex(p) = get!(ids, p, length(ids) + 1)
    edges = [(vertex(segments[k]), vertex(segments[k + 1])) for k in 1:2:length(segments)]
    points = Vector{Point3f}(undef, length(ids))
    for (p, i) in ids
        points[i] = p
    end
    adjacent = [Int[] for _ in points]
    for (e, (a, b)) in enumerate(edges)
        push!(adjacent[a], e)
        push!(adjacent[b], e)
    end
    used = falses(length(edges))
    # follows the chain from vertex `v` via edge `e`, appends the vertices to `chain`
    function walk!(chain, v, e)
        while true
            a, b = edges[e]
            v = a == v ? b : a
            push!(chain, v)
            length(adjacent[v]) == 2 || return chain
            e = adjacent[v][1] == e ? adjacent[v][2] : adjacent[v][1]
            used[e] && return chain
            used[e] = true
        end
    end
    out = Point3f[]
    for e in eachindex(edges)
        used[e] && continue
        used[e] = true
        a, b = edges[e]
        forward = walk!([a], a, e)
        backward = forward[end] == a ? Int[] : walk!(Int[], b, e)
        isempty(out) || push!(out, Point3f(NaN))
        append!(out, points[i] for i in Iterators.flatten((reverse(backward)[1:(end - 1)], forward)))
    end
    return out
end

"""
    _plot_edges!(ax, meshes; kwargs...)

Plots the feature edges of all `meshes` (see `_feature_edges`) as a single `lines!` plot, see
`_polylines`. Of the `kwargs`, i.e. the plot attributes of the object, only `_EDGE_KWARGS` and
`transparency` are used.
"""
function _plot_edges!(ax::_RenderEnv, meshes; transparency::Bool = false, kwargs...)
    pts = Point3f[]
    for m in meshes
        append!(pts, _feature_edges(m))
    end
    isempty(pts) && return nothing
    kw = (; (k => v for (k, v) in pairs(kwargs) if k in _EDGE_KWARGS)...)
    # Polylines instead of line segments, since short segments of thin lines look frayed (GLMakie).
    # The depth shift keeps the lines in front of the faces they bound. The lines of transparent
    # objects are transparent as well, otherwise the faces along the silhouette cover them
    # partially (GLMakie, order independent transparency)
    lines!(ax, _polylines(pts); color = _EDGE_COLOR, linewidth = 1, transparency,
        inspectable = false, depth_shift = -1.0f-5, kw...)
    return nothing
end

#=
Lighting
=#

# Directions along which the light travels, relative to the camera: x right, y up, z towards the viewer
const _KEY_DIRECTION = Vec3f(-0.46, -0.63, -0.63)   # from the upper right front
const _FILL_DIRECTION = Vec3f(1.0, -0.2, -0.5)      # from the left
const _RIM_DIRECTION = Vec3f(0.0, -0.4, 1.0)        # from behind

"""
    _studio_lights(multi::Bool)

Returns the ambient light color and the directional lights of the `:studio` rig. With `multi = true`
(backends with `MultiLightShading`, i.e. GLMakie) a key, fill and rim light, otherwise the key light
only.
"""
function _studio_lights(multi::Bool)
    light(c, dir) = Makie.DirectionalLight(RGBf(c, c, c), dir, true)
    multi || return RGBf(0.45, 0.45, 0.45), [light(0.75, _KEY_DIRECTION)]
    return RGBf(0.35, 0.35, 0.35),
        [light(0.8, _KEY_DIRECTION), light(0.35, _FILL_DIRECTION), light(0.3, _RIM_DIRECTION)]
end

"""Returns `true` if the active Makie backend supports several lights, i.e. `MultiLightShading`."""
function _multi_light_backend()
    backend = Makie.current_backend()
    return !ismissing(backend) && nameof(backend) === :GLMakie
end

"""
    _apply_lighting!(scene, preset::Symbol, multi::Bool)

Sets the lights of the `preset` (`:studio` or `:none`) in the `scene`, see `_studio_lights`.
"""
function _apply_lighting!(scene, preset::Symbol, multi::Bool)
    preset in (:studio, :none) ||
        throw(ArgumentError("unknown lighting preset :$preset, must be :studio or :none"))
    preset === :none && return nothing
    ambient, lights = _studio_lights(multi)
    Makie.set_ambient_light!(scene, ambient)
    Makie.set_lights!(scene, lights)
    return nothing
end

"""
    studio_lighting!(ls::LScene; preset = :studio)

Sets up the lighting of the 3D view `ls`. The `:studio` rig consists of an ambient light, a key
light from the upper right front, a fill light from the left and a rim light from behind, all
relative to the camera. Backends that support a single directional light only (e.g. CairoMakie)
get the ambient and the key light. `preset = :none` leaves the lights unchanged.

[`live_view`](@ref) applies the rig by default. Call it for scenes created via `render!`, e.g.

```julia
fig = Figure()
ax = LScene(fig[1, 1])
studio_lighting!(ax)
render!(ax, system)
```
"""
function studio_lighting!(ls::LScene; preset::Symbol = :studio)
    _apply_lighting!(ls.scene, preset, _multi_light_backend())
    return nothing
end
