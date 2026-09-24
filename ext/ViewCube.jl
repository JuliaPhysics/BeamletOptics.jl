using Makie: Scene, Rect2i, Rect3f

# Names and outward normals of the faces of the view cube, see `view_cube!`
const _CUBE_FACES = (
    ("Top", (0, 0, 1)), ("Bottom", (0, 0, -1)), ("Front", (0, -1, 0)),
    ("Back", (0, 1, 0)), ("Right", (1, 0, 0)), ("Left", (-1, 0, 0))
)

# Distance of the eye of the cube camera from the center of the cube
const _CUBE_EYE_DIST = 4.0
# Half-width of the orthographic view of the cube camera, such that the cube fits in any orientation
const _CUBE_VIEW_RADIUS = 1.9
# Depth range [NDC] of the cube, in front of most of the main scene, since the depth buffer is shared
const _CUBE_DEPTH = (-1.0, -0.9)

"""
    _CubeAnimation

State of the animated transition of the main camera to a standard view, see `_set_region_view!`.
The direction from `lookat` to the eye and the up vector are rotated from `o0`, `u0` to `o1`, `u1`
within `duration` [s], while `lookat` and the distance `dist` are kept: both are rotated by `angle`
about `axis` (slerp of the direction), then the up vector is rolled by `roll` about the direction.
"""
mutable struct _CubeAnimation
    lookat::Vector{Float64}
    dist::Float64
    o0::Vector{Float64}
    u0::Vector{Float64}
    o1::Vector{Float64}
    u1::Vector{Float64}
    axis::Vector{Float64}
    angle::Float64
    roll::Float64
    duration::Float64
    elapsed::Float64
end

"""
    ViewCube

View cube returned by [`view_cube!`](@ref). The cube is drawn into the corner subscene `scene` of the
`LScene` `ls` and rotates with the camera of `ls`. Use `close` to remove it.
"""
mutable struct ViewCube
    ls::LScene
    scene::Scene
    # cube mesh, edges and labels
    faces::Vector{AbstractPlot}
    # highlight of the region under the cursor
    highlight::AbstractPlot
    highlight_mesh::Observable
    hovered::Union{Nothing, NTuple{3, Int}}
    # animated transition to a standard view
    anim::Union{Nothing, _CubeAnimation}
    duration::Float64
    margin::Float64
    # eye and up vector of the cube camera, looking at the origin
    eye::Vector{Float64}
    up::Vector{Float64}
    # a left click on the cube is in progress, its release and drags are consumed
    pressed::Bool
    listeners::Vector{Any}
    size::Int
    corner::Symbol
end

Base.show(io::IO, c::ViewCube) = print(io, "ViewCube(", c.corner, ", ", c.size, " px)")

const _CUBE_CORNERS = (:top_right, :top_left, :bottom_right, :bottom_left)

"""Returns the viewport of the view cube of edge length `size` in the `corner` of the viewport `vp`."""
function _corner_rect(vp, size::Int, corner::Symbol; margin::Int = 10)
    (x0, y0), (w, h) = minimum(vp), GeometryBasics.widths(vp)
    x = corner in (:top_right, :bottom_right) ? x0 + w - size - margin : x0 + margin
    y = corner in (:top_right, :top_left) ? y0 + h - size - margin : y0 + margin
    return Rect2i(round(Int, x), round(Int, y), size, size)
end

"""
    _ray_cube_hit(origin, dir)

Intersects the ray from `origin` along `dir` with the cube `[-1, 1]^3` (slab method). Returns the
nearest hit point with `t > 0`, i.e. the exit point if `origin` is inside of the cube, or `nothing`
if the ray misses the cube.
"""
function _ray_cube_hit(origin, dir)
    tmin, tmax = -Inf, Inf
    for i in 1:3
        if abs(dir[i]) < 1e-12
            abs(origin[i]) <= 1 || return nothing
        else
            t1, t2 = (-1 - origin[i]) / dir[i], (1 - origin[i]) / dir[i]
            tmin, tmax = max(tmin, min(t1, t2)), min(tmax, max(t1, t2))
        end
    end
    tmax < tmin && return nothing
    t = tmin > 0 ? tmin : tmax
    t > 0 || return nothing
    return Vector{Float64}(origin) .+ t .* Vector{Float64}(dir)
end

"""
    _cube_region(p; margin = 0.3)

Returns the region `s ∈ {-1, 0, 1}^3` of the point `p` on the surface of the cube `[-1, 1]^3`: a
face (one nonzero entry), an edge (two) or a corner (three). `s[i]` is the sign of `p[i]` if
`|p[i]| > 1 - margin`, otherwise `0`.
"""
function _cube_region(p; margin = 0.3)
    return ntuple(i -> abs(p[i]) > 1 - margin ? Int(sign(p[i])) : 0, 3)
end

"""
    _region_view(s)

Returns `(offset, up)` of the standard view of the cube region `s`: the unit direction from `lookat`
to the eye, i.e. the normalized sum of the adjacent face directions, and the up vector, i.e. `+z`
projected perpendicular to the `offset`, or `+y` for views along `±z`.
"""
function _region_view(s)
    o = normalize(Vector{Float64}(collect(s)))
    ez = [0.0, 0.0, 1.0]
    c = dot(ez, o)
    up = abs(c) < 0.99 ? normalize(ez .- c .* o) : [0.0, 1.0, 0.0]
    return o, up
end

"""Returns a unit vector perpendicular to the unit vector `v`."""
function _perp(v)
    a = abs(v[1]) < 0.9 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0]
    return normalize(cross(v, a))
end

"""Rotates `v` by the angle `φ` about the unit vector `k` (Rodrigues' rotation formula)."""
function _rotate(v, k, φ)
    return v .* cos(φ) .+ cross(k, v) .* sin(φ) .+ k .* (dot(k, v) * (1 - cos(φ)))
end

"""
    _rotation_between(a, b, fallback)

Returns `(axis, angle)` of the shortest rotation of the unit vector `a` onto the unit vector `b`.
If `a` and `b` are antiparallel, the rotation is about the component of `fallback` perpendicular
to `a`, such that the path does not degenerate.
"""
function _rotation_between(a, b, fallback)
    θ = acos(clamp(dot(a, b), -1, 1))
    k = cross(a, b)
    if norm(k) < 1e-9
        k = something(_orthogonalize(fallback, a), _perp(a))
    else
        k = normalize(k)
    end
    return k, θ
end

"""Returns the unit component of `u` perpendicular to the unit vector `d`, or `nothing`."""
function _orthogonalize(u, d)
    v = u .- dot(u, d) .* d
    n = norm(v)
    return n < 1e-6 ? nothing : v ./ n
end

"""
    _region_mesh(s; margin = 0.3, offset = 0.01)

Returns the mesh of the region `s` on the surface of the cube `[-1, 1]^3`, i.e. one rectangle on
each face with a nonzero `s[i]`, lifted by `offset` outward.
"""
function _region_mesh(s; margin = 0.3, offset = 0.01)
    pts = Point3f[]
    faces = GLTriangleFace[]
    interval(si) = si == 0 ? (-1 + margin, 1 - margin) : si > 0 ? (1 - margin, 1.0) : (-1.0, -1 + margin)
    for i in 1:3
        s[i] == 0 && continue
        j, k = mod1(i + 1, 3), mod1(i + 2, 3)
        (a0, a1), (b0, b1) = interval(s[j]), interval(s[k])
        n = length(pts)
        for (a, b) in ((a0, b0), (a1, b0), (a1, b1), (a0, b1))
            p = zeros(3)
            p[i], p[j], p[k] = s[i] * (1 + offset), a, b
            push!(pts, Point3f(p))
        end
        push!(faces, GLTriangleFace(n + 1, n + 2, n + 3), GLTriangleFace(n + 1, n + 3, n + 4))
    end
    return GeometryBasics.Mesh(pts, faces)
end

"""Returns the ray of the cube camera through the pixel `xy`, relative to the cube viewport."""
function _cube_ray(cube::ViewCube, xy)
    w, h = GeometryBasics.widths(cube.scene.viewport[])
    u_z = normalize(-cube.eye)
    u_x = normalize(cross(u_z, cube.up))
    u_y = cross(u_x, u_z)
    rx, ry = 2 * xy[1] / w - 1, 2 * xy[2] / h - 1
    r = _CUBE_VIEW_RADIUS
    origin = cube.eye .+ r .* (rx .* u_x .+ ry .* u_y)
    return origin, u_z
end

"""Returns the region of the cube under the cursor, or `nothing` if the cursor misses the cube."""
function _region_at_cursor(cube::ViewCube)
    mp = events(cube.ls.scene).mouseposition[]
    vp = cube.scene.viewport[]
    Makie.Vec(mp) in vp || return nothing
    origin, dir = _cube_ray(cube, Makie.Vec2(mp) .- minimum(vp))
    p = _ray_cube_hit(origin, dir)
    return isnothing(p) ? nothing : _cube_region(p; margin = cube.margin)
end

"""Sets the matrices of the cube camera to the eye `4 · o` and the up vector `up`, looking at 0."""
function _set_cube_camera!(cube::ViewCube, o, up)
    cube.eye = _CUBE_EYE_DIST .* o
    cube.up = up
    r = _CUBE_VIEW_RADIUS
    view = Makie.lookat(Makie.Vec3d(cube.eye...), Makie.Vec3d(0, 0, 0), Makie.Vec3d(up...))
    proj = Makie.orthographicprojection(-r, r, -r, r, _CUBE_EYE_DIST - 2.5, _CUBE_EYE_DIST + 2.5)
    # Compresses the depth into a thin slice at the near plane
    (z0, z1) = _CUBE_DEPTH
    a, b = (z1 - z0) / 2, (z1 + z0) / 2
    S = Makie.Mat4d(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, a, 0, 0, 0, b, 1)
    cam = cube.scene.camera
    Makie.set_proj_view!(cam, S * proj, view)
    cam.eyeposition[] = Vec3f(cube.eye...)
    cam.upvector[] = Vec3f(up...)
    cam.view_direction[] = Vec3f(normalize(-cube.eye)...)
    return nothing
end

"""Rotates the cube camera to the orientation of the main camera of the `LScene`."""
function _sync_cube_camera!(cube::ViewCube)
    cam = cameracontrols(cube.ls.scene)
    d = Vector{Float64}(cam.eyeposition[]) .- Vector{Float64}(cam.lookat[])
    norm(d) > 0 || return nothing
    o = normalize(d)
    # Transient states while `set_view` updates the observables one by one are skipped
    up = _orthogonalize(Vector{Float64}(cam.upvector[]), o)
    isnothing(up) && return nothing
    _set_cube_camera!(cube, o, up)
    return nothing
end

"""Highlights the region `s` of the cube, or nothing."""
function _set_hover!(cube::ViewCube, s)
    s == cube.hovered && return nothing
    cube.hovered = s
    if isnothing(s)
        cube.highlight.visible[] = false
    else
        cube.highlight_mesh[] = _region_mesh(s; margin = cube.margin)
        cube.highlight.visible[] = true
    end
    return nothing
end

"""Applies the frame at the fraction `t` of the animation `a` to the main camera."""
function _apply_frame!(cube::ViewCube, a::_CubeAnimation, t)
    if t >= 1
        o, up = a.o1, a.u1
    else
        # Smooth start and stop
        τ = t^2 * (3 - 2t)
        o = normalize(_rotate(a.o0, a.axis, τ * a.angle))
        # The up vector is rotated along, hence it stays perpendicular to the direction
        u = _rotate(a.u0, a.axis, τ * a.angle)
        up = normalize(_rotate(u, o, τ * a.roll))
    end
    set_view(cube.ls, a.lookat .+ a.dist .* o, a.lookat, up)
    return nothing
end

"""Advances the animation of the `cube` by `dt` [s]."""
function _step_animation!(cube::ViewCube, dt)
    a = cube.anim
    isnothing(a) && return nothing
    a.elapsed += dt
    t = a.duration > 0 ? min(a.elapsed / a.duration, 1.0) : 1.0
    t >= 1 && (cube.anim = nothing)
    _apply_frame!(cube, a, t)
    return nothing
end

"""
    _set_region_view!(cube, s)

Moves the main camera to the standard view of the region `s`, keeping `lookat` and the distance of
the eye. Animated over `cube.duration`, driven by `tick` events, instant if the duration is `0`.
"""
function _set_region_view!(cube::ViewCube, s)
    cam = cameracontrols(cube.ls.scene)
    lookat = Vector{Float64}(cam.lookat[])
    d = Vector{Float64}(cam.eyeposition[]) .- lookat
    dist = norm(d)
    o0 = dist > 0 ? d ./ dist : [0.0, 0.0, 1.0]
    u0 = something(_orthogonalize(Vector{Float64}(cam.upvector[]), o0), _perp(o0))
    o1, u1 = _region_view(s)
    # Opposite views rotate via the current up vector, i.e. over the top of the screen
    axis, angle = _rotation_between(o0, o1, cross(o0, u0))
    # Remaining roll about the new direction from the rotated to the new up vector
    ur = _rotate(u0, axis, angle)
    roll = atan(dot(cross(ur, u1), o1), dot(ur, u1))
    # A new click during an animation starts from the current state
    cube.anim = _CubeAnimation(lookat, dist, o0, u0, o1, u1, axis, angle, roll, cube.duration, 0.0)
    cube.duration > 0 || _step_animation!(cube, 0.0)
    return nothing
end

"""
    view_cube!(ls::LScene; size = 110, corner = :top_right, duration = 0.3)

Adds a view cube in the `corner` of the `LScene` `ls`, see the docstring in `Render.jl`. Returns a
`ViewCube`, which can be removed via `close`.
"""
function view_cube!(
        ls::LScene;
        size::Integer = 110,
        corner::Symbol = :top_right,
        duration::Real = 0.3
    )
    corner in _CUBE_CORNERS || throw(ArgumentError("corner must be one of $_CUBE_CORNERS, got :$corner"))
    size > 0 || throw(ArgumentError("size must be positive, got $size"))
    duration >= 0 || throw(ArgumentError("duration must not be negative, got $duration"))
    main = ls.scene
    cameracontrols(main) isa Makie.Camera3D ||
        throw(ArgumentError("view_cube! requires an LScene with a Camera3D"))
    size = Int(size)
    listeners = Any[]

    # Subscene without camera controls, the camera matrices are set directly
    vp = Observable(_corner_rect(main.viewport[], size, corner); ignore_equal_values = true)
    push!(listeners, on(v -> (vp[] = _corner_rect(v, size, corner)), main.viewport))
    scene = Scene(main; viewport = vp, clear = false)

    # Cube with labels on the faces and a colored axis triad at the corner (-1, -1, -1)
    edge_pts = Point3f[]
    edge_colors = RGBAf[]
    dark = RGBAf(0.25, 0.25, 0.25, 1)
    for i in 1:3, a in (-1, 1), b in (-1, 1)
        j, k = mod1(i + 1, 3), mod1(i + 2, 3)
        p, q = zeros(3), zeros(3)
        p[i], q[i] = -1, 1
        p[j], p[k] = a, b
        q[j], q[k] = a, b
        append!(edge_pts, (Point3f(p), Point3f(q)))
        color = (p[j] == -1 && p[k] == -1) ? RGBAf(Makie.to_color((:red, :green, :blue)[i])) : dark
        append!(edge_colors, (color, color))
    end
    names = [f[1] for f in _CUBE_FACES]
    label_pos = [Point3f(1.02 .* f[2]) for f in _CUBE_FACES]
    label_rot = map(_CUBE_FACES) do (_, n)
        # Readable from outside with the up vector of the standard view of the face
        o, up = _region_view(n)
        right = cross(up, o)
        return _quat_from_rotmatrix(Float32.(hcat(right, up, o)))
    end
    # The cube is never clipped by the clip planes of the main scene
    faces = AbstractPlot[
        mesh!(scene, Rect3f(Point3f(-1), Vec3f(2)); color = RGBf(0.88, 0.88, 0.88),
            clip_planes = Plane3f[]),
        linesegments!(scene, edge_pts; color = edge_colors, linewidth = 2, depth_shift = -1f-4,
            clip_planes = Plane3f[]),
        text!(scene, label_pos; text = names, rotation = collect(label_rot), markerspace = :data,
            fontsize = 0.42, align = (:center, :center), color = :black, clip_planes = Plane3f[])
    ]
    highlight_mesh = Observable(_region_mesh((0, 0, 1)))
    highlight = mesh!(scene, highlight_mesh; color = RGBAf(0.2, 0.5, 1.0, 0.6), transparency = true,
        visible = false, clip_planes = Plane3f[])

    cube = ViewCube(ls, scene, faces, highlight, highlight_mesh, nothing, nothing, Float64(duration),
        0.3, zeros(3), zeros(3), false, listeners, size, corner)
    cam = cameracontrols(main)
    for obs in (cam.eyeposition, cam.lookat, cam.upvector)
        push!(listeners, on(_ -> _sync_cube_camera!(cube), obs))
    end
    _sync_cube_camera!(cube)
    iszero(cube.eye) && _set_cube_camera!(cube, normalize([1.0, -1.0, 1.0]), _region_view((1, -1, 1))[2])

    # Before the camera and the kinematic controls (priority 200), see `kinematic_controls!`
    push!(listeners, on(events(main).mousebutton, priority = 300) do event
        event.button == Mouse.left || return Consume(false)
        if event.action == Mouse.press
            s = _region_at_cursor(cube)
            isnothing(s) && return Consume(false)
            cube.pressed = true
            _set_region_view!(cube, s)
            return Consume(true)
        elseif event.action == Mouse.release && cube.pressed
            cube.pressed = false
            return Consume(true)
        end
        return Consume(false)
    end)
    push!(listeners, on(events(main).mouseposition, priority = 300) do _
        _set_hover!(cube, _region_at_cursor(cube))
        # Drags that started on the cube do not reach the camera
        return Consume(cube.pressed)
    end)
    push!(listeners, on(tick -> _step_animation!(cube, tick.delta_time), events(main).tick))
    return cube
end

function Base.close(cube::ViewCube)
    foreach(off, cube.listeners)
    empty!(cube.listeners)
    cube.anim = nothing
    # Closed already
    isnothing(cube.scene.parent) && return nothing
    empty!(cube.faces)
    # Deletes the plots and removes the subscene from the scene and the screens
    Makie.free(cube.scene)
    return nothing
end
