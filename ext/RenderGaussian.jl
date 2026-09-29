#=
Geometry and plot function of the Gaussian beamlets, shared with `RenderAstigmaticGaussian.jl`
=#

# Generating beams of the overlays of `show_beams`, an `AstigmaticGaussianBeamlet` has 4 divergence
# and waist beams
const _GAUSS_BEAM_GROUPS = (:chief, :divergence, :waist)
const _GAUSS_BEAM_COLORS = (:red, :green, :blue)
const _AGB_BEAM_FIELD_GROUPS = (chief = (:c,), divergence = (:dxp, :dxm, :dyp, :dym), waist = (:wxp, :wxm, :wyp, :wym))

"""
    _generating_beams(beamlet, group)

The generating beams of one segment tree node of a Gaussian `beamlet` that make up the `group`
(`:chief`, `:divergence` or `:waist`) of the `show_beams` overlay.
"""
_generating_beams(gauss::BMO.GaussianBeamlet, group::Symbol) = (getfield(gauss, group),)
_generating_beams(agb::BMO.AstigmaticGaussianBeamlet, group::Symbol) =
    map(f -> getfield(agb, f), getfield(_AGB_BEAM_FIELD_GROUPS, group))

"""
    _beamlets(thing, render_every = 1)

The Gaussian beamlets drawn of `thing`: the beamlet itself or every `render_every`-th beamlet of an
[`BMO.AstigmaticBeamGroup`](@ref).
"""
_beamlets(beamlet::Union{BMO.GaussianBeamlet, BMO.AstigmaticGaussianBeamlet}, render_every = 1) = (beamlet,)
_beamlets(bg::BMO.AstigmaticBeamGroup, render_every = 1) = view(BMO.beams(bg), 1:render_every:length(BMO.beams(bg)))

"""
    _foreach_generating_ray(f, thing, group; render_every = 1)

Calls `f(ray)` for all rays of the generating beams of the `group` (see
[`_generating_beams`](@ref)) of all segments of the beamlets of `thing` (see [`_beamlets`](@ref)).
"""
function _foreach_generating_ray(f, thing, group::Symbol; render_every::Int = 1)
    for beamlet in _beamlets(thing, render_every)
        for child in PreOrderDFS(beamlet)
            for beam in _generating_beams(child, group)
                foreach(f, BMO.rays(beam))
            end
        end
    end
    return nothing
end

"""
    _generating_segments(thing, group; flen, render_every = 1) -> Vector{Point3f}

Segments of the generating rays of the `group` of the Gaussian `thing`, see
[`_foreach_generating_ray`](@ref).
"""
function _generating_segments(thing, group::Symbol; flen, render_every::Int = 1)
    pts = Point3f[]
    _foreach_generating_ray(ray -> _push_ray_segment!(pts, ray; flen), thing, group; render_every)
    return pts
end

"""
    _generating_ends(thing, group; render_every = 1) -> Vector{Point3f}

Start points (and end points with intersection) of the generating rays of the `group` of the
Gaussian `thing`, see `show_pos`.
"""
function _generating_ends(thing, group::Symbol; render_every::Int = 1)
    pts = Point3f[]
    _foreach_generating_ray(ray -> _push_ray_ends!(pts, ray), thing, group; render_every)
    return pts
end

"""
    _gaussian_mesh(gauss::GaussianBeamlet; flen, r_res, z_res)

Returns a single mesh of the 1/e² envelope of all chief ray segments of the `gauss`ian beamlet,
see the `render!` method of the [`GaussianBeamlet`](@ref).
"""
function _gaussian_mesh(gauss::BMO.GaussianBeamlet{T}; flen, r_res::Int, z_res::Int) where {T}
    pts = Point3f[]
    faces = GLTriangleFace[]
    vs = _ring_angles(r_res)
    for child in PreOrderDFS(gauss)
        # Length tracking variable
        l = isnothing(child.parent) ? zero(T) : length(child.parent)
        for ray in BMO.rays(child.chief)
            u = LinRange(0, isnothing(BMO.intersection(ray)) ? flen : length(ray), z_res)
            w = BMO.gauss_parameters(child, u .+ l)[1]
            R = BMO.align3d([0, 1, 0], ray.dir)
            p = position(ray)
            offset = length(pts)
            # Beam surface along the local y-axis, transformed into world coords. The rings run
            # counterclockwise about +y, i.e. the faces point outwards, see `_push_grid_faces!`
            for (i, ui) in enumerate(u), v in vs
                push!(pts, Point3f(R * Point3(w[i] * cos(v), ui, -w[i] * sin(v)) + p))
            end
            _push_grid_faces!(faces, offset, r_res, z_res)
            if !isnothing(BMO.intersection(ray))
                l += length(ray)
            end
        end
    end
    return Mesh(pts, faces)
end

"""
    _push_grid_faces!(faces, offset, r_res, z_res)

Pushes the triangles of a tube of `z_res` rings of `r_res` vertices each (radial index fastest),
which start at index `offset + 1`, onto `faces`. Each ring is closed, i.e. its last vertex is
connected to its first, see `_ring_angles`. The faces point outwards if the rings run
counterclockwise about the direction in which the rings follow each other.
"""
function _push_grid_faces!(faces::Vector{GLTriangleFace}, offset, r_res, z_res)
    for i in 1:(z_res - 1), j in 1:r_res
        a = offset + (i - 1) * r_res + j
        a1 = offset + (i - 1) * r_res + mod1(j + 1, r_res)
        c = a + r_res
        c1 = a1 + r_res
        push!(faces, GLTriangleFace(a, a1, c))
        push!(faces, GLTriangleFace(a1, c1, c))
    end
    return nothing
end

"""
    _ring_angles(r_res)

The `r_res` angles [rad] of the vertices of a ring of the envelope mesh, without the duplicate of
the first vertex at 2π: `_push_grid_faces!` closes the ring, such that the vertices on the seam
are shared and its normals are smooth.
"""
_ring_angles(r_res::Int) = 2π .* (0:(r_res - 1)) ./ r_res

"""
    _plot_generating_beams!(bp, axis, thing; flen, render_every, show_pos, transparency)

Adds the `show_beams` overlay of the Gaussian `thing` to `bp`: one `linesegments` plot per group of
generating beams (chief: red, divergence: green, waist: blue) and, with `show_pos`, one `scatter`
plot of their positions per group.
"""
function _plot_generating_beams!(bp::_BeamPlots, axis::_RenderEnv, thing; flen, render_every::Int = 1,
        show_pos, transparency)
    for (group, color) in zip(_GAUSS_BEAM_GROUPS, _GAUSS_BEAM_COLORS)
        segments = _observe!(bp, () -> _generating_segments(thing, group; flen, render_every))
        _add!(bp, linesegments!(axis, segments; color, transparency))
        if show_pos
            ends = _observe!(bp, () -> _generating_ends(thing, group; render_every))
            _add!(bp, scatter!(axis, ends; color))
        end
    end
    return bp
end

"""
    _plot_gaussian!(axis, gauss::GaussianBeamlet; kwargs...) -> _BeamPlots

The plot function of the [`GaussianBeamlet`](@ref): the envelope of all segments in a single
`mesh` plot and the overlay of `show_beams`. See `render!(axis, gauss)` for the kwargs.
"""
function _plot_gaussian!(
        axis::_RenderEnv,
        gauss::BMO.GaussianBeamlet;
        show_beams,
        show_pos,
        r_res::Int,
        z_res::Int,
        flen,
        color,
        transparency,
        kwargs...
    )
    bp = _BeamPlots()
    # the mesh type depends on its size, hence Any
    geometry = () -> _gaussian_mesh(gauss; flen, r_res, z_res)
    envelope = _observe!(bp, geometry, Observable{Any}(geometry()))
    _add!(bp, mesh!(axis, envelope; color, transparency, kwargs...))
    show_beams && _plot_generating_beams!(bp, axis, gauss; flen, show_pos, transparency)
    return bp
end

"""
    render!(ax, gauss::GaussianBeamlet; kwargs...)

Render the 1/e² envelope of the `GaussianBeamlet` into the specified `axis`, all segments as a
single mesh.

With `show_beams = true` the generating rays are overlayed into the axis as follows:

- `chief` beam: red
- `divergence` beam: green
- `waist` beam: blue

# Keyword args

- `show_beams = false`: plot the generating rays of the [`GaussianBeamlet`](@ref)
- `show_pos = false`: with `show_beams`, marks the positions of the generating rays
- `flen = 0.1`: length of the final beam in case of no intersection
- `r_res::Int = 50`: radial resolution of the beam
- `z_res::Int = 100`: resolution along the optical axis of the beam

# Makie kwargs

- `color = :red`
- `transparency = true`

Additional kwargs are passed into the mesh plot of the Gaussian envelope.
"""
function render!(axis::_RenderEnv, gauss::BMO.GaussianBeamlet; r_res::Int = 50, z_res::Int = 100, kwargs...)
    _plot_gaussian!(axis, gauss; _gaussian_defaults(kwargs)..., r_res, z_res)
    return nothing
end

"""The kwargs of `render!` and `live_render!` of a `GaussianBeamlet` except `r_res`, `z_res`."""
_gaussian_defaults(kwargs) = (; show_beams = false, show_pos = false, flen = 0.1, color = :red,
    transparency = true, kwargs...)

"""
    live_render!(axis, gauss::GaussianBeamlet; kwargs...)

Live-renders the 1/e² envelope of the [`GaussianBeamlet`](@ref) like `render!`, as a single mesh,
see [`live_render!`](@ref). The kwargs are those of `render!`, with a coarser default resolution:

- `r_res::Int = 24`: radial resolution of the beam
- `z_res::Int = 40`: resolution along the optical axis of the beam
"""
function live_render!(axis::_RenderEnv, gauss::BMO.GaussianBeamlet; r_res::Int = 24, z_res::Int = 40, kwargs...)
    kw = _gaussian_defaults(kwargs)
    bp = _plot_gaussian!(axis, gauss; kw..., r_res, z_res)
    return BeamRenderHandle(gauss, axis, bp, (; flen = Float64(kw.flen), render_every = 1, r_res, z_res))
end
