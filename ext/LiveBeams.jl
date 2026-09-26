using Makie: linesegments!, Observable, notify
import GeometryBasics
using GeometryBasics: GLTriangleFace

"""
    BeamRenderHandle{B} <: AbstractRenderHandle

Live rendering handle of an `AbstractRay`, a [`Beam`](@ref) or an [`BeamletOptics.AbstractBeamGroup`](@ref).
All ray segments are drawn by a single `linesegments` plot.
"""
mutable struct BeamRenderHandle{B} <: AbstractRenderHandle
    thing::B
    axis::_RenderEnv
    plot::AbstractPlot
    points::Observable{Vector{Point3f}}
    flen::Float64
    render_every::Int
end

function Base.show(io::IO, h::BeamRenderHandle{B}) where {B}
    print(io, "BeamRenderHandle{", B, "}(", length(h.points[]) ÷ 2, " segments)")
end

function _push_ray_segment!(pts::Vector{Point3f}, ray::BMO.AbstractRay; flen)
    isect = BMO.intersection(ray)
    len = isnothing(isect) ? flen : length(isect)
    push!(pts, Point3f(position(ray)))
    push!(pts, Point3f(position(ray) + len * BMO.direction(ray)))
    return nothing
end

_collect_segments!(pts, ray::BMO.AbstractRay; flen, kwargs...) = _push_ray_segment!(pts, ray; flen)

function _collect_segments!(pts, beam::Beam; flen, kwargs...)
    for child in PreOrderDFS(beam)
        for ray in BMO.rays(child)
            _push_ray_segment!(pts, ray; flen)
        end
    end
    return nothing
end

function _collect_segments!(pts, bg::BMO.AbstractBeamGroup; flen, render_every::Int = 1)
    bms = BMO.beams(bg)
    if !isempty(bms) && !(first(bms) isa Beam)
        throw(ArgumentError("live_render! via linesegments is only implemented for beam groups of " *
            "`Beam`s; groups of `AstigmaticGaussianBeamlet`s are rendered via the dedicated " *
            "`live_render!(axis, ::AstigmaticBeamGroup; ...)` mesh method"))
    end
    for i in 1:render_every:length(bms)
        _collect_segments!(pts, bms[i]; flen)
    end
    return nothing
end

"""
    live_render!(axis, beam; kwargs...)

Live-renders a ray, [`Beam`](@ref) or beam group as a single `linesegments` plot, see
[`live_render!`](@ref). The number of segments may change between updates.

# Keyword args

- `flen = 1.0`: length of the final segment in case of no intersection [m]
- `render_every = 5`: renders only every e.g. fifth beam of a beam group

# Makie kwargs

- `color = :blue`
- `linewidth = 1.0`
- `transparency = true`
"""
function live_render!(
        axis::_RenderEnv,
        thing::Union{BMO.AbstractRay, Beam, BMO.AbstractBeamGroup};
        # kwargs
        flen = 1.0,
        render_every::Int = 5,
        # Makie kwargs
        color = :blue,
        linewidth = 1.0,
        transparency = true,
        kwargs...
    )
    pts = Point3f[]
    _collect_segments!(pts, thing; flen, render_every)
    obs = Observable(pts)
    plot = linesegments!(axis, obs; color, linewidth, transparency, kwargs...)
    return BeamRenderHandle(thing, axis, plot, obs, Float64(flen), render_every)
end

function update_render!(h::BeamRenderHandle)
    empty!(h.points[])
    _collect_segments!(h.points[], h.thing; flen = h.flen, render_every = h.render_every)
    notify(h.points)
    return h
end

remove_render!(h::BeamRenderHandle) = (delete!(h.axis, h.plot); return nothing)

"""
    GaussianRenderHandle{G} <: AbstractRenderHandle

Live rendering handle of a [`GaussianBeamlet`](@ref). The 1/e² envelope of all segments is drawn by
a single `mesh` plot. The generating rays can optionally be shown as three additional plots.
"""
mutable struct GaussianRenderHandle{G} <: AbstractRenderHandle
    thing::G
    axis::_RenderEnv
    mesh_plot::AbstractPlot
    # the mesh type depends on its size, hence Any
    mesh_obs::Observable{Any}
    flen::Float64
    r_res::Int
    z_res::Int
    # chief, divergence and waist beam overlays, empty if show_beams = false
    beam_plots::Vector{AbstractPlot}
    beam_obs::Vector{Observable{Vector{Point3f}}}
end

function Base.show(io::IO, h::GaussianRenderHandle{G}) where {G}
    m = h.mesh_obs[]
    print(io, "GaussianRenderHandle{", G, "}(", length(GeometryBasics.coordinates(m)), " vertices, ",
        length(GeometryBasics.faces(m)), " faces)")
end

const _GAUSS_BEAM_FIELDS = (:chief, :divergence, :waist)

function _collect_beamlet_segments!(pts, gauss::BMO.GaussianBeamlet, field::Symbol; flen)
    for child in PreOrderDFS(gauss)
        for ray in BMO.rays(getfield(child, field))
            _push_ray_segment!(pts, ray; flen)
        end
    end
    return nothing
end

"""
    _gaussian_mesh(gauss; flen, r_res, z_res)

Returns a single mesh of the 1/e² envelope of all chief ray segments of the `gauss`ian beamlet,
see the `render!` method of the [`GaussianBeamlet`](@ref).
"""
function _gaussian_mesh(gauss::BMO.GaussianBeamlet{T}; flen, r_res::Int, z_res::Int) where {T}
    pts = Point3f[]
    faces = GLTriangleFace[]
    vs = LinRange(0, 2π, r_res)
    for child in PreOrderDFS(gauss)
        # Length tracking variable
        l = isnothing(child.parent) ? zero(T) : length(child.parent)
        for ray in BMO.rays(child.chief)
            u = LinRange(0, isnothing(BMO.intersection(ray)) ? flen : length(ray), z_res)
            w = BMO.gauss_parameters(child, u .+ l)[1]
            R = BMO.align3d([0, 1, 0], ray.dir)
            p = position(ray)
            offset = length(pts)
            # Beam surface along the local y-axis, transformed into world coords
            for (i, ui) in enumerate(u), v in vs
                push!(pts, Point3f(R * Point3(w[i] * cos(v), ui, w[i] * sin(v)) + p))
            end
            for i in 1:(z_res - 1), j in 1:(r_res - 1)
                a = offset + (i - 1) * r_res + j
                c = offset + i * r_res + j
                push!(faces, GLTriangleFace(a, a + 1, c))
                push!(faces, GLTriangleFace(a + 1, c + 1, c))
            end
            if !isnothing(BMO.intersection(ray))
                l += length(ray)
            end
        end
    end
    return GeometryBasics.Mesh(pts, faces)
end

"""
    live_render!(axis, gauss::GaussianBeamlet; kwargs...)

Live-renders the 1/e² envelope of the [`GaussianBeamlet`](@ref) as a single mesh, see
[`live_render!`](@ref).

With `show_beams = true` the generating rays are overlayed into the axis as follows:

- `chief` beam: red
- `divergence` beam: green
- `waist` beam: blue

# Keyword args

- `show_beams = false`: plot the generating rays of the [`GaussianBeamlet`](@ref)
- `flen = 0.1`: length of the final beam in case of no intersection
- `r_res::Int = 24`: radial resolution of the beam
- `z_res::Int = 40`: resolution along the optical axis of the beam

# Makie kwargs

- `color = :red`
- `transparency = true`
"""
function live_render!(
        axis::_RenderEnv,
        gauss::BMO.GaussianBeamlet;
        # kwargs
        show_beams = false,
        r_res::Int = 24,
        z_res::Int = 40,
        flen = 0.1,
        # Makie kwargs
        color = :red,
        transparency = true,
        kwargs...
    )
    mesh_obs = Observable{Any}(_gaussian_mesh(gauss; flen, r_res, z_res))
    mesh_plot = mesh!(axis, mesh_obs; color, transparency, kwargs...)
    beam_plots = AbstractPlot[]
    beam_obs = Observable{Vector{Point3f}}[]
    if show_beams
        for (field, c) in zip(_GAUSS_BEAM_FIELDS, (:red, :green, :blue))
            pts = Point3f[]
            _collect_beamlet_segments!(pts, gauss, field; flen)
            obs = Observable(pts)
            push!(beam_obs, obs)
            push!(beam_plots, linesegments!(axis, obs; color = c))
        end
    end
    return GaussianRenderHandle(gauss, axis, mesh_plot, mesh_obs, Float64(flen), r_res, z_res,
        beam_plots, beam_obs)
end

function update_render!(h::GaussianRenderHandle)
    h.mesh_obs[] = _gaussian_mesh(h.thing; flen = h.flen, r_res = h.r_res, z_res = h.z_res)
    for (field, obs) in zip(_GAUSS_BEAM_FIELDS, h.beam_obs)
        empty!(obs[])
        _collect_beamlet_segments!(obs[], h.thing, field; flen = h.flen)
        notify(obs)
    end
    return h
end

function remove_render!(h::GaussianRenderHandle)
    delete!(h.axis, h.mesh_plot)
    foreach(p -> delete!(h.axis, p), h.beam_plots)
    empty!(h.beam_plots)
    return nothing
end

# Beams of the overlays, an `AstigmaticGaussianBeamlet` has 4 divergence and waist beams
const _AGB_BEAM_FIELD_GROUPS = (chief = (:c,), divergence = (:dxp, :dxm, :dyp, :dym), waist = (:wxp, :wxm, :wyp, :wym))

function _collect_beamlet_segments!(pts, agb::BMO.AstigmaticGaussianBeamlet, group::Symbol; flen)
    for child in PreOrderDFS(agb)
        for fname in getfield(_AGB_BEAM_FIELD_GROUPS, group)
            for ray in BMO.rays(getfield(child, fname))
                _push_ray_segment!(pts, ray; flen)
            end
        end
    end
    return nothing
end

"""
    _smooth_astigmatic_axes!(bs, cs)

Applies parallel transport to the elliptical basis vectors `bs`, `cs` (one pair per
longitudinal sample) so that they vary smoothly along a segment, preserving handedness
and preventing mesh twisting/flips near a focus, see the static `render!` method of
[`BMO.AstigmaticGaussianBeamlet`](@ref).
"""
function _smooth_astigmatic_axes!(bs, cs)
    for j in 2:length(bs)
        b_new = bs[j]
        c_new = cs[j]
        b_old = bs[j - 1]
        c_old = cs[j - 1]

        # Preserve the handedness, avoids inside-out meshes
        if dot(cross(b_new, c_new), cross(b_old, c_old)) < 0
            b_new = -b_new
        end

        # Align with the previous slice, avoids twisting when the axes swap
        bn, cn, bo, co = normalize(b_new), normalize(c_new), normalize(b_old), normalize(c_old)
        phi = atan(dot(cn, bo) - dot(bn, co), dot(bn, bo) + dot(cn, co))
        cos_phi, sin_phi = cos(phi), sin(phi)

        bs[j] = b_new * cos_phi + c_new * sin_phi
        cs[j] = -b_new * sin_phi + c_new * cos_phi
    end
    return nothing
end

"""
    _append_astigmatic_mesh!(pts, faces, agb; flen, r_res, z_res)

Appends the 1/e² envelope mesh of `agb` (all of its segments) to the shared `pts`/`faces`
buffers, see the `render!` method of [`BMO.AstigmaticGaussianBeamlet`](@ref).
"""
function _append_astigmatic_mesh!(
        pts::Vector{Point3f}, faces::Vector{GLTriangleFace}, agb::BMO.AstigmaticGaussianBeamlet{T};
        flen, r_res::Int, z_res::Int
    ) where {T}
    vs = LinRange(0, 2π, r_res)
    for child in PreOrderDFS(agb)
        p = child.parent
        l = isnothing(p) ? zero(T) : length(p)
        for ray in BMO.rays(child.c)
            l_local = isnothing(BMO.intersection(ray)) ? flen : length(ray)
            us = LinRange(0, l_local, z_res) .+ l
            p0s, bs, cs, = BMO.waist_parameters(child, us)
            _smooth_astigmatic_axes!(bs, cs)
            offset = length(pts)
            for i in eachindex(p0s), v in vs
                push!(pts, Point3f(BMO.ellipse(v, p0s[i], bs[i], cs[i])))
            end
            for i in 1:(z_res - 1), j in 1:(r_res - 1)
                a = offset + (i - 1) * r_res + j
                c2 = offset + i * r_res + j
                push!(faces, GLTriangleFace(a, a + 1, c2))
                push!(faces, GLTriangleFace(a + 1, c2 + 1, c2))
            end
            if !isnothing(BMO.intersection(ray))
                l += length(ray)
            end
        end
    end
    return nothing
end

"""
    _gaussian_mesh(agb::AstigmaticGaussianBeamlet; flen, r_res, z_res)

Returns a single mesh of the 1/e² envelope of all chief ray segments of the `agb`, see
the `render!` method of [`BMO.AstigmaticGaussianBeamlet`](@ref).
"""
function _gaussian_mesh(agb::BMO.AstigmaticGaussianBeamlet; flen, r_res::Int, z_res::Int)
    pts = Point3f[]
    faces = GLTriangleFace[]
    _append_astigmatic_mesh!(pts, faces, agb; flen, r_res, z_res)
    return GeometryBasics.Mesh(pts, faces)
end

"""
    live_render!(axis, agb::AstigmaticGaussianBeamlet; kwargs...)

Live-renders the 1/e² envelope of the [`BMO.AstigmaticGaussianBeamlet`](@ref) as a single
mesh, see [`live_render!`](@ref). Polarization overlay kwargs of the static `render!`
method (`show_polarization`, `pol_*`) are not supported and are silently ignored via
`kwargs...`.

With `show_beams = true` the generating rays are overlayed into the axis as follows:

- chief ray: red
- divergence rays (`dxp`, `dxm`, `dyp`, `dym`): green
- waist rays (`wxp`, `wxm`, `wyp`, `wym`): blue

# Keyword args

- `show_beams = false`: plot the generating rays of the [`BMO.AstigmaticGaussianBeamlet`](@ref)
- `flen = 0.1`: length of the final beam in case of no intersection
- `r_res::Int = 64`: radial (angular) resolution of the beam
- `z_res::Int = 100`: resolution along the optical axis of the beam

# Makie kwargs

- `color = :red`
- `transparency = true`
"""
function live_render!(
        axis::_RenderEnv,
        agb::BMO.AstigmaticGaussianBeamlet;
        # kwargs
        show_beams = false,
        r_res::Int = 64,
        z_res::Int = 100,
        flen = 0.1,
        # Makie kwargs
        color = :red,
        transparency = true,
        kwargs...
    )
    mesh_obs = Observable{Any}(_gaussian_mesh(agb; flen, r_res, z_res))
    mesh_plot = mesh!(axis, mesh_obs; color, transparency, kwargs...)
    beam_plots = AbstractPlot[]
    beam_obs = Observable{Vector{Point3f}}[]
    if show_beams
        for (field, c) in zip(_GAUSS_BEAM_FIELDS, (:red, :green, :blue))
            pts = Point3f[]
            _collect_beamlet_segments!(pts, agb, field; flen)
            obs = Observable(pts)
            push!(beam_obs, obs)
            push!(beam_plots, linesegments!(axis, obs; color = c))
        end
    end
    return GaussianRenderHandle(agb, axis, mesh_plot, mesh_obs, Float64(flen), r_res, z_res,
        beam_plots, beam_obs)
end

"""
    AstigmaticGroupRenderHandle{G} <: AbstractRenderHandle

Live rendering handle of an [`BMO.AstigmaticBeamGroup`](@ref). The 1/e² envelopes of every
`render_every`-th beamlet are merged into a single `mesh` plot.
"""
mutable struct AstigmaticGroupRenderHandle{G} <: AbstractRenderHandle
    thing::G
    axis::_RenderEnv
    mesh_plot::AbstractPlot
    mesh_obs::Observable{Any}
    flen::Float64
    r_res::Int
    z_res::Int
    render_every::Int
end

function Base.show(io::IO, h::AstigmaticGroupRenderHandle{G}) where {G}
    m = h.mesh_obs[]
    print(io, "AstigmaticGroupRenderHandle{", G, "}(", length(GeometryBasics.coordinates(m)), " vertices, ",
        length(GeometryBasics.faces(m)), " faces)")
end

function _astigmatic_group_mesh(bg::BMO.AstigmaticBeamGroup; flen, r_res::Int, z_res::Int, render_every::Int)
    pts = Point3f[]
    faces = GLTriangleFace[]
    bms = BMO.beams(bg)
    for i in 1:render_every:length(bms)
        _append_astigmatic_mesh!(pts, faces, bms[i]; flen, r_res, z_res)
    end
    return GeometryBasics.Mesh(pts, faces)
end

"""
    live_render!(axis, bg::AstigmaticBeamGroup; kwargs...)

Live-renders every `render_every`-th beamlet envelope of the [`BMO.AstigmaticBeamGroup`](@ref)
as a single merged mesh, see [`live_render!`](@ref). The default resolution is deliberately
low (compared to a single [`BMO.AstigmaticGaussianBeamlet`](@ref)) since a group typically
contains many beamlets.

# Keyword args

- `render_every::Int = 5`: renders only every e.g. fifth beamlet of the group
- `flen = 0.1`: length of the final beam segment in case of no intersection [m]
- `r_res::Int = 10`: radial (angular) resolution per beamlet
- `z_res::Int = 8`: longitudinal resolution per beamlet segment

# Makie kwargs

- `color = :red`
- `transparency = true`
"""
function live_render!(
        axis::_RenderEnv,
        bg::BMO.AstigmaticBeamGroup;
        # kwargs
        render_every::Int = 5,
        flen = 0.1,
        r_res::Int = 10,
        z_res::Int = 8,
        # Makie kwargs
        color = :red,
        transparency = true,
        kwargs...
    )
    mesh_obs = Observable{Any}(_astigmatic_group_mesh(bg; flen, r_res, z_res, render_every))
    mesh_plot = mesh!(axis, mesh_obs; color, transparency, kwargs...)
    return AstigmaticGroupRenderHandle(bg, axis, mesh_plot, mesh_obs, Float64(flen), r_res, z_res, render_every)
end

function update_render!(h::AstigmaticGroupRenderHandle)
    h.mesh_obs[] = _astigmatic_group_mesh(h.thing; flen = h.flen, r_res = h.r_res, z_res = h.z_res,
        render_every = h.render_every)
    return h
end

remove_render!(h::AstigmaticGroupRenderHandle) = (delete!(h.axis, h.mesh_plot); return nothing)
