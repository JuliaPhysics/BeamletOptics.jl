import GeometryBasics

#=
Geometry and plot function of the `AstigmaticGaussianBeamlet` and the `AstigmaticBeamGroup`, see
`RenderGaussian.jl` for the parts shared with the `GaussianBeamlet`
=#

"""
    _smooth_astigmatic_axes!(bs, cs)

Applies parallel transport to the elliptical basis vectors `bs`, `cs` (one pair per
longitudinal sample) so that they vary smoothly along a segment, preserving handedness
and preventing mesh twisting/flips near a focus, see the `render!` method of
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

        # Align with the previous slice, avoids twisting when the axes swap. The axes must be
        # normalized, otherwise their magnitudes bias the rotation.
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
    vs = _ring_angles(r_res)
    for child in PreOrderDFS(agb)
        p = child.parent
        l = isnothing(p) ? zero(T) : length(p)
        for ray in BMO.rays(child.c)
            l_local = isnothing(BMO.intersection(ray)) ? flen : length(ray)
            us = LinRange(0, l_local, z_res) .+ l
            p0s, bs, cs, = BMO.waist_parameters(child, us)
            # The rings must run counterclockwise about the ray, such that the faces point
            # outwards, see `_push_grid_faces!`. The smoothing keeps this handedness.
            if dot(cross(bs[1], cs[1]), BMO.direction(ray)) < 0
                bs[1] = -bs[1]
            end
            _smooth_astigmatic_axes!(bs, cs)
            offset = length(pts)
            for i in eachindex(p0s), v in vs
                push!(pts, Point3f(BMO.ellipse(v, p0s[i], bs[i], cs[i])))
            end
            _push_grid_faces!(faces, offset, r_res, z_res)
            if !isnothing(BMO.intersection(ray))
                l += length(ray)
            end
        end
    end
    return nothing
end

"""
    _gaussian_mesh(x::Union{AstigmaticGaussianBeamlet, AstigmaticBeamGroup}; flen, r_res, z_res, render_every = 1)

Returns a single mesh of the 1/e² envelopes of all chief ray segments of the beamlet `x`, or of
every `render_every`-th beamlet of the group `x`, see the `render!` method of
[`BMO.AstigmaticGaussianBeamlet`](@ref).
"""
function _gaussian_mesh(x::Union{BMO.AstigmaticGaussianBeamlet, BMO.AstigmaticBeamGroup};
        flen, r_res::Int, z_res::Int, render_every::Int = 1)
    pts = Point3f[]
    faces = GLTriangleFace[]
    for agb in _beamlets(x, render_every)
        _append_astigmatic_mesh!(pts, faces, agb; flen, r_res, z_res)
    end
    return Mesh(pts, faces)
end

"""
    _plot_astigmatic!(axis, x; kwargs...) -> _BeamPlots

The plot function of the [`BMO.AstigmaticGaussianBeamlet`](@ref) and of the
[`BMO.AstigmaticBeamGroup`](@ref) (every `render_every`-th beamlet): the envelopes of all segments
in a single `mesh` plot, the rings of `show_waist` in a single `scatter` plot, the overlay of
`show_beams` and the polarization curves of `show_polarization` in a single `lines` plot. See
`render!(axis, agb)` for the kwargs.
"""
function _plot_astigmatic!(
        axis::_RenderEnv,
        x::Union{BMO.AstigmaticGaussianBeamlet, BMO.AstigmaticBeamGroup};
        show_beams,
        show_pos,
        r_res::Int,
        z_res::Int,
        flen,
        render_every::Int,
        show_waist,
        show_polarization,
        pol_λ,
        pol_scale,
        pol_focus_exponent,
        pol_gain_max,
        pol_ppl,
        pol_color,
        pol_linewidth,
        color,
        transparency,
        markersize,
        kwargs...
    )
    bp = _BeamPlots()
    coloring = _coloring(color)
    envelope = _envelope!(bp, coloring, () -> _gaussian_mesh(x; flen, r_res, z_res, render_every),
        () -> _envelope_colors(x, coloring; r_res, z_res, render_every))
    _coupled!(bp, mesh!, axis, envelope, coloring; transparency, kwargs...)
    # the vertices of the envelope, updated with it
    show_waist && _coupled!(bp, scatter!, axis, envelope, coloring, m -> Vector{Point3f}(GeometryBasics.coordinates(m));
        markersize)
    show_beams && _plot_generating_beams!(bp, axis, x; flen, render_every, show_pos, transparency)
    if show_polarization
        curve = _observe!(bp,
            () -> _field_curve(x; flen, render_every, λ_vis = pol_λ, scale = pol_scale,
                focus_exponent = pol_focus_exponent, gain_max = pol_gain_max, ppl = pol_ppl))
        _add!(bp, lines!(axis, curve; color = pol_color, linewidth = pol_linewidth))
    end
    return bp
end

"""The kwargs of `render!` and `live_render!` of astigmatic beamlets except `r_res`, `z_res`."""
_astigmatic_defaults(kwargs) = (;
    show_beams = false,
    show_pos = false,
    flen = 0.1,
    render_every = 5,
    show_waist = false,
    show_polarization = false,
    pol_λ = nothing,
    pol_scale = 1.0,
    pol_focus_exponent = 1.0,
    pol_gain_max = 10.0,
    pol_ppl = 32,
    pol_color = :crimson,
    pol_linewidth = 2.0,
    color = :red,
    transparency = true,
    markersize = 10,
    kwargs...)

"""
    render!(ax, agb::AstigmaticGaussianBeamlet; kwargs...)

Render the 1/e² envelope of the [`AstigmaticGaussianBeamlet`](@ref) as a smooth 3D surface, all
segments as a single mesh.

With `show_beams = true` the generating rays are overlayed into the axis as follows:

- `chief` beam: red
- `divergence` beam: green
- `waist` beam: blue

# Keyword args

- `show_beams = false`: plot the generating rays (chief, waist, divergence)
- `show_pos = false`: with `show_beams`, marks the positions of the generating rays
- `show_waist = false`: marks the vertices of the envelope, i.e. the ellipses along the beam
- `flen = 0.1`: length of the final beam segment in case of no intersection [m]
- `z_res::Int = 100`: longitudinal resolution
- `r_res::Int = 64`: radial (angular) resolution

# Polarization kwargs

- `show_polarization = false`: overlay the E-field curve along the chief ray
- `pol_λ = nothing`: visualization wavelength [m], default = total plotted length / 20.
  This is a plotting parameter, not the physical ray wavelength; the curve shows the
  `t = 0` snapshot `Re{E⊥·exp(i·k·s)}` along the accumulated optical path `s`. Gouy
  phase and phase-front curvature are ignored; the curve is a qualitative
  visualization only. Values finer than `plotted length / 2000` are clamped with a
  warning.
- `pol_scale = 1.0`: curve amplitude as a multiple of the mean 1/e² beam radius at the
  start of the beamlet (at the maximum |E⊥|)
- `pol_focus_exponent = 1.0`: the amplitude follows the on-axis field amplitude,
  `(A_ref / A(z))^(pol_focus_exponent/2)` with `A = wx·wy`. `1` is the physical scaling
  `E ∝ √(w0x·w0y / (wx·wy))`, which raises the curve in the focus; values in `(0, 1)`
  compress the gain for tight foci, `0` gives a constant amplitude.
- `pol_gain_max = 10.0`: upper bound on the focus gain, in multiples of the reference
  amplitude at the start of the beamlet. The curve saturates at this value rather than
  leaving the scene when the beam focuses tighter downstream than its input waist.
- `pol_ppl = 32`: sample points per `pol_λ` along the curve
- `pol_color = :crimson`: field curve color
- `pol_linewidth = 2.0`: field curve line width

# Makie kwargs

- `color = :red`: envelope color. `color = :wavelength` (or `(:wavelength, alpha)`) colors each
  beamlet segment by its wavelength, see [`wavelength_color`](@ref). The overlay of `show_beams`
  keeps its colors.
- `transparency = true`
- `markersize = 10`: size of the points of `show_waist`

Additional kwargs are passed into the mesh plot of the envelope.
"""
function render!(axis::_RenderEnv, agb::BMO.AstigmaticGaussianBeamlet; r_res::Int = 64, z_res::Int = 100, kwargs...)
    _plot_astigmatic!(axis, agb; _astigmatic_defaults(kwargs)..., r_res, z_res)
    return axis
end

"""
    render!(ax, bg::AstigmaticBeamGroup; kwargs...)

Renders every `render_every`-th beamlet of the [`BMO.AstigmaticBeamGroup`](@ref) like
`render!(ax, agb::AstigmaticGaussianBeamlet)`; the envelopes of all rendered beamlets are merged
into a single mesh.

# Keyword args

- `render_every::Int = 5`: renders only every e.g. fifth beamlet of the group
- `r_res::Int = 64`, `z_res::Int = 100`: resolution per beamlet

All other kwargs of the [`BMO.AstigmaticGaussianBeamlet`](@ref) method apply with the same
defaults.
"""
function render!(axis::_RenderEnv, bg::BMO.AstigmaticBeamGroup; r_res::Int = 64, z_res::Int = 100, kwargs...)
    _plot_astigmatic!(axis, bg; _astigmatic_defaults(kwargs)..., r_res, z_res)
    return nothing
end

"""
    live_render!(axis, agb::AstigmaticGaussianBeamlet; kwargs...)

Live-renders the 1/e² envelope of the [`BMO.AstigmaticGaussianBeamlet`](@ref) like `render!`, as a
single mesh, see [`live_render!`](@ref). The kwargs and their defaults are those of `render!`.
"""
function live_render!(axis::_RenderEnv, agb::BMO.AstigmaticGaussianBeamlet; r_res::Int = 64, z_res::Int = 100, kwargs...)
    kw = _astigmatic_defaults(kwargs)
    bp = _plot_astigmatic!(axis, agb; kw..., r_res, z_res)
    return BeamRenderHandle(agb, axis, bp, (; flen = Float64(kw.flen), render_every = 1, r_res, z_res))
end

"""
    live_render!(axis, bg::AstigmaticBeamGroup; kwargs...)

Live-renders every `render_every`-th beamlet envelope of the [`BMO.AstigmaticBeamGroup`](@ref)
like `render!`, as a single merged mesh, see [`live_render!`](@ref). The default resolution is
deliberately low (compared to `render!`) since a group typically contains many beamlets.

# Keyword args

- `render_every::Int = 5`: renders only every e.g. fifth beamlet of the group
- `flen = 0.1`: length of the final beam segment in case of no intersection [m]
- `r_res::Int = 10`: radial (angular) resolution per beamlet
- `z_res::Int = 8`: longitudinal resolution per beamlet segment

All other kwargs of `render!` apply with the same defaults.
"""
function live_render!(axis::_RenderEnv, bg::BMO.AstigmaticBeamGroup; r_res::Int = 10, z_res::Int = 8, kwargs...)
    kw = _astigmatic_defaults(kwargs)
    bp = _plot_astigmatic!(axis, bg; kw..., r_res, z_res)
    return BeamRenderHandle(bg, axis, bp, (; flen = Float64(kw.flen), kw.render_every, r_res, z_res))
end
