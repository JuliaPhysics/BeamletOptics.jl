"""
    render!(ax, agb::AstigmaticGaussianBeamlet; kwargs...)

Render the 1/e² envelope of the [`AstigmaticGaussianBeamlet`](@ref) as a smooth 3D surface.

With `show_beams = true` the generating rays are overlayed into the axis as follows:

- `chief` beam: red
- `divergence` beam: green
- `waist` beam: blue

# Keyword args

- `show_beams = false`: plot the generating rays (chief, waist, divergence)
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

- `color = :red`
- `transparency = true`
"""
function render!(
        axis::_RenderEnv,
        agb::BMO.AstigmaticGaussianBeamlet{T};
        # kwargs
        show_beams = false,
        show_pos = false,
        r_res = 64,
        z_res = 100,
        flen = 0.1,
        show_waist = false,
        # Polarization kwargs
        show_polarization = false,
        pol_λ = nothing,
        pol_scale = 1.0,
        pol_focus_exponent = 1.0,
        pol_gain_max = 10.0,
        pol_ppl = 32,
        pol_color = :crimson,
        pol_linewidth = 2.0,
        # Makie kwargs
        color = :red,
        transparency = true,
        markersize = 10,
        kwargs...
) where {T}
    vs = LinRange(0, 2π, r_res)
    for child in PreOrderDFS(agb)
        # Length tracking variable
        p = child.parent
        if isnothing(p)
            l = zero(T)
        else
            l = length(p)
        end

        for ray in BMO.rays(child.c)
            # Calculate local segment length
            if isnothing(BMO.intersection(ray))
                l_local = flen
            else
                l_local = length(ray)
            end
            us = LinRange(0, l_local, z_res) .+ l
            # Precompute waist parameters at each z
            params = [BMO.waist_parameters(child, u) for u in us]

            # sort params b and c
            ps = getindex.(params, 1)
            bs = getindex.(params, 2)
            cs = getindex.(params, 3)

            # Ensure elliptical basis vectors (bs, cs) vary smoothly along the segment
            # to prevent mesh twisting/flips, especially when passing through a focus.
            for j in 2:length(bs)
                b_new = bs[j]
                c_new = cs[j]
                b_old = bs[j - 1]
                c_old = cs[j - 1]

                # preserve handedness (prevent inside-out mesh flips)
                cross_old = cross(b_old, c_old)
                cross_new = cross(b_new, c_new)
                if dot(cross_new, cross_old) < 0
                    b_new = -b_new
                end

                # optimal rotation to align with previous slice (Parallel Transport)
                # MUST normalize to prevent magnitude-bias from twisting the mesh when major/minor axes swap!
                bn = normalize(b_new)
                cn = normalize(c_new)
                bo = normalize(b_old)
                co = normalize(c_old)

                X = dot(bn, bo) + dot(cn, co)
                Y = dot(cn, bo) - dot(bn, co)
                phi = atan(Y, X)

                cos_phi = cos(phi)
                sin_phi = sin(phi)

                bs[j] = b_new * cos_phi + c_new * sin_phi
                cs[j] = -b_new * sin_phi + c_new * cos_phi
            end

            # Build surface mesh matrices
            pts = Matrix{Point3{T}}(undef, length(params), length(vs))

            for i in eachindex(params)
                pts[i, :] = [BMO.ellipse(v, ps[i], bs[i], cs[i]) for v in vs]
            end

            Xt = getindex.(pts, 1)
            Yt = getindex.(pts, 2)
            Zt = getindex.(pts, 3)

            # Render the envelope as a smooth surface
            surface!(axis, Xt, Yt, Zt;
                color = fill(color, size(Xt)),
                transparency,
                kwargs...
            )

            # Optionally, plot waist ellipse
            if show_waist
                scatter!.(axis, pts; color, markersize)
            end

            # Bump length tracker
            if !isnothing(BMO.intersection(ray))
                l += length(ray)
            end
        end

        # Optionally, plot generating rays
        if show_beams
            render!(axis, child.c; show_pos, flen, color = :red)
            render!(axis, child.dxp; show_pos, flen, color = :green)
            render!(axis, child.dyp; show_pos, flen, color = :green)
            render!(axis, child.dxm; show_pos, flen, color = :green)
            render!(axis, child.dym; show_pos, flen, color = :green)
            render!(axis, child.wxp; show_pos, flen, color = :blue)
            render!(axis, child.wyp; show_pos, flen, color = :blue)
            render!(axis, child.wxm; show_pos, flen, color = :blue)
            render!(axis, child.wym; show_pos, flen, color = :blue)
        end
    end

    if show_polarization
        pts = _polarization_points(agb; flen, λ_vis = pol_λ, scale = pol_scale,
            focus_exponent = pol_focus_exponent, gain_max = pol_gain_max, ppl = pol_ppl)
        _render_field_curve!(axis, pts; color = pol_color, linewidth = pol_linewidth)
    end

    return axis
end
