"""
    render!(ax, agb::AstigmaticGaussianBeamlet; kwargs...)

Render the 1/e² envelope of the [`AstigmaticGaussianBeamlet`](@ref) as a smooth 3D surface.

# Keyword args

- `show_beams = false`: plot the generating rays (chief, waist, divergence)
- `flen = 0.1`: length of the final beam segment in case of no intersection [m]
- `z_res::Int = 100`: longitudinal resolution
- `r_res::Int = 64`: radial (angular) resolution

# Makie kwargs

- `color = :red`
- `alpha = 0.3`: transparency of the envelope
"""
function render!(
        axis::_RenderEnv,
        agb::BMO.AstigmaticGaussianBeamlet{T};
        # kwargs
        show_beams = false,
        show_pos = false,
        flen = 0.1,
        z_res = 100,
        r_res = 64,
        # Makie kwargs
        color = :red,
        alpha = 0.3,
        kwargs...
    ) where {T}
    
    # Helper to calculate a point on the elliptical cross-section
    ellipse(t, p0, w1, w2) = p0 + w1 * cos(t) + w2 * sin(t)

    for child in PreOrderDFS(agb)
        l = length(child) + flen
        
        # Longitudinal and angular ranges
        us = LinRange(0, l, z_res)
        vs = LinRange(0, 2π, r_res)
        
        # Precompute waist parameters at each z
        params = [BMO.waist_parameters(child, u) for u in us]
        
        # Build surface mesh matrices
        Xt = [ellipse(v, p[1], p[2], p[3])[1] for p in params, v in vs]
        Yt = [ellipse(v, p[1], p[2], p[3])[2] for p in params, v in vs]
        Zt = [ellipse(v, p[1], p[2], p[3])[3] for p in params, v in vs]

        # Render the envelope as a smooth surface
        surface!(axis, Xt, Yt, Zt; 
            color = fill(color, size(Xt)), 
            alpha = alpha, 
            transparency = true, 
            kwargs...)

        # Optionally, plot generating rays
        if show_beams
            render!(axis, child.c; show_pos, flen, color = :red)
            render!(axis, child.dxp; show_pos, flen, color = :green)
            render!(axis, child.wxp; show_pos, flen, color = :blue)
            render!(axis, child.dyp; show_pos, flen, color = :cyan)
            render!(axis, child.wyp; show_pos, flen, color = :magenta)
        end
    end
    return axis
end
