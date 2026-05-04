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
        # Makie kwargs
        color = :red,
        transparency = true,
        kwargs...
    ) where {T}
    
    for child in PreOrderDFS(agb)
        l = length(child) + flen
        
        # Longitudinal and angular ranges
        us = LinRange(0, l, z_res)
        vs = LinRange(0, 2π, r_res)
        
        # Precompute waist parameters at each z
        params = [BMO.waist_parameters(child, u) for u in us]
        
        # Build surface mesh matrices
        pts = [BMO.ellipse(v, p[1], p[2], p[3]) for p in params, v in vs]

        Xt = getindex.(pts, 1)
        Yt = getindex.(pts, 2)
        Zt = getindex.(pts, 3)

        # Render the envelope as a smooth surface
        surface!(axis, Xt, Yt, Zt; 
            color = fill(color, size(Xt)), 
            transparency, 
            kwargs...)

        # Optionally, plot generating rays
        if show_beams
            render!(axis, child.c; show_pos, flen,   color = :red)
            render!(axis, child.dxp; show_pos, flen, color = :green)
            render!(axis, child.dyp; show_pos, flen, color = :green)
            render!(axis, child.dxm; show_pos, flen, color = :green)
            render!(axis, child.dym; show_pos, flen, color = :green)
            render!(axis, child.wxp; show_pos, flen, color = :blue)
            render!(axis, child.wyp; show_pos, flen, color = :blue)
            render!(axis, child.wxm; show_pos, flen, color = :blue)
            render!(axis, child.wym; show_pos, flen, color = :blue)
        end

        # Optionally, plot waist ellipse
        if show_waist
            scatter!.(axis, pts; color)
        end
    end
    return axis
end
