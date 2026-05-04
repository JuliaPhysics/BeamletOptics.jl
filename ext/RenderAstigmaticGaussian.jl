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

            b0 = first(bs)
            c0 = first(cs)

            # flip ellipse basis vectors if necessary
            bs[findall(dot.(bs, b0) .< 0)] .*= -1
            cs[findall(dot.(cs, c0) .< 0)] .*= -1

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
                scatter!.(axis, pts; color)
            end

            # Bump length tracker
            if !isnothing(BMO.intersection(ray))
                l += length(ray)
            end
        end

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


    end
    return axis
end
