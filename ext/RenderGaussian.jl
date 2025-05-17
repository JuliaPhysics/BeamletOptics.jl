"""
    render!(ax, gauss::GaussianBeamlet; kwargs...)

Render the 1/e² envelope of the `GaussianBeamlet` into the specified `axis`.

With `show_beams = true` the generating rays are overlayed into the axis as follows:

- `chief` beam: red
- `divergence` beam: green
- `waist` beam: blue

# Keyword args

- `show_beams = false`: plot the generating rays of the [`GaussianBeamlet`](@ref)
- `flen = 0.1`: length of the final beam in case of no intersection
- `r_res::Int = 50`: radial resolution of the beam
- `z_res::Int = 100`: resolution along the optical axis of the beam

# Makie kwargs

- `color = :red`
- `transparency = true`

Additional kwargs can be passed into the surface plot of the Gaussian envelope.
"""
function render!(
        axis::_RenderEnv,
        gauss::GaussianBeamlet{T};
        # kwargs
        show_beams = false,
        show_pos = false,
        r_res = 50,
        z_res = 100,
        flen = 0.1,
        # Makie kwargs
        color = :red,
        transparency = true,
        kwargs...
    ) where {T}
    for child in PreOrderDFS(gauss)
        # Length tracking variable
        p = child.parent # FIXME: `AbstractTrees` not defined error
        if isnothing(p)
            l = zero(T)
        else
            l = length(p)
        end
        # Render each ray segment
        for ray in BMO.rays(child.chief)
            # Generate local u, v coords
            if isnothing(BMO.intersection(ray))
                u = LinRange(0, flen, z_res)
            else
                u = LinRange(0, length(ray), z_res)
            end
            v = LinRange(0, 2π, r_res)
            # Calculate beam surface at origin along y-axis
            w = BMO.gauss_parameters(child, u .+ l)[1]
            X = [w[i] * cos(v) for (i, u) in enumerate(u), v in v]
            Y = [u for u in u, v in v]
            Z = [w[i] * sin(v) for (i, u) in enumerate(u), v in v]
            # Transform into world coords
            R = BMO.align3d([0, 1, 0], ray.dir)
            Xt = R[1, 1] * X + R[1, 2] * Y + R[1, 3] * Z .+ BMO.position(ray)[1]
            Yt = R[2, 1] * X + R[2, 2] * Y + R[2, 3] * Z .+ BMO.position(ray)[2]
            Zt = R[3, 1] * X + R[3, 2] * Y + R[3, 3] * Z .+ BMO.position(ray)[3]

            surface!(axis, Xt, Yt, Zt; transparency, colormap = [color, color], kwargs...)

            # Bump length tracker
            if !isnothing(BMO.intersection(ray))
                l += length(ray)
            end
        end
        # Optionally, plot generating rays
        if show_beams
            render!(axis, child.chief; show_pos, transparency, flen, color = :red)
            render!(axis, child.divergence; show_pos, transparency, flen, color = :green)
            render!(axis, child.waist; show_pos, transparency, flen, color = :blue)
        end
    end
    return nothing
end