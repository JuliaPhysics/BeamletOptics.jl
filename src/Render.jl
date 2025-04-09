function render!(::Any, ::Any, kwargs...)
    @warn "It appears no suitable Makie backend is loaded to trigger BeamletOpticsMakieExt"
    return nothing
end

"""
    render_ray!(axis, ray::AbstractRay; color=:blue, flen=1.0)

Renders a `ray` as a 3D line. If the ray has no intersection, the substitute length `flen` is used.
"""
function render_ray!(axis, ray::AbstractRay; color = :blue, flen = 1.0, show_pos=false)
    if isnothing(intersection(ray))
        len = flen
    else
        len = length(ray.intersection)
    end
    temp = ray.pos + len * ray.dir

    _render_ray!(axis, ray, temp; color, show_pos)

    return nothing
end
_render_ray!(::Any, ::AbstractRay, ::AbstractVector; color = :blue) = nothing

"""
    render_beam!(axis, beam::Beam; color=:blue, flen=1.0, show_pos=false)

Render the entire `beam` into the specified 3D-`axis`. A `color` can be specified.
"""
function render_beam!(axis, beam::Beam; color = :blue, flen = 1.0, show_pos=false)
    for child in PreOrderDFS(beam)
        for ray in rays(child)
            render_ray!(axis, ray; color, flen, show_pos)
        end
    end
    return nothing
end

_render_beam!(::Any, ::AbstractBeam; color=:red, flen=1.0) = nothing

"""
    render_beam!(axis, gauss::GaussianBeamlet; color=:red, flen=0.1, show_beams=false)

Render the surface of a `GaussianBeamlet` as `color`. With `show_beams` the generating rays are plotted as follows:

- `chief` ray: red
- `divergence` ray: green
- `waist` ray: blue

# Keyword args

- `color = :red`: color of the beam as per the Makie syntax, i.e. :blue
- `flen = 0.1`: length of the final beam in case of no intersection
- `show_beams = false`: plot the generating rays of the [`GaussianBeamlet`](@ref)
- `r_res::Int = 50`: radial resolution of the beam
- `z_res::Int = 100`: resolution along the optical axis of the beam
"""
function render_beam!(axis,
        gauss::GaussianBeamlet{T};
        color = :red,
        flen = 0.1,
        show_beams = false,
        r_res = 50,
        z_res = 100) where {T}
    for child in PreOrderDFS(gauss)
        # Length tracking variable
        p = AbstractTrees.parent(child)
        if isnothing(p)
            l = zero(T)
        else
            l = length(p)
        end
        # Render each ray segment
        for ray in rays(child.chief)
            # Generate local u, v coords
            if isnothing(intersection(ray))
                u = LinRange(0, flen, z_res)
            else
                u = LinRange(0, length(ray), z_res)
            end
            v = LinRange(0, 2π, r_res)
            # Calculate beam surface at origin along y-axis
            w = gauss_parameters(child, u .+ l)[1]
            X = [w[i] * cos(v) for (i, u) in enumerate(u), v in v]
            Y = [u for u in u, v in v]
            Z = [w[i] * sin(v) for (i, u) in enumerate(u), v in v]
            # Transform into world coords
            R = align3d([0, 1, 0], ray.dir)
            Xt = R[1, 1] * X + R[1, 2] * Y + R[1, 3] * Z .+ ray.pos[1]
            Yt = R[2, 1] * X + R[2, 2] * Y + R[2, 3] * Z .+ ray.pos[2]
            Zt = R[3, 1] * X + R[3, 2] * Y + R[3, 3] * Z .+ ray.pos[3]
            render_surface!(axis, Xt, Yt, Zt; transparency = true, colormap = [color, color])
            # Bump length tracker
            if !isnothing(intersection(ray))
                l += length(ray)
            end
        end
        # Optionally, plot generating rays
        if show_beams
            render_beam!(axis, child.chief, flen = flen, color = :red)
            render_beam!(axis, child.divergence, flen = flen, color = :green)
            render_beam!(axis, child.waist, flen = flen, color = :blue)
        end
    end
    return nothing
end
render_surface!(::Any, X, Y, Z; kwargs...) = nothing

"""
    render_object_normals!(axis, mesh::Mesh; l=0.01)

Helper function to visualize the `mesh` normals of an object. The normals are displayed with a standard `l`ength of 0.01
"""
function render_object_normals!(axis, mesh::Mesh; l = 0.01)
    for fID in 1:size(mesh.faces)[1]
        nml = normal3d(mesh, fID)
        pos = @view mesh.vertices[mesh.faces[fID, 1], :]
        vec = pos + l * nml
        _render_object_normal!(axis, pos, vec; color=:blue)
    end
    return nothing
end
_render_object_normal!(::Any, ::AbstractVector, ::AbstractVector; color = :blue) = nothing

