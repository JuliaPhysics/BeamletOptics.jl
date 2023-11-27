"""
    render_object!(axis, mesh::AbstractMesh)

Render `mesh` into the specified 3D-`axis`.
"""
render_object!(::Any, mesh::AbstractMesh) = nothing

render_object!(::Any, ::AbstractEntity) = nothing
render_object!(axis, object::AbstractObject) = render_object!(axis, shape(object))

"""
    render_system!(axis, system::AbstractSystem)

Render all objects contained in the `system`.
"""
function render_system!(axis, system::AbstractSystem)
    # Avoid use of objects(system)
    for object in system.objects
        render_object!(axis, object)
    end
    return nothing
end

"""
    render_ray!(axis, ray::AbstractRay; color=:blue, flen=1.0)

Renders a `ray` as a 3D line. If the ray has no intersection, the substitute length `flen` is used.
"""
function render_ray!(axis, ray::AbstractRay; color = :blue, flen = 1.0)
    if isnothing(intersection(ray))
        len = flen
    else
        len = length(ray.intersection)
    end
    temp = ray.pos + len * ray.dir

    _render_ray!(axis, ray, temp; color)

    return nothing
end
_render_ray!(::Any, ::AbstractRay, ::AbstractVector; color = :blue) = nothing

"""
    render_beam!(axis, beam::Beam; color=:blue, flen=1.0)

Render the entire `beam` into the specified 3D-`axis`. A `color` can be specified.
"""
function render_beam!(axis, beam::Beam; color = :blue, flen = 1.0)
    for child in PreOrderDFS(beam)
        for ray in rays(child)
            render_ray!(axis, ray, color = color, flen = flen)
        end
    end
    return nothing
end

"""
    render_beam!(axis, gauss::GaussianBeamlet; color=:red, flen=0.1, show_beams=false)

Render the surface of a `GaussianBeamlet` as `color`. With `show_beams` the generating rays are plotted as follows:

- `chief` ray: red
- `divergence` ray: green
- `waist` ray: blue
"""
function render_beam!(axis,
        gauss::GaussianBeamlet{T};
        color = :red,
        flen = 0.1,
        show_beams = false) where {T}
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
                u = LinRange(0, flen, 50)
            else
                u = LinRange(0, length(ray), 50)
            end
            v = LinRange(0, 2π, 20)
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
            render_gaussian_beam_surface!(axis, Xt, Yt, Zt; transparency = true, colormap = [color, color])
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
render_gaussian_beam_surface!(::Any, X, Y, Z; kwargs...) = nothing

"""
    render_object_normals!(axis, mesh::Mesh; l=0.01)

Helper function to visualize the `mesh` normals of an object. The normals are displayed with a standard `l`ength of 0.01
"""
function render_object_normals!(axis, mesh::Mesh; l = 0.01)
    for fID in 1:size(mesh.faces)[1]
        nml = normal3d(mesh, fID)
        pos = @view mesh.vertices[mesh.faces[fID, 1], :]
        vec = pos + l * nml
        _render_object_normal!(axis, pos, vec; color=blue)
    end
    return nothing
end
_render_object_normal!(::Any, ::AbstractVector, ::AbstractVector; color = :blue) = nothing
