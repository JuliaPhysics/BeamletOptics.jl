"""
    render_object!(axis, mesh::AbstractMesh)

Render `mesh` into the specified 3D-`axis`.
"""
function render_object!(axis, mesh::AbstractMesh)
    mesh!(axis, vertices(mesh), faces(mesh), transparency=true)
    return nothing
end

render_object!(axis, ::AbstractEntity) = nothing
render_object!(axis, object::AbstractObject) = render_object!(axis, shape(object))

"""
    render_system!(axis::Axis3, system::System)

Render all objects contained in the `system`.
"""
function render_system!(axis, system::System)
    for object in system.object
        render_object!(axis, object)
    end
    return nothing
end

"""
    render_ray!(axis, ray::AbstractRay; color=:blue, flen=1.0)

Renders a `ray` as a 3D line. If the ray has no intersection, the substitute length `flen` is used.
"""
function render_ray!(axis, ray::AbstractRay; color=:blue, flen=1.0)
    if isnothing(intersection(ray))
        len = flen
    else
        len = ray.intersection.t
    end
    temp = ray.pos + len * ray.dir
    lines!(
        axis,
        [ray.pos[1], temp[1]],
        [ray.pos[2], temp[2]],
        [ray.pos[3], temp[3]],
        color=color
    )
    return nothing
end

"""
    render_beam!(axis, beam::Beam; color=:blue, flen=1.0)

Render the entire `beam` into the specified 3D-`axis`. A `color` can be specified.
"""
function render_beam!(axis, beam::Beam; color=:blue, flen=1.0)
    for ray in beam.rays
        render_ray!(axis, ray, color=color, flen=flen)
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
function render_beam!(axis, gauss::GaussianBeamlet{T}; color=:red, flen=0.1, show_beams=false) where T
    # Length tracking variable
    l::T = 0
    for ray in gauss.chief.rays
        # Generate local u, v coords
        if isnothing(intersection(ray))
            u = LinRange(0 , flen, 50)
        else
            u = LinRange(0, length(ray), 50)
        end
        v = LinRange(0, 2π, 20)
        # Calculate beam surface at origin along y-axis
        w = gauss_parameters(gauss, u .+ l)[1]
        X = [w[i] * cos(v) for (i, u) in enumerate(u), v in v]  
        Y = [u for u in u, v in v]                              
        Z = [w[i] * sin(v) for (i, u) in enumerate(u), v in v]
        # Transform into world coords
        R = align3d([0,1,0], ray.dir)
        Xt = R[1,1]*X + R[1,2]*Y + R[1,3]*Z .+ ray.pos[1]
        Yt = R[2,1]*X + R[2,2]*Y + R[2,3]*Z .+ ray.pos[2]
        Zt = R[3,1]*X + R[3,2]*Y + R[3,3]*Z .+ ray.pos[3]
        surface!(axis, Xt, Yt, Zt, transparency = true, colormap = [color, color])
        # Bump length tracker
        if !isnothing(intersection(ray)) 
            l += length(ray)
        end
    end
    # Optionally, plot generating rays
    if show_beams
        render_beam!(axis, gauss.chief, flen=flen, color=:red)
        render_beam!(axis, gauss.divergence, flen=flen, color=:green)
        render_beam!(axis, gauss.waist, flen=flen, color=:blue)
    end
    return nothing
end

"""
    render_object_normals!(axis::Axis3, mesh::Mesh; l=0.01)

Helper function to visualize the `mesh` normals of an object. The normals are displayed with a standard `l`ength of 0.01
"""
function render_object_normals!(axis, mesh::Mesh; l=0.01)
    for fID = 1:size(mesh.faces)[1]
        nml = normal3d(mesh, fID)
        pos = mesh.vertices[mesh.faces[fID,1],:]
        vec = pos + l * nml
        lines!(
            axis,
            [pos[1], vec[1]],
            [pos[2], vec[2]],
            [pos[3], vec[3]],
            color=:blue
        )
    end
    return nothing
end
