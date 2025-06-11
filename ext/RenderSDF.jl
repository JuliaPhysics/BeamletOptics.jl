"""
    render!(ax, sdf; kwargs...)

Render the surface of the `sdf` based on the marching cubes algorithm into the specified `axis`.

Additional kwargs can be passed into the mesh plot.
"""
function render!(
        ax::_RenderEnv,
        sdf::BMO.AbstractSDF;
        # kwargs
        x_resolution::Int=100,
        y_resolution::Int=100,
        z_resolution::Int=100,
        kwargs...
    )
    # Get object limits
    xmin, xmax, ymin, ymax, zmin, zmax = BMO.bounding_box(sdf)

    x = LinRange(xmin - 1e-4, xmax + 1e-4, x_resolution)
    y = LinRange(ymin - 1e-4, ymax + 1e-4, y_resolution)
    z = LinRange(zmin - 1e-4, zmax + 1e-4, z_resolution)
    sdf_values = Float32.([BMO.sdf(sdf, [i, j, k]) for i in x, j in y, k in z])
    mc = MC(sdf_values; x = Float32.(x), y = Float32.(y), z = Float32.(z))
    march(mc)
    vertices = transpose(reinterpret(reshape, Float32, mc.vertices))
    faces = transpose(reinterpret(reshape, Int64, mc.triangles))
    mesh!(ax, vertices, faces; kwargs...)
    return nothing
end

function render!(ax::_RenderEnv, s::BMO.UnionSDF; kwargs...)
    for sdf in s.sdfs
        render!(ax, sdf; kwargs...)
    end
end
