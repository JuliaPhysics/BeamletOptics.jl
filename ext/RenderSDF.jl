"""
    render_object!(axis, s::AbstractSDF)

Render the surface of `s` based on the marching cubes algorithm.
"""
function render!(ax::_RenderEnv, s::BMO.AbstractSDF; kwargs...)
    # Get object limits
    xmin, xmax, ymin, ymax, zmin, zmax = BMO.bounding_box(s)

    x = LinRange(xmin - 1e-4, xmax + 1e-4, 100)
    y = LinRange(ymin - 1e-4, ymax + 1e-4, 100)
    z = LinRange(zmin - 1e-4, zmax + 1e-4, 100)
    sdf_values = Float32.([BMO.sdf(s, [i, j, k]) for i in x, j in y, k in z])
    mc = MC(sdf_values; x = Float32.(x), y = Float32.(y), z = Float32.(z))
    march(mc)
    vertices = transpose(reinterpret(reshape, Float32, mc.vertices))
    faces = transpose(reinterpret(reshape, Int64, mc.triangles))
    mesh!(ax, vertices, faces; kwargs...)
    return nothing
end