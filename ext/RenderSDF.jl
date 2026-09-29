"""
    render!(ax, sdf; kwargs...)

Render the surface of the `sdf` into the specified `axis`. Shapes with an analytic tessellation
(see `_tessellate`) are rendered as a single mesh. Other shapes are rendered based on the marching
cubes algorithm, with a resolution of `x_resolution`, `y_resolution` and `z_resolution` samples
per axis (defaults: a voxel size of 1/100 of the largest bounding box width, at most 150 samples
per axis).

Additional kwargs can be passed into the mesh plot.
"""
function render!(
        ax::_RenderEnv,
        sdf::BMO.AbstractSDF;
        # kwargs
        x_resolution::Union{Int, Nothing}=nothing,
        y_resolution::Union{Int, Nothing}=nothing,
        z_resolution::Union{Int, Nothing}=nothing,
        kwargs...
    )
    _has_mesh(sdf) && return _render_mesh!(ax, sdf; kwargs...)
    # Get object limits, padded relative to the object size
    xmin, xmax, ymin, ymax, zmin, zmax = BMO.bounding_box(sdf)
    limits = ((xmin, xmax), (ymin, ymax), (zmin, zmax))
    wmax = maximum(hi - lo for (lo, hi) in limits)
    pad = 0.02 * wmax
    voxel = wmax / 100
    function samples((lo, hi), n)
        n = something(n, clamp(ceil(Int, (hi - lo + 2pad) / voxel) + 1, 2, 150))
        return LinRange(lo - pad, hi + pad, n)
    end
    x, y, z = map(samples, limits, (x_resolution, y_resolution, z_resolution))
    sdf_values = Float32.([BMO.sdf(sdf, [i, j, k]) for i in x, j in y, k in z])
    mc = MC(sdf_values; x = Float32.(x), y = Float32.(y), z = Float32.(z))
    march(mc)
    if isempty(mc.vertices)
        return nothing
    end
    pts = [Point3f(v...) for v in mc.vertices]
    fcs = [GLTriangleFace(t...) for t in mc.triangles]
    normals = [let n = BMO.normal3d(sdf, Point3(p...))
                   any(isnan, n) ? Vec3f(0, 1, 0) : Vec3f(n...)
               end for p in pts]
    gb_mesh = Mesh(pts, fcs; normal = normals)
    mesh!(ax, gb_mesh; kwargs...)
    return nothing
end

"""SDF types without a specific `render!` method that are rendered from their analytic mesh."""
const _MeshedSDF = Union{_PrimitiveSDF, BMO.ConvexCylinderSDF, BMO.ConcaveCylinderSDF,
    BMO.MeniscusLensSDF, BMO.UnionSDF, BMO.DifferenceSDF}

# Separate method, such that the dispatch shows which shapes are not rendered by marching cubes
function render!(
        ax::_RenderEnv,
        sdf::_MeshedSDF;
        x_resolution::Union{Int, Nothing}=nothing,
        y_resolution::Union{Int, Nothing}=nothing,
        z_resolution::Union{Int, Nothing}=nothing,
        kwargs...
    )
    _has_mesh(sdf) && return _render_mesh!(ax, sdf; kwargs...)
    return invoke(render!, Tuple{_RenderEnv, BMO.AbstractSDF}, ax, sdf;
        x_resolution, y_resolution, z_resolution, kwargs...)
end
