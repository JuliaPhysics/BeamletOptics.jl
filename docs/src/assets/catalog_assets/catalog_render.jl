#=
Shared tile renderer for the catalogs on the component and beam overview pages (see
`docs/src/components/ComponentCatalog.vue`). Included by `catalog_showcase.jl` and
`beam_catalog_showcase.jl`, each within its own module.

Every tile is a 400×300 `LScene` seen from the same isometric direction. The camera is framed
such that the bounding sphere of the rendered plots fits into the field of view.
=#
using GLMakie, BeamletOptics
using LinearAlgebra: normalize, norm

const BMO = BeamletOptics

GLMakie.activate!(; ssao = true)

"""
Oblique viewing direction from `+x`, `-y`, `+z`. Flatter than isometric, so that the rims and
stacked elements of optics along the `y`-axis remain visible.
"""
tile_view() = normalize([1.0, -0.45, 0.5])

"Bounding box center and diagonal of the given plots."
function plots_bbox(plots)
    lo = fill(Inf, 3)
    hi = fill(-Inf, 3)
    for p in plots
        bb = Makie.data_limits(p)
        o = Makie.origin(bb)
        w = Makie.widths(bb)
        all(isfinite, o) && all(isfinite, w) || continue
        lo .= min.(lo, o)
        hi .= max.(hi, o .+ w)
    end
    all(isfinite, lo) || return zeros(3), 1.0
    return (lo .+ hi) ./ 2, max(norm(hi .- lo), 1e-3)
end

"""
    save_tile(name, draw!; view = tile_view(), zoom = 1.0)

Calls `draw!(ax)` on a fresh `LScene` and saves the figure as `catalog_<name>.png`. `zoom > 1`
moves the camera closer than the fitted bounding sphere distance. A failing tile is logged and
skipped, the catalog then shows the entry without an image.
"""
function save_tile(name::AbstractString, draw!::Function; view = tile_view(), zoom::Real = 1.0)
    try
        fig = Figure(size = (400, 300))
        ax = LScene(fig[1, 1]; show_axis = false)
        n0 = length(ax.scene.plots)
        draw!(ax)
        c, d = plots_bbox(ax.scene.plots[(n0 + 1):end])
        # the automatic centering on display would reset the eye distance to d, which cuts off
        # compact objects like cubes
        cam = ax.scene.camera_controls
        cam.settings[:center] = false
        set_view(ax, c .+ (d / 2) / sind(cam.fov[] / 2) / zoom .* view, c, [0, 0, 1])
        save("catalog_$(name).png", fig; px_per_unit = 2)
    catch e
        @warn "Catalog tile $(name) failed" exception = (e, catch_backtrace())
    end
    return nothing
end

"Names of the entries in the catalog JSON file at `path`."
catalog_names(path) = [m.captures[1] for m in eachmatch(r"\"name\":\s*\"([^\"]+)\"", read(path, String))]

"Warns about catalog entries without a tile and tiles without a catalog entry."
function check_tiles(json_path, tile_names)
    isfile(json_path) || return @warn "Catalog file not found" json_path
    names = catalog_names(json_path)
    missing_tiles = setdiff(names, tile_names)
    missing_entries = setdiff(tile_names, names)
    isempty(missing_tiles) || @warn "Catalog entries without a tile in $(basename(json_path))" missing_tiles
    isempty(missing_entries) || @warn "Tiles without a catalog entry in $(basename(json_path))" missing_entries
    return nothing
end
