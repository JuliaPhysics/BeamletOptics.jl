#=
Tiles of the beam catalog on `docs/src/basics/beams/overview.md`, one per entry of
`docs/src/components/beam_catalog.json`. Beams and beam groups are traced through a small
test system of a bi-convex lens and a detector near its focus. Wrapped in a module since all
showcase scripts are included into the same namespace.
=#
module BeamCatalogTiles

include(joinpath(@__DIR__, "..", "catalog_assets", "catalog_render.jl"))

const mm = 1e-3
const inch = BMO.inch
const λ = 1e-6

"LB1811-like bi-convex lens and a detector near its focus."
function test_system()
    lens = SphericalLens(34.9mm, -34.9mm, 6.8mm, inch, λ -> 1.5168)
    translate3d!(lens, [0, 50mm, 0])
    detector = Detector(inch)
    translate3d!(detector, [0, 85mm, 0])
    return System([lens, detector])
end

"Tile that traces the beam returned by `build` through the test system and renders both."
function traced_tile(build; trace_kwargs = (;), kwargs...)
    return function (ax)
        system = test_system()
        beam = build()
        solve_system!(system, beam; trace_kwargs...)
        render!(ax, system)
        render!(ax, beam; kwargs...)
    end
end

# nearly side-on, so that focusing along the optical axis (+y) is visible. The scenes are long
# along y, so the fitted bounding sphere leaves room to zoom in on the landscape tile.
const BEAM_VIEW = normalize([1.0, -0.25, 0.35])

# kwargs for single beamlets, rays and groups of rays
const beamlet_kwargs = (color = (:red, 0.4), flen = 20mm)
const ray_group_kwargs = (render_every = 1, linewidth = 1.5, flen = 20mm)
# kwargs for groups of astigmatic Gaussian beamlets: every beamlet, coarse translucent surfaces
const beamlet_group_kwargs = (render_every = 1, r_res = 16, z_res = 40, color = (:red, 0.15), flen = 10mm)

"Three astigmatic Gaussian beamlets side by side."
function astigmatic_group()
    beams = [AstigmaticGaussianBeamlet([x, 0, 0], [0, 1, 0], λ, 1mm) for x in (-6mm, 0, 6mm)]
    return AstigmaticBeamGroup(beams, [0, 0, 0], [0, 1, 0])
end

"Gaussian amplitude with a flat phase, sampled on a 5×5 grid."
function wavefront_decomposition()
    x = y = collect(range(-6mm, 6mm; length = 5))
    amplitude = [exp(-(xi^2 + yi^2) / (4mm)^2) for xi in x, yi in y]
    phase = zeros(length(x), length(y))
    return WavefrontBeamletDecomposition(x, y, amplitude, phase, [0, 1, 0], λ)
end

const TILES = [
    # rays
    "Ray" => (ax -> render!(ax, Ray([0, 0, 0], [0, 1, 0]); flen = 60mm, show_pos = true, linewidth = 2)),
    "PolarizedRay" => (ax -> render!(ax, PolarizedRay([0, 0, 0], [0, 1, 0], λ, [1, 0, im] / √2);
        flen = 60mm, linewidth = 2, show_polarization = true, pol_λ = 6mm)),
    # beams
    "Beam" => traced_tile(() -> Beam([0, 0, 5mm], [0, 1, 0], λ); flen = 20mm, linewidth = 2),
    # beamlets
    "GaussianBeamlet" => traced_tile(() -> GaussianBeamlet([0, 0, 0], [0, 1, 0], λ, 2mm); beamlet_kwargs...),
    "AstigmaticGaussianBeamlet" => traced_tile(
        () -> AstigmaticGaussianBeamlet([0, 0, 0], [0, 1, 0], λ, 2mm, 1mm; support = [0, 0, 1]);
        trace_kwargs = (; check_invariant = false), beamlet_kwargs...),
    # beam groups
    "CollimatedSource" => traced_tile(
        () -> CollimatedSource([0, 0, 0], [0, 1, 0], 15mm, λ; num_rings = 3, num_rays = 60);
        ray_group_kwargs...),
    "UniformDiscSource" => traced_tile(
        () -> UniformDiscSource([0, 0, 0], [0, 1, 0], 15mm, λ; num_rays = 60);
        ray_group_kwargs...),
    "PointSource" => traced_tile(
        () -> PointSource([0, 0, 0], [0, 1, 0], deg2rad(8), λ; num_rings = 3, num_rays = 60);
        ray_group_kwargs...),
    "UniformPointSource" => traced_tile(
        () -> UniformPointSource([0, 0, 0], [0, 1, 0], deg2rad(8), λ; num_rays = 60);
        ray_group_kwargs...),
    "AstigmaticBeamGroup" => traced_tile(astigmatic_group; beamlet_group_kwargs...),
    "CollimatedGaussianBeamletSource" => traced_tile(
        () -> CollimatedGaussianBeamletSource([0, 0, 0], [0, 1, 0], 12mm, λ, 3mm; n_grid = 4);
        beamlet_group_kwargs...),
    "GaussianBeamletDecomposition" => traced_tile(
        () -> GaussianBeamletDecomposition([0, 0, 0], [0, 1, 0], λ, 4mm; n_grid = 5);
        beamlet_group_kwargs...),
    "SphericalGaussianBeamletSource" => traced_tile(
        () -> SphericalGaussianBeamletSource([0, 0, 0], [0, 1, 0], deg2rad(8), λ; num_rings = 2, num_rays = 40);
        beamlet_group_kwargs...),
    "EllipticalGaussianBeamletSource" => traced_tile(
        () -> EllipticalGaussianBeamletSource([0, 0, 0], [0, 1, 0], deg2rad(10), deg2rad(4), λ;
            num_rings = 2, num_rays = 40);
        beamlet_group_kwargs...),
    "WavefrontBeamletDecomposition" => traced_tile(wavefront_decomposition; beamlet_group_kwargs...),
]

const CATALOG_JSON = joinpath(@__DIR__, "..", "..", "components", "beam_catalog.json")

let
    check_tiles(CATALOG_JSON, first.(TILES))
    for (name, draw!) in TILES
        save_tile(name, draw!; view = BEAM_VIEW, zoom = 1.3)
    end
end

end # module
