#=
Render gallery of all exported components and beam types.

Renders every catalogue entry into an isometric and a side view, records render time and mesh
statistics and writes `output/index.md` plus `output/overview.png`. See `README.md`.

Run from the repository root with:

    julia --project=docs gallery/render_gallery.jl
=#
using GLMakie, BeamletOptics
using LinearAlgebra: normalize, norm
import CairoMakie # overview only: GLMakie clips figures larger than the screen

const BMO = BeamletOptics

GLMakie.activate!(visible = false)

const mm = 1e-3
const inch = BMO.inch

const GALLERY_DIR = @__DIR__
const OUTPUT_DIR = joinpath(GALLERY_DIR, "output")
const REPO_DIR = normpath(joinpath(GALLERY_DIR, ".."))
const ASSET_DIR = joinpath(REPO_DIR, "docs", "src", "assets")

# refractive indices
const NBK7 = λ -> 1.5168
const λs = [488e-9, 707e-9, 1064e-9]
const NLAK22 = DiscreteRefractiveIndex(λs, [1.6591, 1.6456, 1.6374])
const NSF10 = DiscreteRefractiveIndex(λs, [1.7460, 1.7168, 1.7021])

#=
S1 catalogue
Each entry stores the exported names it covers explicitly, since `nameof(typeof(obj))` does not
resolve constructor functions (e.g. `RoundPlanoMirror` returns a `Mirror`).
=#

"Small test system for the beam entries: LB1811-like bi-convex lens and a detector near the focus."
function beam_test_system()
    lens = SphericalLens(34.9mm, -34.9mm, 6.8mm, inch, NBK7)
    translate3d!(lens, [0, 50mm, 0])
    detector = Detector(inch)
    translate3d!(detector, [0, 85mm, 0])
    return System([lens, detector])
end

function traced(beam; kwargs...)
    system = beam_test_system()
    solve_system!(system, beam; kwargs...)
    return (system, beam)
end

"Two cemented spherical lenses (AC254-100-like), assembled via the `DoubletLens` type directly."
function manual_doublet()
    front = SphericalLens(62.8mm, -45.7mm, 4mm, inch, NLAK22)
    back = SphericalLens(-45.7mm, -128.2mm, 2.5mm, inch, NSF10)
    translate3d!(back, [0, thickness(BMO.shape(front)), 0])
    return DoubletLens(front, back)
end

"Three cemented lenses built from surfaces, assembled via the `TripletLens` type directly."
function manual_triplet()
    front = Lens(SphericalSurface(60mm, inch), SphericalSurface(25mm, inch), 3mm, NSF10)
    middle = Lens(SphericalSurface(25mm, inch), SphericalSurface(-25mm, inch), 8mm, NLAK22)
    back = Lens(SphericalSurface(-25mm, inch), SphericalSurface(-60mm, inch), 3mm, NSF10)
    translate3d!(middle, [0, thickness(BMO.shape(front)), 0])
    translate3d!(back, [0, thickness(BMO.shape(front)) + thickness(BMO.shape(middle)), 0])
    return TripletLens(front, middle, back)
end

"Thorlabs 354710-C-like even asphere (showcase `aspherical_lens_showcase.jl`)."
function asphere()
    return Lens(
        EvenAsphericalSurface(3.618e-3, 3.04e-3, -44.874,
            [0, -0.14756 * (1e3)^3, 0.035194 * (1e3)^5, -0.0032262 * (1e3)^7,
                0.0018592 * (1e3)^9, 0.00036658 * (1e3)^11, -0.00016039 * (1e3)^13,
                -3.1846e-5 * (1e3)^15]),
        EvenAsphericalSurface(2.161e-3, 3.7e-3, -10.719,
            [0, -0.096568 * (1e3)^3, 0.026771 * (1e3)^5, -0.011261 * (1e3)^7,
                0.0019879 * (1e3)^9, 0.00015579 * (1e3)^11, -0.00012433 * (1e3)^13,
                1.5264e-5 * (1e3)^15]),
        0.7e-3,
        n -> 1.580200
    )
end

"Acylindrical lens (showcase `cylindrical_lens_showcase.jl`)."
function acylinder()
    return Lens(
        AcylindricalSurface(-15.538e-3, 25e-3, 50e-3, -1.0,
            [0, 1.1926075e-5 * (1e3)^3, -2.9323497e-9 * (1e3)^5, -1.8718889e-11 * (1e3)^7,
                -1.7009961e-14 * (1e3)^9, 3.5481542e-17 * (1e3)^11, 6.5241296e-20 * (1e3)^13]),
        7.5e-3,
        n -> 1.777
    )
end

function cube_group()
    cbs = CubeBeamsplitter(inch, NBK7)
    mount = MeshDummy(joinpath(ASSET_DIR, "bs_assets", "CBS Mount.stl"))
    return ObjectGroup([cbs, mount])
end

const CATALOGUE = [
    # mirrors
    (category = "mirrors", name = "Mirror", exports = (:Mirror,),
        build = () -> Mirror(BMO.CylinderSDF(inch / 2, 6mm))),
    (category = "mirrors", name = "SquarePlanoMirror2D", exports = (:SquarePlanoMirror2D,),
        build = () -> SquarePlanoMirror2D(inch)),
    (category = "mirrors", name = "RectangularPlanoMirror", exports = (:RectangularPlanoMirror,),
        build = () -> RectangularPlanoMirror(50mm, 25mm, 6mm)),
    (category = "mirrors", name = "SquarePlanoMirror", exports = (:SquarePlanoMirror,),
        build = () -> SquarePlanoMirror(inch, 6mm)),
    (category = "mirrors", name = "RoundPlanoMirror", exports = (:RoundPlanoMirror,),
        build = () -> RoundPlanoMirror(inch, 6mm)),
    (category = "mirrors", name = "SphericalMirror", exports = (:SphericalMirror,),
        build = () -> SphericalMirror(200mm, 6mm, inch)),
    (category = "mirrors", name = "RightAnglePrismMirror", exports = (:RightAnglePrismMirror,),
        build = () -> RightAnglePrismMirror(inch, inch)),
    (category = "mirrors", name = "ConicMirror", exports = (:ConicMirror,),
        build = () -> ConicMirror(200mm, -0.5, 50mm; thickness = 6mm)),
    (category = "mirrors", name = "OffAxisConicMirror", exports = (:OffAxisConicMirror,),
        build = () -> OffAxisConicMirror(200mm, -0.5, 40mm, inch)),
    (category = "mirrors", name = "ParabolicMirror", exports = (:ParabolicMirror,),
        build = () -> ParabolicMirror(200mm, 100mm; hole_diameter = 10mm, thickness = 5mm)),
    (category = "mirrors", name = "OffAxisParabolicMirror", exports = (:OffAxisParabolicMirror,),
        build = () -> OffAxisParabolicMirror(2inch, inch)),
    (category = "mirrors", name = "EllipsoidalMirror", exports = (:EllipsoidalMirror,),
        build = () -> EllipsoidalMirror(40mm, 160mm, 100mm; thickness = 5mm)),
    (category = "mirrors", name = "OffAxisEllipsoidalMirror", exports = (:OffAxisEllipsoidalMirror,),
        build = () -> OffAxisEllipsoidalMirror(100mm, 300mm, 30mm, inch)),
    (category = "mirrors", name = "HyperbolicMirror", exports = (:HyperbolicMirror,),
        build = () -> HyperbolicMirror(-50mm, 200mm, 30mm; thickness = 5mm)),
    (category = "mirrors", name = "OffAxisHyperbolicMirror", exports = (:OffAxisHyperbolicMirror,),
        build = () -> OffAxisHyperbolicMirror(-50mm, 200mm, 20mm, inch)),
    # lenses
    (category = "lenses", name = "Lens", exports = (:Lens,),
        build = () -> Lens(BMO.CylinderSDF(inch / 2, 5mm), NBK7)),
    (category = "lenses", name = "SphericalLens", exports = (:SphericalLens,),
        build = () -> SphericalLens(34.9mm, -34.9mm, 6.8mm, inch, NBK7)),
    (category = "lenses", name = "ThinLens", exports = (:ThinLens,),
        build = () -> ThinLens(34.9mm, -34.9mm, inch, 1.5)),
    (category = "lenses", name = "SphericalDoubletLens", exports = (:SphericalDoubletLens,),
        build = () -> SphericalDoubletLens(87.9mm, -105.6mm, -1000, 6mm, 3mm, inch, NLAK22, NSF10)),
    (category = "lenses", name = "DoubletLens", exports = (:DoubletLens,),
        build = manual_doublet),
    (category = "lenses", name = "SphericalTripletLens", exports = (:SphericalTripletLens,),
        build = () -> SphericalTripletLens(60mm, 25mm, -25mm, -60mm, 3mm, 8mm, 3mm, inch,
            NSF10, NLAK22, NSF10)),
    (category = "lenses", name = "TripletLens", exports = (:TripletLens,),
        build = manual_triplet),
    # surfaces / lens constructors
    (category = "surfaces", name = "Lens_SphericalSurface_CircularFlatSurface",
        exports = (:Lens, :SphericalSurface, :CircularFlatSurface),
        build = () -> Lens(SphericalSurface(25.8mm, inch), CircularFlatSurface(inch), 5.3mm, NBK7)),
    (category = "surfaces", name = "Lens_CircularFlatSurface_window",
        exports = (:Lens, :CircularFlatSurface),
        build = () -> Lens(CircularFlatSurface(inch), CircularFlatSurface(inch), 3mm, NBK7)),
    (category = "surfaces", name = "Lens_EvenAsphericalSurface",
        exports = (:Lens, :EvenAsphericalSurface),
        build = asphere),
    (category = "surfaces", name = "Lens_CylindricalSurface",
        exports = (:Lens, :CylindricalSurface, :RectangularFlatSurface),
        build = () -> Lens(CylindricalSurface(5.2mm, 10mm, 20mm), 5.9mm, NBK7)),
    (category = "surfaces", name = "Lens_AcylindricalSurface",
        exports = (:Lens, :AcylindricalSurface, :RectangularFlatSurface),
        build = acylinder),
    # prisms
    (category = "prisms", name = "Prism", exports = (:Prism,),
        build = () -> Prism(BMO.BoxSDF(20mm, 20mm, 20mm), NBK7)),
    (category = "prisms", name = "RightAnglePrism", exports = (:RightAnglePrism,),
        build = () -> RightAnglePrism(inch, inch, NBK7)),
    # splitters
    (category = "splitters", name = "ThinBeamsplitter", exports = (:ThinBeamsplitter,),
        build = () -> ThinBeamsplitter(inch)),
    (category = "splitters", name = "RoundThinBeamsplitter", exports = (:RoundThinBeamsplitter,),
        build = () -> RoundThinBeamsplitter(inch)),
    (category = "splitters", name = "RectangularPlateBeamsplitter",
        exports = (:RectangularPlateBeamsplitter,),
        build = () -> RectangularPlateBeamsplitter(36mm, 25mm, 5mm, NBK7)),
    (category = "splitters", name = "RoundPlateBeamsplitter", exports = (:RoundPlateBeamsplitter,),
        build = () -> RoundPlateBeamsplitter(inch, 3mm, NBK7)),
    (category = "splitters", name = "CubeBeamsplitter", exports = (:CubeBeamsplitter,),
        build = () -> CubeBeamsplitter(inch, NBK7)),
    (category = "splitters", name = "RectangularCompensatorPlate",
        exports = (:RectangularCompensatorPlate,),
        build = () -> RectangularCompensatorPlate(36mm, 25mm, 5mm, NBK7)),
    # polarizers
    (category = "polarizers", name = "PolarizationFilter", exports = (:PolarizationFilter,),
        build = () -> PolarizationFilter(inch)),
    (category = "polarizers", name = "RoundPolarizationFilter", exports = (:RoundPolarizationFilter,),
        build = () -> RoundPolarizationFilter(inch)),
    (category = "polarizers", name = "RoundLinearPolarizer",
        exports = (:RoundLinearPolarizer, :LinearPolarizer),
        build = () -> RoundLinearPolarizer(inch, 1mm, 1mm, NBK7)),
    # detectors
    (category = "detectors", name = "Detector", exports = (:Detector,),
        build = () -> Detector(inch)),
    # dummies
    (category = "dummies", name = "MeshDummy", exports = (:MeshDummy,),
        build = () -> MeshDummy(joinpath(ASSET_DIR, "Benchy.stl"))),
    (category = "dummies", name = "NonInteractableObject", exports = (:NonInteractableObject,),
        build = () -> NonInteractableObject(BMO.CylinderSDF(6mm, 50mm))),
    (category = "dummies", name = "IntersectableObject", exports = (:IntersectableObject,),
        build = () -> IntersectableObject(BMO.BoxSDF(20mm, 10mm, 20mm))),
    # misc
    (category = "misc", name = "Retroreflector", exports = (:Retroreflector,),
        build = () -> Retroreflector(inch)),
    # groups
    (category = "groups", name = "ObjectGroup", exports = (:ObjectGroup,),
        build = cube_group),
    # beams, traced through a small lens system
    (category = "beams", name = "Beam", exports = (:Beam,),
        build = () -> traced(Beam([0, 0, 5mm], [0, 1, 0], 1e-6))),
    (category = "beams", name = "CollimatedSource", exports = (:CollimatedSource,),
        build = () -> traced(CollimatedSource([0, 0, 0], [0, 1, 0], 15mm; num_rings = 3, num_rays = 60))),
    (category = "beams", name = "GaussianBeamlet", exports = (:GaussianBeamlet,),
        build = () -> traced(GaussianBeamlet([0, 0, 0], [0, 1, 0], 1e-6, 2mm))),
    (category = "beams", name = "AstigmaticGaussianBeamlet", exports = (:AstigmaticGaussianBeamlet,),
        build = () -> traced(AstigmaticGaussianBeamlet([0, 0, 0], [0, 1, 0], 1e-6, 2mm, 1mm;
                support = [0, 0, 1]); check_invariant = false)),
]

#=
S2 coverage check
=#

const COVERAGE_SECTIONS = ("mirrors", "lenses", "prisms", "detectors", "splitters",
    "polarizing components", "dummies")

"Returns true if `name` is an exported component constructor (type or constructor function)."
function is_component_constructor(name::Symbol)
    isdefined(BMO, name) || return false
    isuppercase(first(string(name))) || return false
    f = getfield(BMO, name)
    f isa Type && return f <: BMO.AbstractObject
    return f isa Function
end

"Parses `src/Exports.jl` and returns the exported component constructor names of the relevant sections."
function exported_components()
    names = Symbol[]
    section = ""
    for line in eachline(joinpath(REPO_DIR, "src", "Exports.jl"))
        stripped = strip(line)
        m = match(r"^#\s*(.*)$", stripped)
        if !isnothing(m)
            section = strip(m.captures[1])
            continue
        end
        isempty(stripped) && continue
        for id in eachmatch(r"[A-Za-z_][A-Za-z0-9_!]*", replace(stripped, r"^export\s+" => ""))
            name = Symbol(id.match)
            if section in COVERAGE_SECTIONS || (section == "misc" && name == :Retroreflector)
                is_component_constructor(name) && push!(names, name)
            end
        end
    end
    return unique(names)
end

function coverage_check()
    covered = Set(n for entry in CATALOGUE for n in entry.exports)
    missing_names = [n for n in exported_components() if !(n in covered)]
    println("missing: ", isempty(missing_names) ? "none" : join(missing_names, ", "))
    return missing_names
end

#=
S3/S4 rendering and metrics
=#

const MARCHING_CUBES_METHOD = which(render!, Tuple{LScene, BMO.AbstractSDF})

"Returns true if any shape of `x` is rendered by the generic marching cubes `render!(::_RenderEnv, ::AbstractSDF)`."
uses_fallback(x) = false
uses_fallback(x::Union{Tuple, AbstractVector}) = any(uses_fallback, x)
uses_fallback(x::BMO.AbstractSystem) = any(uses_fallback, x.objects)
uses_fallback(x::BMO.AbstractObject) = uses_fallback(BMO.shape(x))
uses_fallback(x::BMO.UnionSDF) = any(uses_fallback, x.sdfs)
uses_fallback(x::BMO.AbstractSDF) = which(render!, Tuple{LScene, typeof(x)}) === MARCHING_CUBES_METHOD

render_item!(ax, item) = render!(ax, item)
render_item!(ax, items::Tuple) = foreach(item -> render!(ax, item), items)

"Returns (vertices, faces) of mesh-like plots, recursing through child plots."
function mesh_stats(p)
    nv, nf = 0, 0
    if p isa Makie.Mesh
        m = p[1][]
        if m isa Makie.GeometryBasics.Mesh
            nv += length(Makie.GeometryBasics.coordinates(m))
            nf += length(Makie.GeometryBasics.faces(m))
        end
    elseif p isa Makie.Surface
        z = p[3][]
        nv += length(z)
        nf += 2 * max(size(z, 1) - 1, 0) * max(size(z, 2) - 1, 0)
    end
    for child in p.plots
        cv, cf = mesh_stats(child)
        nv += cv
        nf += cf
    end
    return nv, nf
end

"Bounding box center and size of the given plots."
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

const VIEWS = (iso = normalize([1, -1, 1]), side = [1.0, 0.0, 0.0])

first_line(e) = first(split(sprint(showerror, e), '\n'))

"Renders one view of `entry`, returns a NamedTuple with the metrics."
function render_view(entry, view::Symbol, path::String)
    fig = Figure(size = (600, 450))
    ax = LScene(fig[1, 1]; show_axis = false)
    obj = entry.build()
    n0 = length(ax.scene.plots)
    t = @elapsed render_item!(ax, obj)
    new_plots = ax.scene.plots[(n0 + 1):end]
    nv, nf = 0, 0
    for p in new_plots
        pv, pf = mesh_stats(p)
        nv += pv
        nf += pf
    end
    c, d = plots_bbox(new_plots)
    set_view(ax, c .+ d .* VIEWS[view], c, [0, 0, 1])
    save(path, fig)
    return (time = t, plots = length(new_plots), vertices = nv, faces = nf,
        fallback = uses_fallback(obj))
end

function render_entry(entry)
    dir = joinpath(OUTPUT_DIR, entry.category)
    mkpath(dir)
    result = Dict{Symbol, Any}(:time => NaN, :plots => missing, :vertices => missing,
        :faces => missing, :fallback => missing, :error => "", :images => Dict{Symbol, String}())
    # first call compiles, it is not recorded
    try
        render_item!(LScene(Figure()[1, 1]), entry.build())
    catch e
        result[:error] = first_line(e)
        return result
    end
    for view in keys(VIEWS)
        file = "$(entry.name)_$(view).png"
        try
            m = render_view(entry, view, joinpath(dir, file))
            result[:images][view] = "$(entry.category)/$(file)"
            if view == :iso
                result[:time] = m.time
                result[:plots] = m.plots
                result[:vertices] = m.vertices
                result[:faces] = m.faces
                result[:fallback] = m.fallback
            end
        catch e
            isempty(result[:error]) && (result[:error] = "$(view): $(first_line(e))")
        end
    end
    return result
end

#=
S5 output
=#

md_escape(s) = replace(string(s), "|" => "\\|", "\n" => " ")
fmt_bool(b) = ismissing(b) ? "?" : (b ? "**yes**" : "no")
fmt_time(t) = isnan(t) ? "-" : string(round(1e3 * t; digits = 1))

function img_cell(result, view)
    haskey(result[:images], view) || return "-"
    return "<img src=\"$(result[:images][view])\" width=\"250\">"
end

function write_index(results, missing_names, runtime)
    open(joinpath(OUTPUT_DIR, "index.md"), "w") do io
        println(io, "# BeamletOptics render gallery\n")
        println(io, "Entries: $(length(results)), total runtime: $(round(runtime; digits = 1)) s, ",
            "missing exports: ", isempty(missing_names) ? "none" : join(missing_names, ", "), "\n")
        for category in unique(entry.category for (entry, _) in results)
            println(io, "## $(category)\n")
            println(io, "| name | iso | side | time [ms] | plots | vertices | faces | fallback | error |")
            println(io, "|:--|:--|:--|--:|--:|--:|--:|:--|:--|")
            for (entry, r) in results
                entry.category == category || continue
                println(io, "| $(entry.name) | $(img_cell(r, :iso)) | $(img_cell(r, :side)) | ",
                    "$(fmt_time(r[:time])) | $(coalesce(r[:plots], "-")) | ",
                    "$(coalesce(r[:vertices], "-")) | $(coalesce(r[:faces], "-")) | ",
                    "$(fmt_bool(r[:fallback])) | $(md_escape(r[:error])) |")
            end
            println(io)
        end
    end
    return nothing
end

function write_overview(results; ncols = 6)
    fig = Figure()
    for (i, (entry, r)) in enumerate(results)
        row, col = fldmod1(i, ncols)
        ax = Axis(fig[row, col]; title = entry.name, titlesize = 12, width = 280, height = 210)
        hidedecorations!(ax)
        hidespines!(ax)
        if haskey(r[:images], :iso)
            img = Makie.FileIO.load(joinpath(OUTPUT_DIR, r[:images][:iso]))
            image!(ax, rotr90(img))
        else
            text!(ax, 0.5, 0.5; text = "failed", align = (:center, :center), color = :red)
            limits!(ax, 0, 1, 0, 1)
        end
    end
    resize_to_layout!(fig)
    save(joinpath(OUTPUT_DIR, "overview.png"), fig; backend = CairoMakie)
    return nothing
end

function main()
    rm(OUTPUT_DIR; recursive = true, force = true)
    mkpath(OUTPUT_DIR)
    missing_names = coverage_check()
    results = []
    runtime = @elapsed for entry in CATALOGUE
        print("rendering $(entry.category)/$(entry.name) ... ")
        r = render_entry(entry)
        println(isempty(r[:error]) ? "$(fmt_time(r[:time])) ms" : "ERROR: $(r[:error])")
        push!(results, (entry, r))
    end
    write_index(results, missing_names, runtime)
    write_overview(results)
    n_images = count(isfile, (joinpath(OUTPUT_DIR, f) for (_, r) in results for f in values(r[:images])))
    println("done: $(length(results)) entries, $(n_images) images, $(round(runtime; digits = 1)) s")
    return results
end

main()
