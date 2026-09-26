#=
Tiles of the component catalog on `docs/src/basics/components/components.md`, one per entry of
`docs/src/components/catalog.json`. Wrapped in a module since all showcase scripts are included
into the same namespace.
=#
module ComponentCatalogTiles

include(joinpath(@__DIR__, "catalog_render.jl"))

const mm = 1e-3
const inch = BMO.inch

# refractive indices, the exact values do not matter for a catalog tile
const n_crown = λ -> 1.5168
const n_flint = λ -> 1.7

"Two cemented spherical lenses (AC254-100-like), assembled via the `DoubletLens` type directly."
function manual_doublet()
    front = SphericalLens(62.8mm, -45.7mm, 4mm, inch, n_crown)
    back = SphericalLens(-45.7mm, -128.2mm, 2.5mm, inch, n_flint)
    translate3d!(back, [0, thickness(BMO.shape(front)), 0])
    return DoubletLens(front, back)
end

"Three cemented lenses built from surfaces, assembled via the `TripletLens` type directly."
function manual_triplet()
    front = Lens(SphericalSurface(60mm, inch), SphericalSurface(25mm, inch), 3mm, n_flint)
    middle = Lens(SphericalSurface(25mm, inch), SphericalSurface(-25mm, inch), 8mm, n_crown)
    back = Lens(SphericalSurface(-25mm, inch), SphericalSurface(-60mm, inch), 3mm, n_flint)
    translate3d!(middle, [0, thickness(BMO.shape(front)), 0])
    translate3d!(back, [0, thickness(BMO.shape(front)) + thickness(BMO.shape(middle)), 0])
    return TripletLens(front, middle, back)
end

"Tile that renders the object returned by `build`."
object_tile(build; kwargs...) = ax -> render!(ax, build(); kwargs...)

"Rotates `obj` by 180° about `z`, so that the tile shows the reflective face of off-axis mirrors
and the facets of the retroreflector."
flipped(obj) = (zrotate3d!(obj, π); obj)

const TILES = [
    # mirrors
    "Mirror" => object_tile(() -> Mirror(BMO.CylinderSDF(inch / 2, 6mm))),
    "SquarePlanoMirror2D" => object_tile(() -> SquarePlanoMirror2D(inch)),
    "SquarePlanoMirror" => object_tile(() -> SquarePlanoMirror(inch, 6mm)),
    "RectangularPlanoMirror" => object_tile(() -> RectangularPlanoMirror(50mm, 25mm, 6mm)),
    "RoundPlanoMirror" => object_tile(() -> RoundPlanoMirror(inch, 6mm)),
    "SphericalMirror" => object_tile(() -> SphericalMirror(200mm, 6mm, inch)),
    "RightAnglePrismMirror" => object_tile(() -> RightAnglePrismMirror(inch, inch)),
    "ConicMirror" => object_tile(() -> ConicMirror(200mm, -0.5, 50mm; thickness = 6mm)),
    "OffAxisConicMirror" => object_tile(() -> flipped(OffAxisConicMirror(200mm, -0.5, 40mm, inch))),
    "ParabolicMirror" => object_tile(() -> ParabolicMirror(200mm, 100mm; hole_diameter = 10mm, thickness = 5mm)),
    "OffAxisParabolicMirror" => object_tile(() -> flipped(OffAxisParabolicMirror(2inch, inch))),
    "EllipsoidalMirror" => object_tile(() -> EllipsoidalMirror(40mm, 160mm, 100mm; thickness = 5mm)),
    "OffAxisEllipsoidalMirror" => object_tile(() -> flipped(OffAxisEllipsoidalMirror(100mm, 300mm, 30mm, inch))),
    "HyperbolicMirror" => object_tile(() -> HyperbolicMirror(-50mm, 200mm, 30mm; thickness = 5mm)),
    "OffAxisHyperbolicMirror" => object_tile(() -> flipped(OffAxisHyperbolicMirror(-50mm, 200mm, 20mm, inch))),
    "Retroreflector" => object_tile(() -> flipped(Retroreflector(inch))),
    # lenses
    "Lens" => object_tile(() -> Lens(SphericalSurface(25.8mm, inch), CircularFlatSurface(inch), 5.3mm, n_crown)),
    "SphericalLens" => object_tile(() -> SphericalLens(34.9mm, -34.9mm, 6.8mm, inch, n_crown)),
    "ThinLens" => object_tile(() -> ThinLens(34.9mm, -34.9mm, inch, 1.5)),
    "DoubletLens" => object_tile(manual_doublet),
    "SphericalDoubletLens" => object_tile(() -> SphericalDoubletLens(87.9mm, -105.6mm, -1000, 6mm, 3mm, inch, n_crown, n_flint)),
    "TripletLens" => object_tile(manual_triplet),
    "SphericalTripletLens" => object_tile(() -> SphericalTripletLens(60mm, 25mm, -25mm, -60mm, 3mm, 8mm, 3mm, inch,
        n_flint, n_crown, n_flint)),
    # prisms
    "Prism" => object_tile(() -> Prism(BMO.BoxSDF(20mm, 20mm, 20mm), n_crown)),
    "RightAnglePrism" => object_tile(() -> RightAnglePrism(inch, inch, n_crown)),
    # beamsplitters
    "ThinBeamsplitter" => object_tile(() -> ThinBeamsplitter(inch)),
    "RoundThinBeamsplitter" => object_tile(() -> RoundThinBeamsplitter(inch)),
    "RectangularPlateBeamsplitter" => object_tile(() -> RectangularPlateBeamsplitter(36mm, 25mm, 5mm, n_crown)),
    "RoundPlateBeamsplitter" => object_tile(() -> RoundPlateBeamsplitter(inch, 3mm, n_crown)),
    "CubeBeamsplitter" => object_tile(() -> CubeBeamsplitter(inch, n_crown)),
    "RectangularCompensatorPlate" => object_tile(() -> RectangularCompensatorPlate(36mm, 25mm, 5mm, n_crown)),
    # detectors
    "Detector" => object_tile(() -> Detector(inch)),
    # polarizers
    "PolarizationFilter" => object_tile(() -> PolarizationFilter(inch)),
    "RoundPolarizationFilter" => object_tile(() -> RoundPolarizationFilter(inch)),
    "LinearPolarizer" => object_tile(() -> RoundLinearPolarizer(inch, 1mm, 1mm, n_crown)),
    "RoundLinearPolarizer" => object_tile(() -> RoundLinearPolarizer(inch, 1mm, 1mm, n_crown)),
    # dummies
    "MeshDummy" => object_tile(() -> MeshDummy(joinpath(@__DIR__, "..", "detector_assets", "FDS010.stl"))),
    "NonInteractableObject" => object_tile(() -> NonInteractableObject(BMO.CylinderSDF(6mm, 50mm))),
    "IntersectableObject" => object_tile(() -> IntersectableObject(BMO.BoxSDF(20mm, 10mm, 20mm))),
]

#=
Coverage check: every exported component constructor must have a catalog entry
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

"Exported component constructors of the relevant sections of `src/Exports.jl`."
function exported_components()
    names = String[]
    section = ""
    for line in eachline(joinpath(pkgdir(BeamletOptics), "src", "Exports.jl"))
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
                is_component_constructor(name) && push!(names, id.match)
            end
        end
    end
    return unique(names)
end

const CATALOG_JSON = joinpath(@__DIR__, "..", "..", "components", "catalog.json")

let names = first.(TILES)
    check_tiles(CATALOG_JSON, names)
    not_listed = setdiff(exported_components(), catalog_names(CATALOG_JSON))
    isempty(not_listed) || @warn "Exported components missing from catalog.json" not_listed
    for (name, draw!) in TILES
        save_tile(name, draw!)
    end
end

end # module
