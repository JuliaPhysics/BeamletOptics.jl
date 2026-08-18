using CairoMakie
using GLMakie
using BeamletOptics
using Documenter
using DocumenterCitations

include(joinpath(@__DIR__, "DocUtils.jl"))

DocMeta.setdocmeta!(
    BeamletOptics,
    :DocTestSetup,
    :(using BeamletOptics);
    recursive = true
)

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(;
    modules = [BeamletOptics],
    authors = "Hugo Uittenbosch <hugo.uittenbosch@dlr.de>, Oliver Kliebisch <oliver.kliebisch@dlr.de> and contributors",
    sitename = "BeamletOptics.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://JuliaPhysics.github.io/BeamletOptics.jl",
        edit_link = "master",
        assets = String[],
        size_threshold_ignore = ["reference.md"],
        sidebar_sitename = false
    ),
    pagesonly = true,
    pages = [
        "Home" => "index.md",
        "Getting started" => Any[
            "Tutorials" => Any[
                "Beam expander" => joinpath("tutorials", "expander.md"),
                "Miniature microscope" => joinpath("tutorials", "microscope.md"),
                "Michelson interferometer" => joinpath("tutorials", "michelson.md"),
                "Wedge trap" => joinpath("tutorials", "wolit_trap.md"),
                "Fabry-Pérot cavity" => joinpath("tutorials", "fabry_perot.md")
            ],
            "Examples" => Any[
                "Spherical lenses" => joinpath("examples", "spherical_lenses.md"),
                "Aspherical lenses" => joinpath("examples", "aspherical_lenses.md"),
                "Double Gauss lens" => joinpath("examples", "double_gauss.md"),
                "Lens groups" => joinpath("examples", "lens_groups.md"),
                "Double slit" => joinpath("examples", "double_slit.md"),
                "Optical coatings" => joinpath("examples", "coatings.md"),
                "Polarizing optics" => joinpath("examples", "polarizing_optics.md")
            ]
        ],
        "Basics" => Any[
            "Introduction" => joinpath("basics", "intro.md"),
            "Rays" => joinpath("basics", "rays.md"),
            "Beams" => Any[
                "Basic beam" => joinpath("basics", "beams", "beams.md"),
                "Stigmatic Gaussian" => joinpath("basics", "beams", "stigmatic_beam.md"),
                "Astigmatic Gaussian" => joinpath("basics", "beams", "astigmatic_beam.md"),
                "Beam groups" => joinpath("basics", "beams", "beam_groups.md")
            ],
            "Optical components" => Any[
                "Overview" => joinpath("basics", "components", "components.md"),
                "Mirrors" => joinpath("basics", "components", "mirrors.md"),
                "Lenses" => joinpath("basics", "components", "lenses.md"),
                "Beamsplitters" => joinpath("basics", "components", "beamsplitters.md"),
                "Detectors" => joinpath("basics", "components", "detectors.md"),
                "Polarizing components" => joinpath("basics", "components", "polarizers.md"),
                "Coatings" => joinpath("basics", "components", "coatings.md"),
                "Materials & Absorption" => joinpath("basics", "components", "materials.md")
            ],
            "Optical systems" => joinpath("basics", "systems.md"),
            "Visualization" => joinpath("basics", "render.md")
        ],
        "Developer Documentation" => Any[
            "Developer guide" => Any[
                "Contributing" => joinpath("api", "contribute.md"),
                "Documentation development" => joinpath("api", "docdev.md")
            ],
            "API design" => Any[
                "Introduction" => joinpath("api", "api.md"),
                "Conventions" => joinpath("api", "conventions.md"),
                "Core design" => joinpath("api", "core.md"),
                "Geometry" => Any[
                    "Geometry representation" => joinpath("api", "geometry.md"),
                    "Meshes" => joinpath("api", "meshes.md"),
                    "SDFs" => joinpath("api", "sdfs.md")
                ]
            ]
        ],
        "Reference" => "reference.md"
    ],
    plugins = [bib]
)

deploydocs(;
    repo = "github.com/JuliaPhysics/BeamletOptics.jl.git",
    devbranch = "master",
    push_preview = false
)
