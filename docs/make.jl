using CairoMakie
using GLMakie
using BeamletOptics
using Documenter
using DocumenterCitations

DocMeta.setdocmeta!(BeamletOptics, :DocTestSetup, :(using BeamletOptics); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(;
    modules=[BeamletOptics],
    authors="Hugo Uittenbosch <hugo.uittenbosch@dlr.de>, Oliver Kliebisch <oliver.kliebisch@dlr.de> and contributors",
    sitename="BeamletOptics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaPhysics.github.io/BeamletOptics.jl",
        edit_link="master",
        assets=String[],
        size_threshold_ignore=["reference.md"],
        sidebar_sitename = false,
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => Any[
            "Beam expander" => "tutorials/expander.md",
            "Miniature microscope" => "tutorials/microscope.md",
            "Michelson interferometer" => "tutorials/michelson.md"
        ],
        "Examples" => Any[
            "Spherical lenses" => "examples/spherical_lenses.md",
            "Aspherical lenses" => "examples/aspherical_lenses.md",
            "Double Gauss lens" => "examples/double_gauss.md",
            "Lens groups" => "examples/lens_groups.md",
        ],
        "Basics" => Any[
            "Introduction" => "basics/intro.md",
            "Rays" => "basics/rays.md",
            "Beams" => "basics/beams.md",
            "Optical elements" => "basics/elements.md",
            "Optical systems" => "basics/systems.md",
            "Visualization" => "basics/render.md"
        ],
        "Components" => Any[
            "Overview" => "components/components.md",
            "Mirrors" => "components/mirrors.md",
            "Lenses" => "components/lenses.md",
            "Beamsplitters" => "components/beamsplitters.md",
            "Detectors" => "components/detectors.md",
            "Polarizing components" => "components/polarizers.md",
        ],
        "Developer Documentation" => Any[
            "Developer guide" => Any[
                "Contributing" => "api/contribute.md",
                "Documentation development" => "api/docdev.md"
            ],
            "API design" => Any[
                "Introduction" => "api/api.md",
                "Conventions" => "api/conventions.md",
                "Core design" => "api/core.md",
                "Geometry" => "api/geometry.md",
            ],
        ],
        "Reference" => "reference.md"
    ],
    plugins=[bib],
)

deploydocs(;
    repo="github.com/JuliaPhysics/BeamletOptics.jl.git",
    devbranch="master",
    push_preview=true, # set to false before merging PR into master
)
