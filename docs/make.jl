using CairoMakie
using GLMakie
using BeamletOptics
using Documenter
using DocumenterCitations

CairoMakie.activate!()

DocMeta.setdocmeta!(BeamletOptics, :DocTestSetup, :(using BeamletOptics); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(;
    modules=[BeamletOptics],
    authors="Hugo Uittenbosch <hugo.uittenbosch@dlr.de> and contributors",
    sitename="BeamletOptics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://StackEnjoyer.github.io/BeamletOptics.jl",
        edit_link="master",
        assets=String[],
        size_threshold_ignore=["reference.md"],
        sidebar_sitename = false,
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => Any[
            "Beam expander" => "tutorials/expander.md",
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
        ],
        "Components" => Any[
            "Overview" => "components/components.md",
            "Mirrors" => "components/mirrors.md",
            "Lenses" => "components/lenses.md",
            "Beamsplitters" => "components/beamsplitters.md",
            "Detectors" => "components/detectors.md",
        ],
        "Developer Documentation" => Any[
            "Dev. guide" => "guide.md",
            "API design" => "design.md",
        ],
        "Reference" => "reference.md"
    ],
    plugins=[bib],
)

deploy_flag = false

if deploy_flag
    deploydocs(;
    repo="github.com/StackEnjoyer/BeamletOptics.jl.git",
    devbranch="master",
    )
end