using CairoMakie
using BeamletOptics
using Documenter
using DocumenterCitations

try
    rm(joinpath(@__DIR__, "build"), recursive=true)
    @info "Deleted build folder..."
catch e
    if isa(e, Base.IOError)
        @info "Can't delete build folder, does not exist..."
    else
        rethrow(e)
    end
end

CairoMakie.activate!()

DocMeta.setdocmeta!(BeamletOptics, :DocTestSetup, :(using BeamletOptics); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(;
    modules=[BeamletOptics],
    authors="Hugo Uittenbosch <hugo.uittenbosch@dlr.de> and contributors",
    repo="https://gitlab.dlr.de/optical-air-data/dispersionsinterferometer/scdi-sim/-/blob/{commit}{path}#L{line}",
    sitename="BeamletOptics",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
        size_threshold_ignore=["reference.md"],
        sidebar_sitename = false,
    ),
    pages=[
        "Home" => "index.md",
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
        "Tutorials" => Any[
            "Beam expander" => "tutorials/expander.md",
            "Michelson interferometer" => "tutorials/michelson.md"
        ],
        "API design" => "design.md",
        "Examples" => Any[
            "Spherical lenses" => "examples/spherical_lenses.md",
            "Aspherical lenses" => "examples/aspherical_lenses.md",
            "Double Gauss lens" => "examples/double_gauss.md",
            "Lens groups" => "examples/lens_groups.md",
        ],
        "Dev. guide" => "guide.md",
        "Reference" => "reference.md"
    ],
    plugins=[bib],
)
