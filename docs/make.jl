using CairoMakie
using SCDI
using Documenter
using DocumenterCitations

CairoMakie.activate!()

try
    mkdir(joinpath(@__DIR__, "src", "assets"))
    @info "Created docs assets folder"
catch
    @info "Assets folder in docs already exists"
end

DocMeta.setdocmeta!(SCDI, :DocTestSetup, :(using SCDI); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(;
    modules=[SCDI],
    authors="Hugo Uittenbosch <hugo.uittenbosch@dlr.de> and contributors",
    repo="https://gitlab.dlr.de/optical-air-data/dispersionsinterferometer/scdi-sim/-/blob/{commit}{path}#L{line}",
    sitename="SCDI Sim",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
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
        "Tutorials" => Any[
            "Beam expander" => "tutorials/expander.md",
        ],
        "API design" => "design.md",
        "Examples" => Any[
            "Spherical lenses" => "examples/spherical_lenses.md",
            "Aspherical lenses" => "examples/aspherical_lenses.md",
            "Double Gauss Lens" => "examples/double_gauss.md",
            "Lens groups" => "examples/lens_groups.md",
        ],
        "Reference" => "reference.md"
    ],
    plugins=[bib],
)
