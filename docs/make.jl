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
    pagesonly=true,
    pages=[
        "Home" => "index.md",
        "Getting started" => Any[
            "Tutorials" => Any[
                "Beam expander"             => joinpath("tutorials", "expander.md"),
                "Miniature microscope"      => joinpath("tutorials", "microscope.md"),
                "Michelson interferometer"  => joinpath("tutorials", "michelson.md"),
            ],
            "Examples" => Any[
                "Spherical lenses"          => joinpath("examples", "spherical_lenses.md"),
                "Aspherical lenses"         => joinpath("examples", "aspherical_lenses.md"),
                "Double Gauss lens"         => joinpath("examples", "double_gauss.md"),
                "Lens groups"               => joinpath("examples", "lens_groups.md"),
            ],
        ],
        "Basics" => Any[
            "Introduction"                  => joinpath("basics", "intro.md"),
            "Rays"                          => joinpath("basics", "rays.md"),
            "Beams"                         => joinpath("basics", "beams.md"),
            "Optical components" => Any[
                "Overview"                  => joinpath("basics", "components", "components.md"),
                "Mirrors"                   => joinpath("basics", "components", "mirrors.md"),
                "Lenses"                    => joinpath("basics", "components", "lenses.md"),
                "Beamsplitters"             => joinpath("basics", "components", "beamsplitters.md"),
                "Detectors"                 => joinpath("basics", "components", "detectors.md"),
                "Polarizing components"     => joinpath("basics", "components", "polarizers.md"),
            ],
            "Optical systems"               => joinpath("basics", "systems.md"),
            "Visualization"                 => joinpath("basics", "render.md"),
        ],
        "Developer Documentation" => Any[
            "Developer guide" => Any[
                "Contributing"              => joinpath("api", "contribute.md"),
                "Documentation development" => joinpath("api", "docdev.md"),
            ],
            "API design" => Any[
                "Introduction"              => joinpath("api", "api.md"),
                "Conventions"               => joinpath("api", "conventions.md"),
                "Core design"               => joinpath("api", "core.md"),
                "Geometry" => Any[
                    "Geometry representation"   => joinpath("api", "geometry.md"),
                    "Meshes"                    => joinpath("api", "meshes.md"),
                    "SDFs"                      => joinpath("api", "sdfs.md"),
                ],
            ],
        ],
        "Reference" => "reference.md"
    ],
    plugins=[bib],
)

deploydocs(;
    repo="github.com/JuliaPhysics/BeamletOptics.jl.git",
    devbranch="master",
    push_preview=false,
)
