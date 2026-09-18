using CairoMakie
using GLMakie
using BeamletOptics
using Documenter
using DocumenterCitations
using DocumenterVitepress
import NodeJS_20_jll

include(joinpath(@__DIR__, "DocUtils.jl"))

# DocumenterVitepress runs `npm install` through NodeJS_20_jll, but JLLWrappers does not
# put the artifact's `bin` directory on PATH. npm postinstall scripts that spawn `node`
# themselves (esbuild) then fail, so put it there ourselves.
if Sys.iswindows()
    ENV["PATH"] = string(dirname(NodeJS_20_jll.node_path), ";", ENV["PATH"])
end

# DocumenterCitations 1.5 wraps every in-text citation in a `CitationSiteNode`, an HTML
# anchor the bibliography backlinks point at. DocumenterVitepress only handles the
# `BibliographyNode`, so without this method the node itself ends up in the markdown as
# `DocumenterCitations.CitationSiteNode("...")`. Emit the anchor, then the citation link,
# mirroring what the LaTeX writer of DocumenterCitations does.
function DocumenterVitepress.render(
    io::IO,
    mime::MIME"text/plain",
    node::Documenter.MarkdownAST.Node,
    citation_site::DocumenterCitations.CitationSiteNode,
    page,
    doc;
    kwargs...
)
    print(io, "<a id=\"", citation_site.id, "\"></a>")
    DocumenterVitepress.render(io, mime, node, node.children, page, doc; kwargs...)
    return nothing
end

DocMeta.setdocmeta!(
    BeamletOptics,
    :DocTestSetup,
    :(using BeamletOptics);
    recursive=true
)

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(;
    modules=[BeamletOptics],
    authors="Hugo Uittenbosch <hugo.uittenbosch@dlr.de>, Oliver Kliebisch <oliver.kliebisch@dlr.de> and contributors",
    sitename="BeamletOptics.jl",
    format=DocumenterVitepress.MarkdownVitepress(;
        repo="github.com/JuliaPhysics/BeamletOptics.jl",
        devbranch="master",
        devurl="dev",
    ),
    pagesonly=true,
    pages=[
        "Home" => "index.md",
        "Getting started" => Any[
            "Tutorials" => Any[
                "Beam expander"             => joinpath("tutorials", "expander.md"),
                "Miniature microscope"      => joinpath("tutorials", "microscope.md"),
                "Michelson interferometer"  => joinpath("tutorials", "michelson.md"),
                "Raman spectroscopy"        => joinpath("tutorials", "openraman.md"),
            ],
            "Examples" => Any[
                "Spherical lenses"          => joinpath("examples", "spherical_lenses.md"),
                "Aspherical lenses"         => joinpath("examples", "aspherical_lenses.md"),
                "Double Gauss lens"         => joinpath("examples", "double_gauss.md"),
                "Lens groups"               => joinpath("examples", "lens_groups.md"),
                "Double slit"               => joinpath("examples", "double_slit.md"),
                "Point spread functions"    => joinpath("examples", "psf.md"),
            ],
        ],
        "Basics" => Any[
            "Introduction"                  => joinpath("basics", "intro.md"),
            "Rays"                          => joinpath("basics", "rays.md"),
            "Beams" => Any[
                "Basic beam"                => joinpath("basics", "beams", "beams.md"),
                "Stigmatic Gaussian"        => joinpath("basics", "beams", "stigmatic_beam.md"),
                "Astigmatic Gaussian"       => joinpath("basics", "beams", "astigmatic_beam.md"),
                "Beam groups"               => joinpath("basics", "beams", "beam_groups.md"),
            ],
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

# On Windows DocumenterVitepress only runs `npm install` and tells the user to install
# Node.js system-wide instead of building the site. Do the build here with the JLL's Node.
if Sys.iswindows()
    # The drive letter must be upper case: VS Code starts Julia in `c:\...`, and with a
    # lower case drive letter the VitePress SSR build fails with ERR_MODULE_NOT_FOUND
    # for chunks in `.vitepress/.temp`.
    cd(uppercasefirst(@__DIR__)) do
        run(`$(NodeJS_20_jll.node()) node_modules/vitepress/bin/vitepress.js build build/.documenter`)
    end
end

DocumenterVitepress.deploydocs(;
    repo="github.com/JuliaPhysics/BeamletOptics.jl.git",
    target=joinpath(@__DIR__, "build"),
    devbranch="master",
    push_preview=false,
)
