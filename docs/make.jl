using CairoMakie
using GLMakie
using BeamletOptics
using Documenter
using DocumenterCitations
using DocumenterVitepress
using Literate
import NodeJS_20_jll

include(joinpath(@__DIR__, "DocUtils.jl"))

# Literate-based tutorials: the source of truth is the `.jl` file in `docs/literate`.
# `Literate.markdown` regenerates the corresponding page in `docs/src/tutorials` (see
# `.gitignore`), and `Literate.script`/`Literate.notebook` produce the downloadable
# `.jl`/`.ipynb` files served via the `vitepress_assets` hook in `DocUtils.jl`.
let
    lit_src = joinpath(@__DIR__, "literate", "laser_alignment.jl")
    out_md  = joinpath(@__DIR__, "src", "tutorials")
    out_dl  = joinpath(@__DIR__, "src", "assets", "downloads")
    Literate.markdown(lit_src, out_md; documenter=true, credit=false)
    Literate.script(lit_src, out_dl; credit=false)
    Literate.notebook(lit_src, out_dl; execute=false, credit=false)
end

# DocumenterVitepress runs `npm install` through NodeJS_20_jll, but JLLWrappers does not
# put the artifact's `bin` directory on PATH. npm postinstall scripts that spawn `node`
# themselves (esbuild) then fail, so put it there ourselves.
if Sys.iswindows()
    ENV["PATH"] = string(dirname(NodeJS_20_jll.node_path), ";", ENV["PATH"])
end

# DocumenterVitepress only merges plugin `vitepress_dependencies` (used below for mermaid)
# into `docs/package.json` when that file already exists; otherwise it first `cp`s its own
# read-only template into place. On Windows, `cp` preserves the source's read-only
# attribute, so the subsequent write to merge in the mermaid deps fails with
# `IOError: ... Permission denied`. Pre-seed a writable copy ourselves so DV finds an
# existing (non-read-only) `package.json` and skips its own `cp`. `docs/package.json` is
# gitignored and regenerated on every build.
if Sys.iswindows()
    pkg_json = joinpath(@__DIR__, "package.json")
    if !isfile(pkg_json)
        template = joinpath(dirname(pathof(DocumenterVitepress)), "..", "template", "package.json")
        write(pkg_json, read(template))
    end
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

# On Windows, `@contents` listings of pages in subfolders are broken twice. Documenter
# matches `Pages` against `relpath`s with backslashes, so pages must be given as
# `joinpath("beams", "beams.md")` (see `basics/intro.md`) rather than "beams/beams.md".
# DocumenterVitepress then writes that relpath into the link, and VitePress rejects
# `beams\beams` as a dead link. This is DocumenterVitepress' method with the path
# separators normalised to `/`.
if Sys.iswindows()
    function DocumenterVitepress.render(
        io::IO,
        ::MIME"text/plain",
        node::Documenter.MarkdownAST.Node,
        contents::Documenter.ContentsNode,
        page,
        doc;
        kwargs...
    )
        current_path = nothing
        for (count, path, anchor) in contents.elements
            path = replace(DocumenterVitepress.mdext(path), '\\' => '/')
            header = anchor.object
            anchor_frag = DocumenterVitepress.vitepress_anchor(Documenter.anchor_fragment(anchor))
            url = replace(string(path, anchor_frag), " " => "%20")
            link = DocumenterVitepress.Markdown.Link(replace(anchor.id, "-" => " "), url)
            level = header.level
            if path != current_path
                level = 1
                current_path = path
            end
            print(io, "    "^(level - 1), "- ")
            println(io, replace(DocumenterVitepress.Markdown.plaininline(link), ".md#" => "#"))
        end
        return println(io)
    end
end

DocMeta.setdocmeta!(
    BeamletOptics,
    :DocTestSetup,
    :(using BeamletOptics);
    recursive=true
)

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(;
    modules=[BeamletOptics, Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)],
    authors="Hugo Uittenbosch <hugo.uittenbosch@dlr.de>, Oliver Kliebisch <oliver.kliebisch@dlr.de> and contributors",
    sitename="BeamletOptics.jl",
    format=DocumenterVitepress.MarkdownVitepress(;
        repo="github.com/JuliaPhysics/BeamletOptics.jl",
        devbranch="master",
        devurl="dev",
    ),
    pagesonly=true,
    warnonly=[:missing_docs],
    pages=[
        "Home" => "index.md",
        "Getting started" => Any[
            "Overview" => joinpath("tutorials", "index.md"),
            "Tutorials" => Any[
                "Laser alignment"           => joinpath("tutorials", "laser_alignment.md"),
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
            "Visualization" => Any[
                "Overview"                  => joinpath("basics", "visualization", "overview.md"),
                "Rays and beams"            => joinpath("basics", "visualization", "beams.md"),
                "Gaussian beamlets"         => joinpath("basics", "visualization", "gaussian.md"),
                "Components and systems"    => joinpath("basics", "visualization", "components.md"),
                "Scene and camera"          => joinpath("basics", "visualization", "camera.md"),
            ],
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
            "Roadmap" => "roadmap.md",
        ],
        "Reference" => "reference.md"
    ],
    plugins=[bib, DocUtils.BMODocsExtras()],
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
