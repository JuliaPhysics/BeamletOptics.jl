using SCDI
using Documenter

DocMeta.setdocmeta!(SCDI, :DocTestSetup, :(using SCDI); recursive=true)

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
    ],
)
