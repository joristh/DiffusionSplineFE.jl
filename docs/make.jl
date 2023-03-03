using DiffusionSplineFE
using Documenter

DocMeta.setdocmeta!(DiffusionSplineFE, :DocTestSetup, :(using DiffusionSplineFE); recursive=true)

makedocs(;
    modules=[DiffusionSplineFE],
    authors="Joris Thiel <joris.thiel@tum.de>",
    repo="https://github.com/joristh/DiffusionSplineFE.jl/blob/{commit}{path}#{line}",
    sitename="DiffusionSplineFE.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://joristh.github.io/DiffusionSplineFE.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/joristh/DiffusionSplineFE.jl",
    devbranch="main",
)
