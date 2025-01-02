using OpenSCI
using Documenter

DocMeta.setdocmeta!(OpenSCI, :DocTestSetup, :(using OpenSCI); recursive=true)

makedocs(;
    modules=[OpenSCI],
    authors="Nick Mayhall",
    sitename="OpenSCI.jl",
    format=Documenter.HTML(;
        canonical="https://nmayhall-vt.github.io/OpenSCI.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/nmayhall-vt/OpenSCI.jl",
    devbranch="main",
)
