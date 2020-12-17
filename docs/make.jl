using IOSystems
using Documenter

makedocs(;
    modules=[IOSystems],
    authors="Hans Würfel <git@wuerfel.io> and contributors",
    repo="https://github.com/hexaeder/IOSystems.jl/blob/{commit}{path}#L{line}",
    sitename="IOSystems.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://hexaeder.github.io/IOSystems.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/hexaeder/IOSystems_prototype",
    devbranch="main"
)
