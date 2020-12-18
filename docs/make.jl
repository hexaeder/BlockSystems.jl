using IOSystems
using Documenter
using Literate

# generate examples
EXAMPLE = joinpath(@__DIR__, "..", "examples", "spacecraft.jl")
OUTPUT = joinpath(@__DIR__, "src/generated")

Literate.markdown(EXAMPLE, OUTPUT)
Literate.script(EXAMPLE, OUTPUT)

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
        "Example" => "generated/spacecraft.md"
    ],
)

deploydocs(;
    repo="github.com/hexaeder/IOSystems_prototype",
    devbranch="main"
)
