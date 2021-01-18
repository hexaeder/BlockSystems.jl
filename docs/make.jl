using IOSystems
using Documenter
using Literate

using Plots
using DifferentialEquations
using ModelingToolkit

# generate examples
examples = [
    joinpath(@__DIR__, "..", "examples", "spacecraft.jl"),
    joinpath(@__DIR__, "..", "examples", "kuramoto_network.jl"),
    joinpath(@__DIR__, "..", "examples", "kuramoto_without_nd.jl"),
]
OUTPUT = joinpath(@__DIR__, "src/generated")

for ex in examples
    Literate.markdown(ex, OUTPUT)
end

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
        "Examples" => ["Spacecraft" => "generated/spacecraft.md",
                       "Kuramoto Network" => "generated/kuramoto_network.md",
                       "Kuramoto without ND.jl" => "generated/kuramoto_without_nd.md"]
    ],
)

deploydocs(;
    repo="github.com/hexaeder/IOSystems_prototype",
    devbranch="main"
)
