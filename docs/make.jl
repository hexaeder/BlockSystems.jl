using IOSystems
using Documenter
using Literate

using Plots
using DifferentialEquations
using ModelingToolkit

# generate examples
spacecraft = joinpath(@__DIR__, "..", "examples", "spacecraft.jl")
kuramoto = joinpath(@__DIR__, "..", "examples", "kuramoto_network.jl")
OUTPUT = joinpath(@__DIR__, "src/generated")

for ex in [spacecraft, kuramoto]
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
                       "Kuramoto Network" => "generated/kuramoto_network.md"]
    ],
)

deploydocs(;
    repo="github.com/hexaeder/IOSystems_prototype",
    devbranch="main"
)
