#! julia

using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=dirname(@__DIR__))) # adds the package this script is called from
Pkg.instantiate()
Pkg.update()

run = true
while run
    include("make.jl")

    using LiveServer; serve(dir=joinpath(@__DIR__, "build"))

    println("Run again? Enter! Exit witn 'q'.")
    if readline() == "q"
        global run = false
    end
end
