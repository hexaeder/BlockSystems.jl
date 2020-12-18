@info "Testes of IOSystems.jl"

@testset "IOSystems.jl" begin
    @testset "collect and resolve namespace" begin
        @parameters t
        @variables x(t) i(t) o(t)
        @derivatives D'~t
        eqs  = [D(x) ~ i, o~i]

        @variables x2(t) i2(t) o(t)
        eqs2  = [D(x) ~ i2, o~i2]

        iob1 = IOBlock(eqs, [i], [o], name=:iob1)
        iob2 = IOBlock(eqs2, [i2], [o], name=:iob2)

        @test Set(iob1.inputs) == Set(i)
        @test Set(iob2.inputs) == Set(i2)
        @test Set(iob1.outputs) == Set(o)
        @test Set(iob2.outputs) == Set(o)

        in = collect_and_resolve_namespace([iob1, iob2], :inputs)
        out = collect_and_resolve_namespace([iob1, iob2], :outputs)
        var = collect_and_resolve_namespace([iob1, iob2], :interns)
        @test Set(keys(in)) == Set([i, i2])
        @test Set(keys(out)) == Set([iob1.o, iob2.o])
        @test Set(keys(var)) == Set([iob1.x, iob2.x])

        ios = IOSystem([iob1.o => iob1.i], [iob1, iob2], name=:ios)
        ios.name
        ios.inputs
        ios.interns
        ios.outputs

        ios.iob1₊x
    end
end

using ModelingToolkit
@parameters t
@variables x(t) i(t) o(t)
@derivatives D'~t

iob1 = IOBlock([D(x) ~ i, o ~ x], # define the equation
               [i], # inputs of the system
               [o], # outputs of the system
               name = :block1)
@show iob1.inputs
@show iob1.outputs
@show iob1.istates
print(iob1.inputs)

iob1 = IOBlock([D(x) ~ i, o ~ x], # define the equation
               [i], # inputs of the system
               [o], # outputs of the system
               name = :block1)


using IOSystems
using ModelingToolkit
@parameters t, m
@variables x(t) F(t) v(t)
@derivatives D'~t

spacecraft = IOBlock([D(v) ~ F/m, o ~ x], # define the equation
                     [F], # inputs of the system
                     [x], # outputs of the system
                     name = :spacecraft)
@show spacecraft.inputs
@show spacecraft.outputs
@show spacecraft.istates
