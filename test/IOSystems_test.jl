@info "Testes of IOSystems.jl"

@testset "IOSystems.jl" begin
    @testset "namespaced accessors" begin
        @parameters t
        @variables x1(t) x2(t) i1(t) i2(t) o1(t) o2(t)
        @derivatives D'~t
        eqs  = [D(x1) ~ i1, o1~i1, D(x2) ~ i2, o2~i2]
        iob = IOBlock(eqs, [i1, i2], [o1, o2], name=:ns)
        @test Set(IOSystems.namespace_inputs(iob)) == Set([iob.i1, iob.i2])
        @test Set(IOSystems.namespace_outputs(iob)) == Set([iob.o1, iob.o2])
        @test Set(IOSystems.namespace_istates(iob)) == Set([iob.x1, iob.x2])
    end

    # @testset "collect and resolve namespace" begin
    #     import IOSystems.collect_and_resolve_namespace
    #     @parameters t
    #     @variables x(t) i(t) o(t)
    #     @derivatives D'~t
    #     eqs  = [D(x) ~ i, o~i]

    #     @variables x2(t) i2(t) o(t)
    #     eqs2  = [D(x) ~ i2, o~i2]

    #     iob1 = IOBlock(eqs, [i], [o], name=:iob1)
    #     iob2 = IOBlock(eqs2, [i2], [o], name=:iob2)

    #     IOSystems.namespace_inputs(iob1)

    #     @test Set(iob1.inputs) == Set(i)
    #     @test Set(iob2.inputs) == Set(i2)
    #     @test Set(iob1.outputs) == Set(o)
    #     @test Set(iob2.outputs) == Set(o)

    #     in = collect_and_resolve_namespace([iob1, iob2], :inputs)
    #     out = collect_and_resolve_namespace([iob1, iob2], :outputs)
    #     var = collect_and_resolve_namespace([iob1, iob2], :interns)
    #     @test Set(keys(in)) == Set([i, i2])
    #     @test Set(keys(out)) == Set([iob1.o, iob2.o])
    #     @test Set(keys(var)) == Set([iob1.x, iob2.x])

    #     ios = IOSystem([iob1.o => iob1.i], [iob1, iob2], name=:ios)
    #     ios.name
    #     ios.inputs
    #     ios.interns
    #     ios.outputs

    #     ios.iob1₊x
    # end
end

