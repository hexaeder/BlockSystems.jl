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

    @testset "test creation of namespace map" begin
        import IOSystems.create_namespace_map
        @parameters t
        @variables x(t) i(t) o(t)
        @derivatives D'~t
        eqs  = [D(x) ~ i, o~i]

        @variables x2(t) i2(t) o(t)
        eqs2  = [D(x) ~ i2, D(x2) ~ i2, o~i2]

        iob1 = IOBlock(eqs, [i], [o], name=:iob1)
        iob2 = IOBlock(eqs2, [i2], [o], name=:iob2)

        @test Set(iob1.inputs) == Set(i)
        @test Set(iob2.inputs) == Set(i2)
        @test Set(iob1.outputs) == Set(o)
        @test Set(iob2.outputs) == Set(o)
        @test Set(iob1.istates) == Set(x)
        @test Set(iob2.istates) == Set([x, x2])

        in  = create_namespace_map([iob1, iob2], :inputs)
        in_ex = [iob1.i => i, iob2.i2 => i2]
        @test Set(in) == Set(in_ex)

        out = create_namespace_map([iob1, iob2], :outputs)
        out_ex = [iob1.o => iob1.o, iob2.o => iob2.o]
        @test Set(out) == Set(out_ex)

        var = create_namespace_map([iob1, iob2], :istates)
        var_ex = [iob1.x => iob1.x, iob2.x => iob2.x, iob2.x2 => x2]
        @test Set(var) == Set(var_ex)
    end
end
