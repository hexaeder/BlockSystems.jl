using Test
using IOSystems
using IOSystems: isunique
using ModelingToolkit
using ModelingToolkit: vars, value
using LightGraphs

@info "Testes of IOSystems.jl"

@testset "IOSystems.jl" begin
    @testset "namespaced accessors" begin
        @parameters t i1(t) i2(t) a
        @variables x1(t) x2(t) o1(t) o2(t)
        @derivatives D'~t
        eqs  = [D(x1) ~ a*i1, o1~i1, D(x2) ~ i2, o2~i2]
        iob = IOBlock(eqs, [i1, i2], [o1, o2], name=:ns)
        @test Set(IOSystems.namespace_inputs(iob)) == Set([iob.i1, iob.i2])
        @test Set(IOSystems.namespace_outputs(iob)) == Set([iob.o1, iob.o2])
        @test Set(IOSystems.namespace_istates(iob)) == Set([iob.x1, iob.x2])
        @test Set(IOSystems.namespace_iparams(iob)) == Set([iob.a])
    end

    @testset "creation of IOBlocks" begin
        @parameters t i1(t) i2(t) a b
        @variables x1(t) x2(t) o1(t) o2(t)
        @derivatives D'~t
        eqs = [D(x1) ~ a*i1,
               D(x2) ~ b*i2,
               o1 ~ a*x1,
               o2 ~ b*x2]

        iob = IOBlock(eqs, [i1, i2], [o1, o2], name=:iob)
        @test Set(iob.inputs) == Set([i1, i2])
        @test Set(iob.iparams) == Set([a, b])
        @test Set(iob.istates) == Set([x1, x2])
        @test Set(iob.outputs) == Set([o1, o2])

        @test_throws AssertionError IOBlock(eqs, [x1], [o1,o2])
        @test_throws AssertionError IOBlock(eqs, [i1,i2], [i1,o1,o2])
    end

    @testset "test of create_namespace_promotions" begin
        using ModelingToolkit: value
        using IOSystems: create_namespace_promotions
        @parameters t
        Aa, Ba, Ab, Bc = value.(@parameters A₊a B₊a(t) A₊b B₊c)
        Ax, Bx, Ay, Bz = value.(@variables A₊x(t) B₊x(t) A₊y(t) B₊z(t))
        b, c = value.(@parameters b c)
        y, z = value.(@variables y(t) z(t))

        prom = create_namespace_promotions([Aa, Ba, Ab, Bc, Ax, Bx, Ay, Bz], [])
        @test Set(values(prom)) == Set([Aa, Ba, b, c, Ax, Bx, y, z])

        prom = create_namespace_promotions([Aa, Ba, Ab, Bc, Ax, Bx, Ay, Bz], [b, y])
        @test Set(values(prom)) == Set([Aa, Ba, Ab, c, Ax, Bx, Ay, z])
    end

    @testset "IOBlock from other IOBlock" begin
        @parameters t i1(t) i2(t) a b
        @variables x1(t) x2(t) o1(t) o2(t)
        @derivatives D'~t
        eqs = [D(x1) ~ a*i1,
               D(x2) ~ b*i2,
               o1 ~ a*x1,
               o2 ~ b*x2]
        iob1 = IOBlock(eqs, [i1, i2], [o1, o2], name=:iob1)
        iob2 = IOBlock(iob1, name=:iob2)
        iob3 = IOBlock(iob1)
        @test iob1.name != iob2.name != iob3.name
        @test Set(iob1.inputs) == Set(iob2.inputs) == Set(iob3.inputs)
        @test Set(iob1.iparams) == Set(iob2.iparams) == Set(iob3.iparams)
        @test Set(iob1.istates) == Set(iob2.istates) == Set(iob3.istates)
        @test Set(iob1.outputs) == Set(iob2.outputs) == Set(iob3.outputs)
        @test iob1.system.eqs == iob2.system.eqs == iob3.system.eqs
        @test iob1.system.name == iob1.name
        @test iob2.system.name == iob2.name
        @test iob3.system.name == iob3.name
    end

    @testset "test creation of namespace map" begin
        @parameters t i(t) a b
        @variables x(t) o(t)
        @derivatives D'~t
        eqs  = [D(x) ~ a*i, o ~ b*i]

        @parameters i2(t) a
        @variables x2(t) o(t)
        eqs2  = [D(x) ~ a*i2, D(x2) ~ i2, o~i2]

        iob1 = IOBlock(eqs, [i], [o], name=:iob1)
        iob2 = IOBlock(eqs2, [i2], [o], name=:iob2)

        @test Set(iob1.inputs) == Set(i)
        @test Set(iob2.inputs) == Set(i2)
        @test Set(iob1.outputs) == Set(o)
        @test Set(iob2.outputs) == Set(o)
        @test Set(iob1.istates) == Set(x)
        @test Set(iob2.istates) == Set([x, x2])
        @test Set(iob1.iparams) == Set([a, b])
        @test Set(iob2.iparams) == Set(a)

        sys = IOSystem([], [iob1, iob2])
        @test Set(sys.inputs) == Set([i, i2])
        @test Set(sys.istates) == Set([iob1.x, iob2.x, x2])
        @test Set(sys.iparams) == Set([iob1.a, iob2.a, b])
        @test Set(sys.outputs) == Set([iob1.o, iob2.o])
    end

    @testset "test isunique" begin
        using IOSystems: uniquenames
        @test isunique([1,2,3])
        @test !isunique([1,1,3])

        (t, a, b, mp) = value.(@parameters t a b(t) m)
        (x, y, mv) = value.(@variables x y(t) m)
        @test uniquenames([a, b, x, y, mp])
        @test !uniquenames([a, b, x, y, mp, mv])
    end

    @testset "iosystem asserts" begin
        @parameters t i(t)
        @variables o(t)
        iob1 = IOBlock([o~i],[i],[o],name=:name)
        iob2 = IOBlock([o~i],[i],[o],name=:name)

        # namespace collision
        @test_throws AssertionError IOSystem([iob2.i=>iob1.o], [iob1, iob2])
        iob2 = IOBlock([o~i],[i],[o])
        # mulitiple conneections to same input
        @test_throws AssertionError IOSystem([iob1.i=>iob1.o, iob1.i=>iob2.o], [iob1, iob2])
        # make sure that alle of form input => output
        @test_throws AssertionError IOSystem([iob1.o=>iob1.o], [iob1, iob2])
        @test_throws AssertionError IOSystem([iob1.i=>iob1.i], [iob1, iob2])

        # assert that input maps refere to open inputs
        @test_throws AssertionError IOSystem([iob2.i=>iob1.o], [iob1, iob2], inputs_map = [iob2.i => i])
        # assert that rhs of input map is unique
        iob3 = IOBlock([o~i],[i],[o])
        @test_throws AssertionError IOSystem([iob2.i=>iob1.o],
                                             [iob1, iob2, iob3],
                                             inputs_map = [iob1.i => i, iob3.i => i])

        # test assertions for iparams and istats map
        @parameters t a i(t) b c
        @variables x(t) o(t) y(t)
        @derivatives D'~t
        iob1 = IOBlock([D(x)~ i, o~a*x], [i], [o], name=:iob1)
        iob2 = IOBlock([D(x)~ i, o~a*x], [i], [o], name=:iob2)
        IOSystem([iob2.i=>iob1.o], [iob1, iob2])
        # rhs unique
        @test_throws AssertionError IOSystem([iob2.i=>iob1.o], [iob1, iob2],
                 iparams_map = [iob1.a=>b, iob2.a=>b])
        @test_throws AssertionError IOSystem([iob2.i=>iob1.o], [iob1, iob2],
                 istates_map = [iob1.x=>y, iob2.x=>y])
        # keys in right set
        @test_throws AssertionError IOSystem([iob2.i=>iob1.o], [iob1, iob2],
                 iparams_map = [iob1.x=>b])
        @test_throws AssertionError IOSystem([iob2.i=>iob1.o], [iob1, iob2],
                 istates_map = [iob1.a=>y])

        # tests for outputs
        # rhs unique
        @test_throws AssertionError IOSystem([iob2.i=>iob1.o], [iob1, iob2],
                 outputs_map = [iob1.o=>y, iob2.o=>y])
        @test_throws AssertionError IOSystem([iob2.i=>iob1.o], [iob1, iob2],
                 outputs_map = [iob1.a=>y])
    end

    function test_complete_namespace_promotions(ios)
        eqs = vcat([ModelingToolkit.namespace_equations(iob.system) for iob in ios.systems]...)
        allvars = [(vars(eq.lhs) ∪ vars(eq.rhs)) for eq in eqs]
        allvars = union(allvars...) |> unique
        allvars = setdiff(allvars, [ ios.systems[1].system.iv ])
        allkeys = vcat(collect.(keys.([ios.inputs_map, ios.iparams_map, ios.istates_map, ios.outputs_map]))...)
        allkeys = Set(allkeys) ∪ Set(keys(ios.connections))
        @test isunique(allkeys)
        @test Set(allkeys) == Set(allvars)
    end

    @testset "test creation of systems" begin
        #=
                    *------------*
        in1 --> i1 -|    iob1    |
        in2 --> i2 -|(x1, x2)(a) |-o--*     *-----*
                    *------------*    *-ina-| add |- add ---> out
                    *------------*    *-inb-|     |
        in3 --> i1 -|    iob2    |-o--*     *-----*
        in4 --> i2 -|(x1, x2)(b) |
                    *------------*
        =#
        @parameters t i1(t) i2(t) a b ina(t) inb(t)
        @variables x1(t) x2(t) o(t) add(t)
        @derivatives D'~t
        eqs1  = [D(x1)~a*i1, D(x2)~i2, o~x1+x2]
        iob1 = IOBlock(eqs1, [i1, i2], [o], name=:iob1)

        eqs2  = [D(x1)~b*i1, D(x2)~i2, o~x1+x2]
        iob2 = IOBlock(eqs2, [i1, i2], [o], name=:iob2)

        ioadd = IOBlock([add ~ ina + inb], [ina, inb], [add], name=:add)

        # try with auto namespacing
        sys = IOSystem([ioadd.ina => iob1.o, ioadd.inb => iob2.o],
                       [iob1, iob2, ioadd],
                       name=:sys)
        @test Set(sys.inputs) == Set([iob1.i1, iob1.i2, iob2.i1, iob2.i2])
        @test Set(sys.iparams) == Set([a, b])
        @test Set(sys.istates) == Set([iob1.x1, iob1.x2, iob2.x1, iob2.x2])
        @test Set(sys.outputs) == Set([iob1.o, iob2.o, add])
        test_complete_namespace_promotions(sys)

        # provide maps
        @parameters in1(t) in2(t) in3(t) in4(t) p1 p2
        @variables out(t) y1(t) y2(t)
        sys = IOSystem([ioadd.ina => iob1.o, ioadd.inb => iob2.o],
                       [iob1, iob2, ioadd],
                       inputs_map = Dict(iob1.i1 => in1,
                                         iob1.i2 => in2,
                                         iob2.i1 => in3,
                                         iob2.i2 => in4),
                       iparams_map = Dict(iob1.a => p1,
                                          iob2.b => p2),
                       istates_map = Dict(iob1.x1 => y1,
                                          iob1.x2 => y2,
                                          iob2.x1 => x1,
                                          iob2.x2 => x2),
                       outputs_map = Dict(ioadd.add => out),
                       name=:sys)
        @test Set(sys.inputs) == Set([in1, in2, in3, in4])
        @test Set(sys.iparams) == Set([p1, p2])
        @test Set(sys.istates) == Set([y1, y2, x1, x2, iob1.o, iob2.o])
        @test Set(sys.outputs) == Set([out])
        test_complete_namespace_promotions(sys)

        # provide partial maps
        sys = IOSystem([ioadd.ina => iob1.o, ioadd.inb => iob2.o],
                       [iob1, iob2, ioadd],
                       inputs_map = Dict(iob1.i1 => in1,
                                         iob1.i2 => in2),
                       iparams_map = Dict(iob1.a => p1),
                       istates_map = Dict(iob1.x1 => y1,
                                          iob1.x2 => y2),
                       outputs_map = Dict(ioadd.add => out),
                       name=:sys)
        @test Set(sys.inputs) == Set([in1, in2, i1, i2])
        @test Set(sys.iparams) == Set([p1, b])
        @test Set(sys.istates) == Set([y1, y2, x1, x2, iob1.o, iob2.o])
        @test Set(sys.outputs) == Set([out])
        test_complete_namespace_promotions(sys)
    end

    @testset "is_explicit_algebraic" begin
        using IOSystems: is_explicit_algebraic
        @parameters t
        @variables x(t) y
        @derivatives D'~t
        @test !is_explicit_algebraic(D(x) ~ 0)
        @test !is_explicit_algebraic(D(y) ~ 0)
        @test is_explicit_algebraic(x ~ 0)
        @test is_explicit_algebraic(y ~ 0)
        @test !is_explicit_algebraic(x ~ x^2)
        @test !is_explicit_algebraic(y ~ y^2)
        @test is_explicit_algebraic(y ~ x^2)
        @test is_explicit_algebraic(y ~ x^2)
        @test is_explicit_algebraic(x ~ y^2)
    end

    @testset "substitute" begin
        using IOSystems: eqsubstitute
        @parameters t p i(t) pn inew(t)
        @variables x(t) y xn(t) yn(t)
        @derivatives D'~t
        eq = D(x) ~ x + y + p + i
        @test isequal(eqsubstitute(eq, x=>xn), D(xn) ~ xn + y + p + i)
        @test isequal(eqsubstitute(eq, y=>yn), D(x) ~ x + yn + p + i)
        @test isequal(eqsubstitute(eq, p=>pn), D(x) ~ x + y + pn + i)
        @test isequal(eqsubstitute(eq, i=>inew), D(x) ~ x + y + p + inew)
    end

    @testset "remove namespace of symbols" begin
        using IOSystems: remove_namespace
        using ModelingToolkit: renamespace, to_symbolic, rename
        @parameters t a b(t)
        a = to_symbolic(a)
        b = to_symbolic(b)
        an = rename(a, renamespace(:ns, a.name))
        bn = rename(b, renamespace(:ns, b.op.name))
        @test remove_namespace(:ns, :ns₊n) == :n
        @test remove_namespace("ns", "ns₊n") == "n"
        @test isequal(remove_namespace(:ns, an), a)
        @test isequal(remove_namespace(:ns, a), a)
        @test isequal(remove_namespace(:ns, bn), b)
        @test isequal(remove_namespace(:ns, b), b)

        @test remove_namespace(:ns₊n) == :n
        @test remove_namespace("ns₊n") == "n"
        @test remove_namespace(:ns₊ns2₊n) == :ns2₊n
        @test remove_namespace("ns₊ns2₊n") == "ns2₊n"
        @test isequal(remove_namespace(an), a)
        @test isequal(remove_namespace(a), a)
        @test isequal(remove_namespace(bn), b)
        @test isequal(remove_namespace(b), b)
    end
end
