using Test
using IOSystems
using ModelingToolkit
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

    @testset "test creation of namespace map" begin
        using IOSystems: create_namespace_map
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

        in  = create_namespace_map([iob1, iob2], :inputs)
        in_ex = [iob1.i => i, iob2.i2 => i2]
        @test Set(in) == Set(in_ex)

        params  = create_namespace_map([iob1, iob2], :iparams)
        params_ex = [iob1.a => iob1.a, iob1.b => b, iob2.a => iob2.a]
        @test Set(params) == Set(params_ex)

        out = create_namespace_map([iob1, iob2], :outputs)
        out_ex = [iob1.o => iob1.o, iob2.o => iob2.o]
        @test Set(out) == Set(out_ex)

        var = create_namespace_map([iob1, iob2], :istates)
        var_ex = [iob1.x => iob1.x, iob2.x => iob2.x, iob2.x2 => x2]
        @test Set(var) == Set(var_ex)
    end

    @testset "test isunique" begin
        using IOSystems: isunique
        @test isunique([1,2,3])
        @test !isunique([1,1,3])
    end

    @testset "iosystem asserts" begin
        @parameters t i
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

        # error of namespace clashes
        @test_throws AssertionError IOSystem([],[iob1, iob2],
                                             inputs_map=[iob1.i => i])
        @test_throws AssertionError IOSystem([],[iob1, iob2],
                                             iparams_map=[iob1.a => a])
        @test_throws AssertionError IOSystem([],[iob1, iob2],
                                             istates_map=[iob1.x => x])
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
        @show iob1.inputs iob1.istates iob1.iparams iob1.outputs;

        eqs2  = [D(x1)~b*i1, D(x2)~i2, o~x1+x2]
        iob2 = IOBlock(eqs2, [i1, i2], [o], name=:iob2)
        @show iob2.inputs iob2.istates iob2.iparams iob2.outputs;

        ioadd = IOBlock([add ~ ina + inb], [ina, inb], [add], name=:add)
        @show ioadd.inputs ioadd.istates ioadd.iparams ioadd.outputs;

        # try with auto namespacing
        sys = IOSystem([ioadd.ina => iob1.o, ioadd.inb => iob2.o],
                       [iob1, iob2, ioadd],
                       name=:sys)
        @test Set(sys.inputs) == Set([iob1.i1, iob1.i2, iob2.i1, iob2.i2])
        @test Set(sys.iparams) == Set([a, b])
        @test Set(sys.istates) == Set([iob1.x1, iob1.x2, iob2.x1, iob2.x2])
        @test Set(sys.outputs) == Set([iob1.o, iob2.o, add])

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
        @test Set(sys.istates) == Set([y1, y2, x1, x2])
        @test Set(sys.outputs) == Set([out])

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
        @test Set(sys.istates) == Set([y1, y2, x1, x2])
        @test Set(sys.outputs) == Set([out])
    end

    @testset "test of connect_system" begin
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
        # same system as in testset before
        @parameters t i1(t) i2(t) a b ina(t) inb(t)
        @variables x1(t) x2(t) o(t) add(t)
        @derivatives D'~t
        eqs1  = [D(x1)~a*i1, D(x2)~i2, o~x1+x2]
        iob1 = IOBlock(eqs1, [i1, i2], [o], name=:A)
        eqs2  = [D(x1)~b*i1, D(x2)~i2, o~x1+x2]
        iob2 = IOBlock(eqs2, [i1, i2], [o], name=:B)
        ioadd = IOBlock([add ~ ina + inb], [ina, inb], [add], name=:add)
        @parameters in1(t) in2(t) in3(t) in4(t)
        @variables out
        sys = IOSystem([ioadd.ina => iob1.o, ioadd.inb => iob2.o],
                       [iob1, iob2, ioadd],
                       inputs_map = Dict(iob1.i1 => in1,
                                         iob1.i2 => in2,
                                         iob2.i1 => in3,
                                         iob2.i2 => in4),
                       outputs_map = Dict(ioadd.add => out),
                       name=:sys)
        @show sys.inputs sys.iparams sys.istates sys.outputs

        eqs = connect_system(sys)
        ode = ODESystem(eqs)
        states(ode)
        parameters(ode)

        equation_dependencies(ode)
        variable_dependencies(ode)
        asgraph(ode)

        graph = eqeq_dependencies(asgraph(ode), variable_dependencies(ode))
        edges(graph) |> collect

        @parameters t i(t)
        @variables x(t) y(t) z(t)
        eqs = [D(x) ~ i,
               D(y) ~ z,
               D(z) ~ z]
        ode = ODESystem(eqs)
        graph = eqeq_dependencies(asgraph(ode), variable_dependencies(ode))
        edges(graph) |> collect

        has_self_loops(graph)
        is_connected(graph)
        connected_components(graph)

        # get partions whit no leaving edges
        # this means, this subset of equations is NOT used by the others
        attracting_components(graph)

        # test simplecycles
        g = SimpleDiGraph(5)
        add_edge!(g, 3=>1)
        add_edge!(g, 3=>2)
        add_edge!(g, 1=>1)
        add_edge!(g, 1=>2)
        add_edge!(g, 2=>2)
        add_edge!(g, 4=>3)
        add_edge!(g, 5=>3)
        add_edge!(g, 4=>1)
        cycles = simplecycles(g)

    end

    @testset "isalgebraic" begin
        using IOSystems: isalgebraic
        @parameters t
        @variables x(t) y
        @derivatives D'~t
        @test !isalgebraic(D(x) ~ 0)
        @test !isalgebraic(D(y) ~ 0)
        @test isalgebraic(x ~ 0)
        @test isalgebraic(y ~ 0)
        @test !isalgebraic(x ~ x^2)
        @test !isalgebraic(y ~ y^2)
        @test isalgebraic(y ~ x^2)
        @test isalgebraic(y ~ x^2)
        @test isalgebraic(x ~ y^2)
    end

    @testset "pairwise cycle free" begin
        using IOSystems: pairwise_cycle_free
        g = SimpleDiGraph(5)
        add_edge!(g, 1=>2)
        add_edge!(g, 2=>4)
        add_edge!(g, 4=>5)
        add_edge!(g, 3=>5)
        add_edge!(g, 3=>1)
        add_edge!(g, 2=>3)
        add_edge!(g, 5=>4)
        cycles = simplecycles(g)

        @test pairwise_cycle_free([1,2,3,4,5], cycles) == [3,5]
    end

    @testset "reduce algebraic states" begin
        using IOSystems: reduce_algebraic_states
        @parameters t a b i(t)
        @variables x(t) y(t) o(t) o1(t) o2(t)
        @derivatives D'~t
        # variable o can be reduced
        eqs = [D(x) ~ x + o,
               D(y) ~ y + o,
               o ~ a^b]
        reqs = reduce_algebraic_states(eqs)
        @test reqs == [D(x) ~ x + a^b,
                       D(y) ~ y + a^b]

        eqs = [D(x) ~ i + o,
               o ~ x + i]
        reqs = reduce_algebraic_states(eqs)
        @test reqs == [D(x) ~ i + (x + i)]

        eqs = [D(x) ~ i,
               o1 ~ x + o2,
               D(y) ~ i,
               o2 ~ y + o1]
        reqs = reduce_algebraic_states(eqs)
        @test reqs == [D(x) ~ i,
                       o1 ~ x + (y + o1),
                       D(y) ~ i]

        @test reduce_algebraic_states(reqs) == reqs
    end
end
