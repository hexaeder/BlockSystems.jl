using Test
using BlockSystems
using ModelingToolkit
using LightGraphs

# ordering and simplification should not matter for equity of equations!
import Base.==
function ==(A::Vector{Equation}, B::Vector{Equation})
    Set(simplify.(A)) == Set(simplify.(B))
end

@info "Testes of transformations.jl"

@testset "transformations.jl" begin
    @testset "pairwise cycle free" begin
        using BlockSystems: pairwise_cycle_free
        g = SimpleDiGraph(5)
        add_edge!(g, 1=>2)
        add_edge!(g, 2=>4)
        add_edge!(g, 4=>5)
        add_edge!(g, 3=>5)
        add_edge!(g, 3=>1)
        add_edge!(g, 2=>3)
        add_edge!(g, 5=>4)
        cycles = simplecycles(g)
        @test pairwise_cycle_free(g) == [3,5]
    end

    @testset "remove superfluous states" begin
        using BlockSystems: remove_superfluous_states
        @parameters t a b i(t)
        @variables x(t) y(t) o(t) o1(t) o2(t)
        D = Differential(t)
        eqs = [D(x) ~ x,
               D(y) ~ y,
               o ~ a^b]
        reqs = remove_superfluous_states(eqs, t, [x])
        @test reqs == eqs[[1]]
        reqs = remove_superfluous_states(eqs, t, [y])
        @test reqs == eqs[[2]]
        reqs = remove_superfluous_states(eqs, t, [o])
        @test reqs == eqs[[3]]
        eqs = [D(x) ~ y,
               D(y) ~ x+o,
               o ~ a^b]
        reqs = remove_superfluous_states(eqs, t, [y])
        @test reqs == eqs
        reqs = remove_superfluous_states(eqs, t, [x])
        @test reqs == eqs
        reqs = remove_superfluous_states(eqs, t, [o])
        @test reqs == eqs[[3]]
        eqs = [D(x) ~ x+o,
               D(y) ~ y,
               o ~ a^b]
        reqs = remove_superfluous_states(eqs, t, [x])
        @test reqs == eqs[[1,3]]
        reqs = remove_superfluous_states(eqs, t, [y])
        @test reqs == eqs[[2]]
        reqs = remove_superfluous_states(eqs, t, [o])
        @test reqs == eqs[[3]]
    end

    @testset "remove algebraic states" begin
        using BlockSystems: remove_algebraic_states
        @parameters t a b i(t)
        @variables x(t) y(t) o(t) o1(t) o2(t)
        D = Differential(t)
        # variable o can be removed
        eqs = [D(x) ~ x + o,
               D(y) ~ y + o,
               o ~ a^b]
        reqs = remove_algebraic_states(eqs)
        @test reqs[1] == [D(x) ~ x + a^b,
                          D(y) ~ y + a^b]
        @test reqs[2] == [o ~ a^b]

        # variable o1 and o2 can be removed
        eqs = [D(x) ~ x + o1+o2,
               D(y) ~ y + o2,
               o1 ~ a^b + o2,
               o2 ~ a-b]
        reqs = remove_algebraic_states(eqs)
        @test reqs[1] == [D(x) ~ x + a^b + a-b + a-b,
                          D(y) ~ y + a-b]
        @test reqs[2] == [o1 ~ a^b + a - b,
                          o2 ~ a - b]

        eqs = [D(x) ~ i + o,
               o ~ x + i]
        reqs = remove_algebraic_states(eqs)
        @test reqs[1] == [D(x) ~ i + (x + i)]
        @test reqs[2] == [o ~ x + i]

        eqs = [D(x) ~ i,
               o1 ~ x + o2,
               D(y) ~ i,
               o2 ~ y + o1]
        reqs = remove_algebraic_states(eqs)
        @test reqs[1] == [D(x) ~ i,
                          D(y) ~ i,
                          o1 ~ x + (y + o1)]
        @test reqs[2] == [o2 ~ y + o1]

        rreqs = remove_algebraic_states(reqs[1])
        @test rreqs[1] == reqs[1]
        @test rreqs[2] == []

        # test skip condition
        eqs = [D(x) ~ i + o,
               o ~ x + i]
        reqs = remove_algebraic_states(eqs, skip=[o])
        @test reqs[1] == eqs
        @test reqs[2] == []

        eqs = [D(x) ~ i,
               o1 ~ x + o2,
               D(y) ~ i,
               o2 ~ y + o1]
        reqs = remove_algebraic_states(eqs, skip=[o1])
        @test reqs[1] == [D(x) ~ i,
                          D(y) ~ i,
                          o1 ~ x + (y + o1)]
        @test reqs[2] == [o2 ~ y + o1]

        reqs = remove_algebraic_states(eqs, skip=[o2])
        @test reqs[1] == [D(x) ~ i,
                          D(y) ~ i,
                          o2 ~ y + (x + o2)]
        @test reqs[2] == [o1 ~ x + o2]
    end

    @testset "test of connect_system" begin
        #=
                    +------------+
        in1 --> i1 -|    iob1    |
        in2 --> i2 -|(x1, x2)(a) |-o--+     +-----+
                    +------------+    +-ina-| add |- add ---> out
                    +------------+    +-inb-|     |
        in3 --> i1 -|    iob2    |-o--+     +-----+
        in4 --> i2 -|(x1, x2)(b) |
                    +------------+
        =#
        # same system as in testset before
        @parameters t i1(t) i2(t) a b ina(t) inb(t)
        @variables x1(t) x2(t) o(t) add(t)
        D = Differential(t)
        eqs1  = [D(x1) ~ a*i1, D(x2)~i2, o~x1+x2]
        iob1 = IOBlock(eqs1, [i1, i2], [o], name=:A)
        eqs2  = [D(x1) ~ b*i1, D(x2)~i2, o~x1+x2]
        iob2 = IOBlock(eqs2, [i1, i2], [o], name=:B)
        ioadd = IOBlock([add ~ ina + inb], [ina, inb], [add], name=:add)
        @parameters in1(t) in2(t) in3(t) in4(t)
        @variables out(t)
        sys = IOSystem([iob1.o => ioadd.ina, iob2.o => ioadd.inb],
                       [iob1, iob2, ioadd],
                       namespace_map = Dict(iob1.i1 => in1,
                                            iob1.i2 => in2,
                                            iob2.i1 => in3,
                                            iob2.i2 => in4,
                                            ioadd.add => out),
                       outputs = [out],
                       name=:sys)
        iob = connect_system(sys)

        @test Set(sys.inputs) == Set(iob.inputs)
        @test Set(sys.iparams) == Set(iob.iparams)
        # the connection will delete the states iob1.o and iob2.o
        @test Set(sys.istates) == Set(iob.istates) ∪ Set([iob1.o, iob2.o])
        @test Set(sys.outputs) == Set(iob.outputs)
        @test iob.name == sys.name

        @variables A₊x1(t) A₊x2(t) B₊x1(t) B₊x2(t) A₊o(t) B₊o(t)
        eqs = [D(A₊x1) ~ a * in1,
               D(A₊x2) ~ in2,
               D(B₊x1) ~ b * in3,
               D(B₊x2) ~ in4,
               out ~ (A₊x1 + A₊x2) + (B₊x1 + B₊x2)]
        @test eqs == iob.system.eqs
        @test iob.removed_eqs == [A₊o ~ A₊x1 + A₊x2, B₊o ~ B₊x1 + B₊x2]

        @testset "test rename_vars" begin
            @test_throws ArgumentError new = rename_vars(iob, in2=:in1)
            new = rename_vars(iob, in2=:in2N, out=:outN, A₊x1=:y, a=:c)

            @parameters in2N(t) c
            @variables outN(t) y(t)
            eqs = [D(y) ~ c * in1,
                   D(A₊x2) ~ in2N,
                   D(B₊x1) ~ b * in3,
                   D(B₊x2) ~ in4,
                   outN ~ (y + A₊x2) + (B₊x1 + B₊x2)]
            @test eqs == new.system.eqs
            @test new.removed_eqs == [A₊o ~ y + A₊x2, B₊o ~ B₊x1 + B₊x2]
        end
    end

    @testset "system with removed equations in subsystem" begin
        #=
           +------+         +---------+
        +->| x'=i |--+------| o=a1+a2 |--o
        |  +------+  |   +--|         |
        |  +------+  |   |  +---------+
        +--| o = i|<-+   |
           +------+      |
           +------+      |
        +->| x'=i |--+---+
        |  +------+  |
        |  +------+  |
        +--| o = i|<-+
           +------+
        =#
        @parameters t i(t)
        @variables x(t) o(t)
        D = Differential(t)

        b1 = IOBlock([D(x) ~ i], [i], [x], name=:a)
        b2 = IOBlock([o ~ i], [i], [o], name=:b)
        subsys = IOSystem([b2.o => b1.i, b1.x => b2.i], [b1, b2], name=:A,
                          namespace_map = [b1.x=>x],
                          outputs = [x])
        @test Set(subsys.removed_states) == Set()

        subsysA = connect_system(subsys)
        @test subsysA.system.eqs == [D(x) ~ x]
        @test Set(subsysA.removed_states) == Set([o])
        @test subsysA.removed_eqs == [o ~ x]

        subsysB = IOBlock(subsysA, name=:B)
        @test subsysB.name == :B
        @test subsysB.system.name == :B
        @test subsysB.system.eqs == [D(x) ~ x]
        @test Set(subsysB.removed_states) == Set([o])
        @test subsysB.removed_eqs == [o ~ x]

        @parameters a1(t) a2(t)
        adder = IOBlock([o ~ a1 + a2], [a1, a2], [o], name=:add)

        system = IOSystem([subsysA.x => adder.a1, subsysB.x => adder.a2],
                          [adder, subsysA, subsysB])

        @test Set(system.removed_states) == Set([subsysA.o, subsysB.o])

        systemblock = connect_system(system)
        @test systemblock.system.name == systemblock.name == system.name
        @test Set(systemblock.removed_states) == Set(system.removed_states)

        using BlockSystems: namespace_rem_eqs
        @test systemblock.removed_eqs == vcat(namespace_rem_eqs(subsysA), namespace_rem_eqs(subsysB))

        # psst, i'm putting a test for function_generation here, don't tell ma
        gen = generate_io_function(systemblock,
                                   f_states=[systemblock.A₊x, systemblock.B₊x],
                                   f_rem_states=[systemblock.A₊o, systemblock.B₊o])
        @test Set(gen.states) == Set(Sym{Real}.([:B₊x, :A₊x, :add₊o]))
        @test Set(gen.inputs) == Set()
        @test Set(gen.params) == Set()
        @test Set(gen.rem_states) == Set(Sym{Real}.([:B₊o, :A₊o]))

        out = zeros(2)
        st = rand(3)
        gen.g_ip(out, st, (), (), (), 0.0)
        @test out == st[1:2]
        @test gen.g_oop(st, (), (), (), 0.0) == st[1:2]
    end
end
