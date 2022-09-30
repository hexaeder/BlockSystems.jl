using Test
using BlockSystems
using ModelingToolkit
using ModelingToolkit: get_iv, get_eqs, get_states
using Graphs

# ordering and simplification should not matter for equity of equations!
import Base.==
function ==(A::Vector{Equation}, B::Vector{Equation})
    Set(simplify.(A)) == Set(simplify.(B))
end

@info "Testes of transformations.jl"

@testset "transformations.jl" begin
    @testset "cyclebreaking vertices" begin
        using BlockSystems: cyclebreaking_vertices
        g = SimpleDiGraph(5)
        add_edge!(g, 1=>2)
        add_edge!(g, 2=>4)
        add_edge!(g, 4=>5)
        add_edge!(g, 3=>5)
        add_edge!(g, 3=>1)
        add_edge!(g, 2=>3)
        add_edge!(g, 5=>4)
        cycles = simplecycles(g)
        @test cyclebreaking_vertices(g) == [[2,3,1], [5,4]]
    end

    @testset "remove superfluous states" begin
        using BlockSystems: remove_superfluous_states
        @parameters t a b i(t)
        @variables x(t) y(t) o(t) o1(t) o2(t)
        D = Differential(t)
        eqs = [D(x) ~ x,
               D(y) ~ y,
               o ~ a^b]
        reqs = remove_superfluous_states(IOBlock(eqs, [], [x]; iv=t)) |> equations
        @test reqs == eqs[[1]]
        reqs = remove_superfluous_states(IOBlock(eqs, [], [y]; iv=t)) |> equations
        @test reqs == eqs[[2]]
        reqs = remove_superfluous_states(IOBlock(eqs, [], [o]; iv=t)) |> equations
        @test reqs == eqs[[3]]
        eqs = [D(x) ~ y,
               D(y) ~ x+o,
               o ~ a^b]
        reqs = remove_superfluous_states(IOBlock(eqs, [], [y]; iv=t)) |> equations
        @test reqs == eqs
        reqs = remove_superfluous_states(IOBlock(eqs, [], [x]; iv=t)) |> equations
        @test reqs == eqs
        reqs = remove_superfluous_states(IOBlock(eqs, [], [o]; iv=t)) |> equations
        @test reqs == eqs[[3]]
        eqs = [D(x) ~ x+o,
               D(y) ~ y,
               o ~ a^b]
        reqs = remove_superfluous_states(IOBlock(eqs, [], [x]; iv=t)) |> equations
        @test reqs == eqs[[1,3]]
        reqs = remove_superfluous_states(IOBlock(eqs, [], [y]; iv=t)) |> equations
        @test reqs == eqs[[2]]
        reqs = remove_superfluous_states(IOBlock(eqs, [], [o]; iv=t)) |> equations
        @test reqs == eqs[[3]]
    end

    @testset "remove algebraic states" begin
        @parameters t a b i(t)
        @variables x(t) y(t) o(t) o1(t) o2(t)
        D = Differential(t)
        # variable o can be removed
        eqs = [D(x) ~ x + o,
               D(y) ~ y + o,
               o ~ a^b]
        reqs = substitute_algebraic_states(IOBlock(eqs,[],[x,y]))

        @test equations(reqs) == [D(x) ~ x + a^b,
                          D(y) ~ y + a^b]
        @test reqs.removed_eqs == [o ~ a^b]

        # variable o1 and o2 can be removed
        eqs = [D(x) ~ x + o1+o2,
               D(y) ~ y + o2,
               o1 ~ a^b + o2,
               o2 ~ a-b]
        reqs = substitute_algebraic_states(IOBlock(eqs,[],[x,y]))
        @test equations(reqs) == [D(x) ~ x + a^b + a-b + a-b,
                          D(y) ~ y + a-b]
        @test reqs.removed_eqs == [o1 ~ a^b + a - b,
                          o2 ~ a - b]

        eqs = [D(x) ~ i + o,
               o ~ x + i]
        reqs = substitute_algebraic_states(IOBlock(eqs,[],[x]))
        @test equations(reqs) == [D(x) ~ i + (x + i)]
        @test reqs.removed_eqs == [o ~ x + i]

        eqs = [D(x) ~ i,
               o1 ~ x + o2,
               D(y) ~ i,
               o2 ~ y + 2*o1]
        reqs = substitute_algebraic_states(IOBlock(eqs,[],[x,y]))
        equations(reqs)
        @test equations(reqs) == [D(x) ~ i,
                          D(y) ~ i,]
        @test reqs.removed_eqs == [o2 ~ y + 2*(-x + -y),
                                   o1 ~ -x + -y]

        # TODO: multiple applications should not remove more equations!
        rreqs = substitute_algebraic_states(reqs)
        @test equations(rreqs) == equations(reqs)
        @test rreqs.removed_eqs == reqs.removed_eqs

        # test skip condition
        eqs = [D(x) ~ i + o,
               o ~ x + i]
        reqs = substitute_algebraic_states(IOBlock(eqs,[],[x,o]))
        @test equations(reqs) == eqs
        @test reqs.removed_eqs == []

        eqs = [D(x) ~ i,
               o1 ~ x + o2,
               D(y) ~ i,
               o2 ~ y + 2*o1]
        reqs = substitute_algebraic_states(IOBlock(eqs,[],[x,y,o1]))
        @test equations(reqs) == [D(x) ~ i,
                          D(y) ~ i,
                          o1 ~ -x + -y]
        @test reqs.removed_eqs == [o2 ~ y + 2*o1]

        reqs = substitute_algebraic_states(IOBlock(eqs,[],[x,y,o2]))
        @test equations(reqs) == [D(x) ~ i,
                          D(y) ~ i,
                          o2 ~ -y + -2x]
        @test reqs.removed_eqs == [o1 ~ x + o2]
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
        @variables add₊ina(t) add₊inb(t)
        eqs = [D(A₊x1) ~ a * in1,
               D(A₊x2) ~ in2,
               D(B₊x1) ~ b * in3,
               D(B₊x2) ~ in4,
               out ~ (A₊x1 + A₊x2) + (B₊x1 + B₊x2)]
        @test eqs == get_eqs(iob.system)
        @test iob.removed_eqs == [add₊ina ~ A₊x1 + A₊x2,
                                  add₊inb ~ B₊x1 + B₊x2,
                                  A₊o ~ A₊x1 + A₊x2,
                                  B₊o ~ B₊x1 + B₊x2]

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
            @test eqs == get_eqs(new.system)
            @test new.removed_eqs == [add₊ina ~ y + A₊x2,
                                      add₊inb ~ B₊x1 + B₊x2,
                                      A₊o ~ y + A₊x2,
                                      B₊o ~ B₊x1 + B₊x2]
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
        @test get_eqs(subsysA.system) == [D(x) ~ x]
        @test Set(subsysA.removed_states) ⊇ Set([o])
        @test subsysA.removed_eqs == [b1.i ~ x,
                                      b2.i ~ x,
                                      o ~ x]

        subsysB = IOBlock(subsysA, name=:B)
        @test subsysB.name == :B
        @test subsysB.system.name == :B
        @test get_eqs(subsysB.system) == [D(x) ~ x]
        @test Set(subsysB.removed_states) == Set([b1.i, b2.i, o])
        subsysB.removed_states
        @test subsysB.removed_eqs == [b1.i~x, b2.i~x, o ~ x]

        @parameters a1(t) a2(t)
        adder = IOBlock([o ~ a1 + a2], [a1, a2], [o], name=:add)

        system = IOSystem([subsysA.x => adder.a1, subsysB.x => adder.a2],
                          [adder, subsysA, subsysB])

        @test Set(system.removed_states) ⊇ Set([subsysA.o, subsysB.o])

        systemblock = connect_system(system)
        @test systemblock.system.name == systemblock.name == system.name
        @test Set(systemblock.removed_states) ⊇ Set(system.removed_states)

        using BlockSystems: namespace_rem_eqs
        @test systemblock.removed_eqs ⊇ vcat(namespace_rem_eqs(subsysA), namespace_rem_eqs(subsysB))

        # psst, i'm putting a test for function_generation here, don't tell ma
        gen = generate_io_function(systemblock,
                                   f_states=[systemblock.A₊x, systemblock.B₊x],
                                   f_rem_states=[systemblock.A₊o, systemblock.B₊o],
                                   warn=false)
        @test Set(gen.states) == Set(Sym{Real}.([:B₊x, :A₊x, :add₊o]))
        @test Set(gen.inputs) == Set()
        @test Set(gen.params) == Set()
        @test Set(gen.rem_states) ⊇ Set(Sym{Real}.([:B₊o, :A₊o]))

        out = zeros(8)
        st = rand(3)
        gen.g_ip(out, st, (), (), (), 0.0)
        @test out[1:2] == st[1:2]
        @test gen.g_oop(st, (), (), (), 0.0)[1:2] == st[1:2]
    end

    @testset "substitution of differential states" begin
        # first block
        @parameters t i(t)
        @variables o(t)
        D = Differential(t)
        A = IOBlock([D(o) ~ 1 + D(i)], [i], [o])
        @test length(BlockSystems.rhs_differentials(A)) == 1

        # second block, does not contain derivative of output
        @variables x(t)
        @parameters a(t)
        B = IOBlock([x ~ a + 1], [], [x])
        @test length(BlockSystems.rhs_differentials(B)) == 0

        # third block, contains derivative of output
        @variables y(t)
        @parameters b(t)
        C = IOBlock([D(y) ~ b + 1], [], [y])
        @test length(BlockSystems.rhs_differentials(B)) == 0

        sysAB = IOSystem([B.x => A.i], [A, B]) |> connect_system
        @test length(BlockSystems.rhs_differentials(sysAB)) == 1

        sysAC = IOSystem([C.y => A.i], [A, C]) |> connect_system
        @test length(BlockSystems.rhs_differentials(sysAC)) == 0

        @test equations(sysAC.system) == [D(o) ~ 2 + b, D(y) ~ 1 + b]

        # fourth block with more complex stuff
        @variables x(t) y(t)
        B = IOBlock([D(x) ~ a, y ~ x + 3], [], [y])

        sysAB = IOSystem([B.y => A.i], [A, B]) |> connect_system
        @test length(BlockSystems.rhs_differentials(sysAB)) == 0
    end

    @testset "subs of differential states" begin
        @parameters t a(t) b
        @variables x(t) y(t) z(t)
        D = Differential(t)

        blk = IOBlock([D(x) ~ 1 + D(y),
                       D(y) ~ a + b], [], [])
        blk2 = substitute_derivatives(blk, verbose=true)
        @test equations(blk2) == [D(x) ~ 1 + a + b,
                                  D(y) ~ a + b]

        blk = IOBlock([D(x) ~ 1 + D(y),
                       y ~ 2 + b], [], [])
        blk2 = substitute_derivatives(blk, verbose=true)
        @test equations(blk2) == [D(x) ~ 1 ,
                                  y ~ 2 + b]

        blk = IOBlock([D(x) ~ 1 + D(y),
                       y ~ 2 + a], [], [])
        blk2 = substitute_derivatives(blk, verbose=true)
        @test equations(blk2) == [D(x) ~ 1 + D(a),
                                  y ~ 2 + a]

        blk = IOBlock([D(x) ~ 1 + D(y),
                       y ~ 2 + z,
                       D(z) ~ 5*x + b], [], [])
        blk2 = substitute_derivatives(blk, verbose=true)
        @test equations(blk2) == [D(x) ~ 1 + 5*x + b,
                                  y ~ 2 + z,
                                  D(z) ~ 5*x + b]

        @parameters i(t)
        @variables o(t)
        @named b1 = IOBlock([o ~ sin(t)], [], [o])
        @named b2 = IOBlock([o ~ cos(t)], [], [o])
        @named b3 = IOBlock([o ~ D(i)], [i], [o])
        c = IOSystem([b1.o => b3.i],
                     [b1, b2, b3], outputs=[b3.o]) |> connect_system
        @test equations(c) == [b3.o ~ cos(t)]
        c = IOSystem([b1.o + b2.o => b3.i],
                     [b1, b2, b3], outputs=[b3.o]) |> connect_system
        @test equations(c) == [b3.o ~ cos(t) - sin(t)]
        c = IOSystem([b1.o * b2.o => b3.i],
                     [b1, b2, b3], outputs=[b3.o]) |> connect_system
        @test equations(c) == [b3.o ~ (cos(t))^2 - (sin(t))^2]
    end

    @testset "Don't remove implicit differential equations" begin
        using BlockSystems:rhs_differentials

        @parameters t y(t)
        D = Differential(t)
        @variables x(t)

        blk1 = IOBlock([x ~ y + D(y)], [y], [x]) # implicit differential equation

        @parameters x(t)
        @variables y(t)
        blk2 = IOBlock([D(y) ~ x], [x], [y])

        sys = IOSystem([blk1.x => blk2.x, blk2.y => blk1.y], [blk1, blk2], outputs = [blk1.x, blk2.y])
        sys = connect_system(sys)

        @test isequal(rhs_differentials(sys), Set{SymbolicUtils.Symbolic}()) # all rhs differentials have been resolved
        @test isequal(equations(sys), [D(y) ~ x, 0 ~ y])
    end
   
    @testset "set p" begin
        @parameters t a(t) b
        @variables x(t) y(t) z(t)
        D = Differential(t)

        blk = IOBlock([D(x) ~ 1 + D(y) + a*b], [], [])

        blk2 = set_p(blk, :a=>0)
        @test equations(blk2) == [D(x) ~ 1 + D(y)]
        blk2 = set_p(blk, blk.a=>0)
        @test equations(blk2) == [D(x) ~ 1 + D(y)]
        blk2 = set_p(blk, a=>0)
        @test equations(blk2) == [D(x) ~ 1 + D(y)]

        blk2 = set_p(blk, :b=>1)
        @test equations(blk2) == [D(x) ~ 1 + D(y) + a]
        blk2 = set_p(blk, blk.b=>1)
        @test equations(blk2) == [D(x) ~ 1 + D(y) + a]
        blk2 = set_p(blk, b=>1)
        @test equations(blk2) == [D(x) ~ 1 + D(y) + a]

        blk2 = set_p(blk, Dict(blk.a=>2, b=>4))
        @test equations(blk2) == [D(x) ~ 9 + D(y)]

        blk2 = set_p(blk, blk.a=>2, b=>4; warn=false)
        @test equations(blk2) == [D(x) ~ 9 + D(y)]

        @test_throws ArgumentError blk2 = set_p(blk, :a=>:bla)
        @test_throws ArgumentError blk2 = set_p(blk, x=>2.0)
    end

    @testset "set p for inputs" begin
        @parameters t a(t) b(t)
        @variables x(t) y(t) z(t)
        D = Differential(t)

        blk = IOBlock([D(x) ~ 1 + D(y) + a + b], [a, b], [])

        blk2 = set_p(blk, :a=>0)
        @test equations(blk2) == [D(x) ~ 1 + D(y) + b]
        blk2 = set_p(blk, blk.a=>0)
        @test equations(blk2) == [D(x) ~ 1 + D(y) + b]
        blk2 = set_p(blk, a=>0)
        @test equations(blk2) == [D(x) ~ 1 + D(y) + b]

        blk2 = set_p(blk, :b=>1)
        @test equations(blk2) == [D(x) ~ 1 + D(y) + a + 1]
        blk2 = set_p(blk, blk.b=>1)
        @test equations(blk2) == [D(x) ~ 1 + D(y) + a + 1]
        blk2 = set_p(blk, b=>1)
        @test equations(blk2) == [D(x) ~ 1 + D(y) + a + 1]

        blk2 = set_p(blk, Dict(blk.a=>2, b=>4))
        @test equations(blk2) == [D(x) ~ 7 + D(y)]

        blk2 = set_p(blk, blk.a=>2, b=>4; warn=false)
        @test equations(blk2) == [D(x) ~ 7 + D(y)]

        @test_throws ArgumentError blk2 = set_p(blk, :a=>:bla)
        @test_throws ArgumentError blk2 = set_p(blk, x=>2.0)
    end

    @testset "simplify eqs" begin
        @parameters t a
        @variables x(t)
        D = Differential(t)
        blk = IOBlock([D(x) ~ x/a], [], [])

        # this errors on MTK@0.8 currently
        ModelingToolkit.namespace_equations(simplify_eqs(blk).system)

        @variables a(t) b(t) c(t) d(t)
        term = substitute(a - b, b=>c + d) |> Symbolics.unwrap
        terms = simplify(term)
        @test_broken !any(Symbolics.metadata.(Symbolics.get_variables(terms)).==nothing)
    end

    @testset "set input" begin
        @variables t o(t)
        @parameters i(t)
        blk = IOBlock([o ~ 2*i], [i], [o])

        blk2 = set_input(blk, :i=>4)
        @test equations(blk2) == [o ~ 8]
        blk2 = set_input(blk, i=>4)
        @test equations(blk2) == [o ~ 8]
        blk2 = set_input(blk, blk.i=>4)
        @test equations(blk2) == [o ~ 8]

        @test blk.name == blk2.name
        @test blk2.removed_eqs ⊇ [i~4]

        @parameters a(t)
        blk2 = set_input(blk, :i=>a)
        @test equations(blk2) == [o ~ 2a]
        blk2 = set_input(blk, i=>a)
        @test equations(blk2) == [o ~ 2a]
        blk2 = set_input(blk, blk.i=>a)
        @test equations(blk2) == [o ~ 2a]
    end

    @testset "direct connect" begin
        function safe()
            @variables t o(t) o1(t) o2(t)
            @parameters i(t) i1(t) i2(t)
            @named blkA = IOBlock([o1 ~ 1 + i, o2 ~ 1 - i], [i], [o1, o2])
            @named blkB = IOBlock([o ~ 2 + i1 + 2*i2], [i1, i2], [o])
            return blkA, blkB
        end
        blkA, blkB = safe()

        # first we try the macro as long as those symbols are not in the namespace
        sysA = @connect blkA.o1 => blkB.i1
        sysB = @connect blkA.(o1, o2) => blkB.(i1, i2)
        sysC = @connect blkA.(o1,o2) => blkB.(i1 ,i2) autopromote=false
        sysD = @connect blkA.(o1,o2) => blkB.(i1 ,i2) autopromote=false outputs=[blkB.o]
        sysE = @connect blkA.(o1,o2) => blkB.(i1 ,i2) outputs=:remaining
        sysF = @connect blkA.(o1,o2) => blkB.(i1 ,i2) outputs=:remaining name=:foo

        # create symbols in namespace
        @variables t o(t) o1(t) o2(t)
        @parameters i(t) i1(t) i2(t)

        @test Set(sysA.outputs) == Set([o1, o2, o])
        @test Set(sysA.inputs) == Set([i, i2])

        @test Set(sysB.outputs) == Set([o1, o2, o])
        @test Set(sysB.inputs) == Set([i])

        @test Set(sysC.outputs) == Set([blkA.o1, blkA.o2, blkB.o])

        @test Set(sysD.outputs) == Set([blkB.o])

        @test Set(sysF.outputs) == Set([o])
        @test Set(sysE.outputs) == Set([o])

        @test sysF.name == :foo


        @named blkC = IOBlock([o ~ π*i], [i], [o])
        sys = @connect blkA.(o1,o2)=>blkB.(i1,i2) blkB.o=>blkC.i outputs=:remaining name=:blk autopromote=:false
        @test sys.name == :blk
    end
end
