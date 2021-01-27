using Test
using IOSystems
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
        @test pairwise_cycle_free(g) == [3,5]
    end

    @testset "remove superfluous states" begin
        using IOSystems: remove_superfluous_states
        @parameters t a b i(t)
        @variables x(t) y(t) o(t) o1(t) o2(t)
        @derivatives D'~t
        eqs = [D(x) ~ x,
               D(y) ~ y,
               o ~ a^b]
        reqs = remove_superfluous_states(eqs, [x])
        @test reqs == eqs[[1]]
        reqs = remove_superfluous_states(eqs, [y])
        @test reqs == eqs[[2]]
        reqs = remove_superfluous_states(eqs, [o])
        @test reqs == eqs[[3]]
        eqs = [D(x) ~ y,
               D(y) ~ x+o,
               o ~ a^b]
        reqs = remove_superfluous_states(eqs, [y])
        @test reqs == eqs
        reqs = remove_superfluous_states(eqs, [x])
        @test reqs == eqs
        reqs = remove_superfluous_states(eqs, [o])
        @test reqs == eqs[[3]]
        eqs = [D(x) ~ x+o,
               D(y) ~ y,
               o ~ a^b]
        reqs = remove_superfluous_states(eqs, [x])
        @test reqs == eqs[[1,3]]
        reqs = remove_superfluous_states(eqs, [y])
        @test reqs == eqs[[2]]
        reqs = remove_superfluous_states(eqs, [o])
        @test reqs == eqs[[3]]
    end

    @testset "remove algebraic states" begin
        using IOSystems: remove_algebraic_states
        @parameters t a b i(t)
        @variables x(t) y(t) o(t) o1(t) o2(t)
        @derivatives D'~t
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
        @variables out(t)
        sys = IOSystem([ioadd.ina => iob1.o, ioadd.inb => iob2.o],
                       [iob1, iob2, ioadd],
                       inputs_map = Dict(iob1.i1 => in1,
                                         iob1.i2 => in2,
                                         iob2.i1 => in3,
                                         iob2.i2 => in4),
                       outputs_map = Dict(ioadd.add => out),
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
    end
end
