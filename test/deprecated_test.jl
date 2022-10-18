@testset "deprecated" begin
    @testset "test of rename_vars" begin
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
        @variables A₊x1(t) A₊x2(t) B₊x1(t) B₊x2(t) A₊o(t) B₊o(t)
        @variables add₊ina(t) add₊inb(t)

        @test_throws ArgumentError new = rename_vars(iob, in2=:in1)
        new = rename_vars(iob, in2=:in2N, out=:outN, A₊x1=:y, a=:c)

        @parameters in2N(t) c
        @variables outN(t) y(t)
        eqs = [D(y) ~ c * in1,
               D(A₊x2) ~ in2N,
               D(B₊x1) ~ b * in3,
               D(B₊x2) ~ in4,
               outN ~ (y + A₊x2) + (B₊x1 + B₊x2)]
        @test eqs == equations(new)
        @test new.removed_eqs == [add₊ina ~ y + A₊x2,
                                  add₊inb ~ B₊x1 + B₊x2,
                                  A₊o ~ y + A₊x2,
                                  B₊o ~ B₊x1 + B₊x2]
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
end
