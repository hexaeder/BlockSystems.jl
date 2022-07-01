using Test
using BlockSystems
using ModelingToolkit
using ModelingToolkit: get_iv, get_eqs, get_states, value
using Graphs

@info "Tests of BlockSystems.jl"

@testset "BlockSystems.jl" begin
    @testset "namespaced accessors" begin
        @parameters t i1(t) i2(t) a
        @variables x1(t) x2(t) o1(t) o2(t)
        D = Differential(t)
        eqs  = [D(x1) ~ a*i1, o1~i1, D(x2) ~ i2, o2~i2]
        iob = IOBlock(eqs, [i1, i2], [o1, o2], name=:ns)
        @test Set(BlockSystems.namespace_inputs(iob)) == Set([iob.i1, iob.i2])
        @test Set(BlockSystems.namespace_outputs(iob)) == Set([iob.o1, iob.o2])
        @test Set(BlockSystems.namespace_istates(iob)) == Set([iob.x1, iob.x2])
        @test Set(BlockSystems.namespace_iparams(iob)) == Set([iob.a])
    end

    @testset "getproperty of symbol, num and symbolic" begin
        @parameters t i1(t) i2(t) a
        @variables x1(t) x2(t) o1(t) o2(t)
        D = Differential(t)
        eqs  = [D(x1) ~ a*i1, o1~i1, D(x2) ~ i2, o2~i2]
        iob = IOBlock(eqs, [i1, i2], [o1, o2], name=:ns)

        ga = getproperty(iob, :i1)
        gb = getproperty(iob, i1)
        gc = getproperty(iob, iob.i1)
        @test length(Set([ga,gb,gc])) == 1

        ga = getproperty(iob, :i2)
        gb = getproperty(iob, i2)
        gc = getproperty(iob, iob.i2)
        @test length(Set([ga,gb,gc])) == 1

        ga = getproperty(iob, :a)
        gb = getproperty(iob, a)
        gc = getproperty(iob, iob.a)
        @test length(Set([ga,gb,gc])) == 1

        ga = getproperty(iob, :x1)
        gb = getproperty(iob, x1)
        gc = getproperty(iob, iob.x1)
        @test length(Set([ga,gb,gc])) == 1

        ga = getproperty(iob, :x2)
        gb = getproperty(iob, x2)
        gc = getproperty(iob, iob.x2)
        @test length(Set([ga,gb,gc])) == 1

        ga = getproperty(iob, :o1)
        gb = getproperty(iob, o1)
        gc = getproperty(iob, iob.o1)
        @test length(Set([ga,gb,gc])) == 1

        ga = getproperty(iob, :o2)
        gb = getproperty(iob, o2)
        gc = getproperty(iob, iob.o2)
        @test length(Set([ga,gb,gc])) == 1
    end

    @testset "creation of IOBlocks" begin
        @parameters t i1(t) i2(t) a b
        @variables x1(t) x2(t) o1(t) o2(t)
        D = Differential(t)
        eqs = [D(x1) ~ a*i1,
               D(x2) ~ b*i2,
               o1 ~ a*x1,
               o2 ~ b*x2]

        iob = IOBlock(eqs, [i1, i2], [o1, o2], name=:iob)
        @test Set(iob.inputs) == Set([i1, i2])
        @test Set(iob.iparams) == Set([a, b])
        @test Set(iob.istates) == Set([x1, x2])
        @test Set(iob.outputs) == Set([o1, o2])
        @test Set(iob.removed_states) == Set()

        # x1 is a state not an input. Assign it as input will elad to warn + name conflict between inputs and states
        @test_throws ArgumentError IOBlock(eqs, [x1], [o1,o2], warn=false)

        @test_throws ArgumentError IOBlock(eqs, [i1,i2], [i1,o1,o2], warn=false)

        @parameters i a(t)
        @variables x(t) o(t)
        sys = ODESystem( [D(x) ~ a * i], name=:foo)
        aeq = [i ~ 2*a + i1]
        # name conflict betwen ios and odes
        @test_throws ArgumentError IOBlock(:name, [i.val], [a.val], [], [x.val], sys, aeq; warn=true)

        @parameters t a(t) p
        @variables x(t) y(t)
        D = Differential(t)
        IOBlock([x ~ 2 + a], [], [x])

        # implecit output, warn
        @test_logs (:warn, ) IOBlock([D(y) ~ 2 + a], [], [x])

        @test_throws ArgumentError IOBlock([a ~ 2 + x], [], [x])
        @test_throws ArgumentError IOBlock([p ~ 2 + x], [p], [x], iv=t)
        @test_throws ArgumentError IOBlock([a ~ 2 + x], [a], [x], iv=t)
        @test_logs (:warn, ) IOBlock([a ~ 2 + x + 2*a], [a], [x], iv=t)
    end

    @testset "test of create_namespace_promotions" begin
        using BlockSystems: create_namespace_promotions
        @parameters t
        Aa, Ba, Ab, Bc = value.(@parameters A₊a B₊a(t) A₊b B₊c)
        Ax, Bx, Ay, Bz = value.(@variables A₊x(t) B₊x(t) A₊y(t) B₊z(t))
        b, c = value.(@parameters b c)
        y, z = value.(@variables y(t) z(t))

        prom = create_namespace_promotions([Aa, Ba, Ab, Bc, Ax, Bx, Ay, Bz], [], true)
        @test Set(values(prom)) == Set([Aa, Ba, b, c, Ax, Bx, y, z])

        prom = create_namespace_promotions([Aa, Ba, Ab, Bc, Ax, Bx, Ay, Bz], [b, y], true)
        @test Set(values(prom)) == Set([Aa, Ba, Ab, c, Ax, Bx, Ay, z])

        prom = create_namespace_promotions([Aa, Ba, Ab, Bc, Ax, Bx, Ay, Bz], [], false)
        @test Set(values(prom)) == Set([Aa, Ba, Ab, Bc, Ax, Bx, Ay, Bz])

        prom = create_namespace_promotions([Aa, Ba, Ab, Bc, Ax, Bx, Ay, Bz], [b, y], false)
        @test Set(values(prom)) == Set([Aa, Ba, Ab, Bc, Ax, Bx, Ay, Bz])

        ABa, Ba, a = value.(@parameters A₊B₊a B₊a(t) a)
        prom = create_namespace_promotions([ABa, Ba, a], [], true)
        @test Set(values(prom)) == Set([ABa, Ba, a])
 end

    @testset "IOBlock from other IOBlock" begin
        @parameters t i1(t) i2(t) a b
        @variables x1(t) x2(t) o1(t) o2(t)
        D = Differential(t)
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
        @test get_eqs(iob1.system) == get_eqs(iob2.system) == get_eqs(iob3.system)
        @test iob1.system.name == iob1.name == :iob1
        @test iob2.system.name == iob2.name == :iob2
        @test iob3.system.name == iob3.name
    end

    @testset "test creation of namespace map" begin
        @parameters t i(t) a b
        @variables x(t) o(t)
        D = Differential(t)
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

    @testset "iosystem asserts" begin
        @parameters t i(t)
        @variables o(t)
        iob1 = IOBlock([o ~ i],[i],[o],name=:name)
        iob2 = IOBlock([o ~ i],[i],[o],name=:name)

        # namespace collision
        @test_throws ArgumentError IOSystem([iob1.o => iob2.i], [iob1, iob2])

        iob2 = IOBlock([o ~ i],[i],[o])
        # multiple connections to same input
        @test_throws ArgumentError IOSystem([iob1.o => iob1.i, iob2.o => iob1.i], [iob1, iob2])
        # make sure that alle of form input => output
        @test_throws ArgumentError IOSystem([iob1.o => iob1.o], [iob1, iob2])
        @test_throws ArgumentError IOSystem([iob1.i => iob1.i], [iob1, iob2])

        # assert that input maps refer to open inputs
        @test_throws ArgumentError IOSystem([iob1.o => iob2.i], [iob1, iob2],
                                            namespace_map = [iob2.i => i])
        # assert that rhs of input map is unique
        iob3 = IOBlock([o ~ i],[i],[o])
        @test_throws ArgumentError IOSystem([iob1.o => iob2.i],
                                             [iob1, iob2, iob3],
                                             namespace_map = [iob1.i => i, iob3.i => i])

        # test assertions for iparams and istats map
        @parameters t a i(t) b c
        @variables x(t) o(t) y(t)
        D = Differential(t)
        iob1 = IOBlock([D(x)~ i, o~a*x], [i], [o], name=:iob1)
        iob2 = IOBlock([D(x)~ i, o~a*x], [i], [o], name=:iob2)
        IOSystem([iob1.o => iob2.i], [iob1, iob2])
        # rhs unique
        @test_throws ArgumentError IOSystem([iob1.o => iob2.i], [iob1, iob2],
                                            namespace_map = [iob1.a=>b, iob2.a=>b])
        @test_throws ArgumentError IOSystem([iob1.o => iob2.i], [iob1, iob2],
                                            namespace_map = [iob1.x=>y, iob2.x=>y])

        # rhs unique
        @test_throws ArgumentError IOSystem([iob1.o => iob2.i], [iob1, iob2],
                                            namespace_map = [iob1.o=>y, iob2.o=>y])
    end

    function test_complete_namespace_promotions(ios)
        eqs = vcat([ModelingToolkit.namespace_equations(iob.system) for iob in ios.systems]...)
        allvars = [(get_variables(eq.lhs) ∪ get_variables(eq.rhs)) for eq in eqs]
        allvars = union(allvars...) |> unique
        allvars = setdiff(allvars, [BlockSystems.get_iv(ios)])
        allkeys = keys(ios.namespace_map)
        # closed inputs should not appear in allkeys
        @test isempty(Set(allkeys) ∩ Set(last.(ios.connections)))
        # but there should be no namespace collision with them either..
        allkeys = Set(allkeys) ∪ Set(last.(ios.connections))
        @test allunique(allkeys)
        @test Set(allkeys) == Set(allvars)
    end

    @testset "test creation of systems" begin
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
        @parameters t i1(t) i2(t) a b ina(t) inb(t)
        @variables x1(t) x2(t) o(t) add(t)
        D = Differential(t)
        eqs1  = [D(x1) ~ a*i1, D(x2)~i2, o~x1+x2]
        iob1 = IOBlock(eqs1, [i1, i2], [o], name=:iob1)

        eqs2  = [D(x1) ~ b*i1, D(x2)~i2, o~x1+x2]
        iob2 = IOBlock(eqs2, [i1, i2], [o], name=:iob2)

        ioadd = IOBlock([add ~ ina + inb], [ina, inb], [add], name=:add)

        # try with auto namespacing
        sys = IOSystem([iob1.o => ioadd.ina, iob2.o => ioadd.inb],
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
        sys1 = IOSystem([iob1.o => ioadd.ina, iob2.o => ioadd.inb],
                       [iob1, iob2, ioadd],
                       namespace_map = Dict(iob1.i1 => in1,
                                            iob1.i2 => in2,
                                            iob2.i1 => in3,
                                            iob2.i2 => in4,
                                            iob1.a => p1,
                                            iob2.b => p2,
                                            iob1.x1 => y1,
                                            iob1.x2 => y2,
                                            iob2.x1 => x1,
                                            iob2.x2 => x2,
                                            ioadd.add => out),
                       outputs = [ioadd.add],
                       name=:sys)
        sys2 = IOSystem([iob1.o => ioadd.ina, iob2.o => ioadd.inb],
                       [iob1, iob2, ioadd],
                       namespace_map = Dict(iob1.i1 => :in1,
                                            iob1.i2 => :in2,
                                            iob2.i1 => :in3,
                                            iob2.i2 => :in4,
                                            iob1.a => :p1,
                                            iob2.b => :p2,
                                            iob1.x1 => :y1,
                                            iob1.x2 => :y2,
                                            iob2.x1 => :x1,
                                            iob2.x2 => :x2,
                                            ioadd.add => :out),
                       outputs = [:out],
                       name=:sys)
        @test Set(sys1.inputs) == Set([in1, in2, in3, in4])
        @test Set(sys1.iparams) == Set([p1, p2])
        @test Set(sys1.istates) == Set([y1, y2, x1, x2, iob1.o, iob2.o])
        @test Set(sys1.outputs) == Set([out])
        test_complete_namespace_promotions(sys1)

        @test Set(sys2.inputs) == Set([in1, in2, in3, in4])
        @test Set(sys2.iparams) == Set([p1, p2])
        @test Set(sys2.istates) == Set([y1, y2, x1, x2, iob1.o, iob2.o])
        @test Set(sys2.outputs) == Set([out])
        test_complete_namespace_promotions(sys2)

        # provide partial maps
        sys = IOSystem([iob1.o => ioadd.ina, iob2.o => ioadd.inb],
                       [iob1, iob2, ioadd],
                       namespace_map = Dict(iob1.i1 => in1,
                                            iob1.i2 => in2,
                                            iob1.a => p1,
                                            iob1.x1 => y1,
                                            iob1.x2 => y2,
                                            ioadd.add => out),
                       outputs = [out], # provide outputs as namespaced variable
                       name=:sys)
        @test Set(sys.inputs) == Set([in1, in2, i1, i2])
        @test Set(sys.iparams) == Set([p1, b])
        @test Set(sys.istates) == Set([y1, y2, x1, x2, iob1.o, iob2.o])
        @test Set(sys.outputs) == Set([out])
        test_complete_namespace_promotions(sys)

        # check for argument error for bad outputs argument
        sys1 = IOSystem([iob1.o => ioadd.ina, iob2.o => ioadd.inb],
                       [iob1, iob2, ioadd],
                       namespace_map = Dict(ioadd.add => out),
                       outputs = [ioadd.add], # provide outputs as namespaced variable
                       name=:sys)
        sys2 = IOSystem([iob1.o => ioadd.ina, iob2.o => ioadd.inb],
                       [iob1, iob2, ioadd],
                       namespace_map = Dict(ioadd.add => out),
                       outputs = [out], # provide outputs as namespaced variable
                       name=:sys)
        @test Set(sys1.istates) == Set(sys2.istates)
        @test Set(sys1.outputs) == Set(sys2.outputs)
        @test_throws ArgumentError IOSystem([iob1.o => ioadd.ina, iob2.o => ioadd.inb],
                                            [iob1, iob2, ioadd],
                                            namespace_map = Dict(ioadd.add => out),
                                            outputs = [p1], # provide outputs as namespaced variable
                                            name=:sys)
    end

    @testset "algebraic connections" begin
        @parameters t p1 p2 i(t)
        @variables o1(t) o2(t) o3(t)
        dt = Differential(t)

        @named b1 = IOBlock([o1 ~ p1], [], [o1])
        @named b2 = IOBlock([o2 ~ p2], [], [o2])
        @named b3 = IOBlock([dt(o3) ~ i], [i], [o3])

        sys = IOSystem([b1.o1 + b2.o2 => b3.i],
                       [b1, b2, b3], outputs=[b3.o3])
        con = connect_system(sys, verbose=false)
        @test equations(con) == [dt(o3) ~ p1 + p2]
    end

    @testset "test BlockSpec" begin
        # test the constructors
        @parameters t
        @parameters uᵢ(t) uᵣ(t)
        @variables x(t) iᵢ(t) iᵣ(t)
        bs1 = BlockSpec([:uᵣ, :uᵢ], [:iᵣ, :iᵢ])
        bs2 = BlockSpec(value.([uᵣ, uᵢ]), value.([iᵣ, iᵢ]))
        bs3 = BlockSpec([uᵣ, uᵢ], [iᵣ, iᵢ])
        @test bs1.inputs == bs2.inputs == bs3.inputs
        @test bs1.outputs == bs2.outputs == bs3.outputs

        D = Differential(t)
        eqs  = [D(x) ~ x, iᵢ ~ uᵢ, iᵣ ~ uᵣ]
        iob1 = IOBlock(eqs, [uᵢ, uᵣ], [iᵢ, iᵣ])
        iob2 = IOBlock(eqs, [uᵢ], [iᵢ, iᵣ])
        iob3 = IOBlock(eqs, [uᵢ, uᵣ], [iᵢ])

        @test fulfills(iob1, bs1) == bs1(iob1) == true
        @test fulfills(iob2, bs1) == bs1(iob2) == false
        @test fulfills(iob3, bs1) == bs1(iob3) == false

        sys1 = IOSystem([], [iob1, iob2])
        sys2 = IOSystem([], [iob1, iob2], outputs=[iᵢ, iᵣ],
                        namespace_map = [iob1.iᵢ => iᵢ,
                                         iob1.iᵣ => iᵣ,
                                         iob1.uᵢ => uᵢ,
                                         iob1.uᵣ => uᵣ])

        @test BlockSpec([:uᵣ, :uᵢ], [:iᵣ, :iᵢ])(sys1) == false # no name promotion
        @test BlockSpec([:uᵣ, :uᵢ], [:iᵣ, :iᵢ])(sys2) == false # to much inputs

        # test the strict stuff
        @parameters t
        @parameters i1(t) i2(t)
        @variables o1(t) o2(t)
        iob = IOBlock([o1 ~ i1, o2 ~ i2], [i1, i2], [o1, o2], iv=t)

        @test BlockSpec([:i1, :i2], [:o1])(iob) == true
        @test BlockSpec([:i1, :i2], [:o1], in_strict=1, out_strict=1)(iob) == false
        @test BlockSpec([:i1], [:o1], in_strict=1, out_strict=1)(iob) == false
        @test BlockSpec([:i1], [:o1], in_strict=0, out_strict=0)(iob) == true

        @test BlockSpec([:i1], [:o1, :o2])(iob) == false
        @test BlockSpec([:i1], [:o1, :o2], in_strict=0, out_strict=0)(iob) == true
        @test BlockSpec([:i1], [:o1, :o2], in_strict=0, out_strict=1)(iob) == true
    end

    @testset "test auto connections" begin
        using BlockSystems: autoconnections

        @variables t a(t) b(t) c(t)
        @named blk1 = IOBlock([a ~ 1, b ~ 2, c ~ 3],
                              [], [a,b,c])
        @named blk2 = IOBlock([a ~ 1],
                              [], [a])

        @variables o1(t)
        @parameters a(t) b(t) c(t)
        @named blk3 = IOBlock([o1 ~ a + b + c],
                              [a, b, c],
                              [o1])
        @variables o2(t)
        @named blk4 = IOBlock([o2 ~ a + b],
                              [a, b],
                              [o2])

        @test Set(autoconnections([blk1, blk3])) ==
            Set([blk1.a => blk3.a,
                 blk1.b => blk3.b,
                 blk1.c => blk3.c])

        @test Set(autoconnections([blk1, blk2, blk3])) ==
            Set([blk1.b => blk3.b,
                 blk1.c => blk3.c])

        @test Set(autoconnections([blk1, blk3, blk4])) ==
            Set([blk1.a => blk3.a,
                 blk1.a => blk4.a,
                 blk1.b => blk3.b,
                 blk1.b => blk4.b,
                 blk1.c => blk3.c])

        IOSystem(:autocon, [blk1, blk3])
        IOSystem(:autocon, [blk1, blk2, blk3])
        IOSystem(:autocon, [blk1, blk3, blk4])
    end
end
