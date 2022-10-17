#!julia --startup-file=no

using Pkg
Pkg.activate(temp=true)
Pkg.add("TimerOutputs")

@info "Add BlockSystems at rev $(ARGS[1])"
Pkg.add(url=expanduser("~/.julia/dev/BlockSystems/"), rev=ARGS[1])
Pkg.precompile()

using TimerOutputs

println("using time")
@timeit "using time" begin
    using BlockSystems
end

BlockSystems.WARN[] = false

println("create blocks")
@timeit "create blocks" begin
    @parameters t i1(t) i2(t) a b ina(t) inb(t)
    @variables x1(t) x2(t) o(t) add(t)
    D = Differential(t)
    eqs1  = [D(x1) ~ a*i1, D(x2)~i2, o~x1+x2]
    iob1 = IOBlock(eqs1, [i1, i2], [o], name=:iob1)

    eqs2  = [D(x1) ~ b*i1, D(x2)~i2, o~x1+x2]
    iob2 = IOBlock(eqs2, [i1, i2], [o], name=:iob2)

    ioadd = IOBlock([add ~ ina + inb], [ina, inb], [add], name=:add)
end

println("create systems")
@timeit "system creation" begin
    @timeit "create sys" begin
        sys = IOSystem([iob1.o => ioadd.ina, iob2.o => ioadd.inb],
                       [iob1, iob2, ioadd],
                       name=:sys)
    end

    @timeit "create sys1" begin
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
    end
    @timeit "create sys2" begin
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
    end
end

println("connect systems")
@timeit "connect systems" begin
    con = connect_system(sys)
    con1 = connect_system(sys1)
    con2 = connect_system(sys2)
end
println("generate functions")
@timeit "generate functions" begin
    generate_io_function(con);
    generate_io_function(con1);
    generate_io_function(con2);
end

BlockSystems.WARN[] = true

print_timer()
