using Test
using IOSystems
using ModelingToolkit
using LinearAlgebra

@info "Testes of function_generation.jl"

@testset "function_generation.jl" begin
    @testset "transform algebraic equations" begin
        using IOSystems: transform_algebraic_equations
        @parameters t a b
        @variables x(t) y(t)
        @derivatives D'~t
        eqs = [D(x)~y,
               y ~ x,
               2*b + 8 ~ a + a,
               1 ~ x + y]
        eqs = transform_algebraic_equations(eqs)
        @test eqs == [D(x)~y,
                      0 ~ x - y,
                      0 ~ (a + a) - (2*b + 8),
                      0 ~ (x + y) - 1]
    end

    @testset "reorder of equations" begin
        using IOSystems: reorder_by_states
        @parameters t
        @variables x(t) y(t) a(t) b(t)
        @derivatives D'~t
        eqs = [D(x)~x, D(y)~y, 0 ~ x+y+a, 0 ~ x+y+b]
        @test reorder_by_states(eqs, [x,y,a,b]) == eqs[[1,2,3,4]]
        @test reorder_by_states(eqs, [y,x,a,b]) == eqs[[2,1,3,4]]
        @test reorder_by_states(eqs, [a,x,y,b]) == eqs[[3,1,2,4]]
        @test reorder_by_states(eqs, [x,a,y,b]) == eqs[[1,3,2,4]]
        @test reorder_by_states(eqs, [y,a,x,b]) == eqs[[2,3,1,4]]
        @test reorder_by_states(eqs, [a,y,x,b]) == eqs[[3,2,1,4]]
        @test reorder_by_states(eqs, [a,y,b,x]) == eqs[[3,2,4,1]]
        @test reorder_by_states(eqs, [y,a,b,x]) == eqs[[2,3,4,1]]
        @test reorder_by_states(eqs, [b,a,y,x]) == eqs[[3,4,2,1]]
        @test reorder_by_states(eqs, [a,b,y,x]) == eqs[[3,4,2,1]]
        @test reorder_by_states(eqs, [y,b,a,x]) == eqs[[2,3,4,1]]
        @test reorder_by_states(eqs, [b,y,a,x]) == eqs[[3,2,4,1]]
        @test reorder_by_states(eqs, [b,x,a,y]) == eqs[[3,1,4,2]]
        @test reorder_by_states(eqs, [x,b,a,y]) == eqs[[1,3,4,2]]
        @test reorder_by_states(eqs, [a,b,x,y]) == eqs[[3,4,1,2]]
        @test reorder_by_states(eqs, [b,a,x,y]) == eqs[[3,4,1,2]]
        @test reorder_by_states(eqs, [x,a,b,y]) == eqs[[1,3,4,2]]
        @test reorder_by_states(eqs, [a,x,b,y]) == eqs[[3,1,4,2]]
        @test reorder_by_states(eqs, [y,x,b,a]) == eqs[[2,1,3,4]]
        @test reorder_by_states(eqs, [x,y,b,a]) == eqs[[1,2,3,4]]
        @test reorder_by_states(eqs, [b,y,x,a]) == eqs[[3,2,1,4]]
        @test reorder_by_states(eqs, [y,b,x,a]) == eqs[[2,3,1,4]]
        @test reorder_by_states(eqs, [x,b,y,a]) == eqs[[1,3,2,4]]
        @test reorder_by_states(eqs, [b,x,y,a]) == eqs[[3,1,2,4]]

        eqs = [0 ~ x+y]
        @test reorder_by_states(eqs, [x]) == eqs
    end

    @testset "generate_massmatrix" begin
        using IOSystems: generate_massmatrix
        @parameters t
        @variables x(t) y(t) a(t) b(t)
        @derivatives D'~t
        eqs = [D(x)~x, D(y)~y, 0 ~ x+y+a, 0 ~ x+y+b]
        @test generate_massmatrix(eqs) == Diagonal([1,1,0,0])
        @test generate_massmatrix(eqs[[1,3,2,4]]) == Diagonal([1,0,1,0])
        eqs = [D(x)~x, D(y)~y]
        @test generate_massmatrix(eqs) == I
    end

    @testset "all_static" begin
        using IOSystems: all_static
        @parameters t a b i1(t) i2(t)
        @variables o1(t) o2(t)
        @derivatives D'~t
        eqs = [o1 ~ a*i1 + i2,
               o2 ~ b*i1 + i2]
        @test all_static(eqs)
        eqs = [D(o1) ~ a*i1 + i2,
               o2 ~ b*i1 + i2]
        @test !all_static(eqs)
        eqs = [o1 ~ a*i1 + i2 + o2,
               o2 ~ b*i1 + i2]
        @test !all_static(eqs)
    end

    @testset "ode function generation" begin
        @parameters t a i(t)
        @variables x(t) y(t)
        @derivatives D'~t
        iob = IOBlock([D(x)~a*i, D(y)~x], [i], [x, y])

        iof = generate_io_function(iob);
        @test iof.massm == I
        @test isequal(iof.inputs, [Sym{Real}(:i)])
        @test isequal(iof.states, Sym{Real}.([:x, :y]))
        @test isequal(iof.params, [Sym{Real}(:a)])
        @test isequal(iof.params, [Sym{Real}(:a)])
        @test iof.g_ip === nothing
        @test iof.g_oop === nothing
        @test isempty(iof.rem_states)

        @test iof.f_oop([-1,1], [4], [-1], 0) == [-4, -1]
        a = Vector{Float64}(undef, 2)
        iof.f_ip(a, [-1,1], [4], [-1], 0)
        @test a == [-4, -1]

        iof = generate_io_function(iob, f_states=[y])
        @test isequal(iof.states, Sym{Real}.([:y, :x]))
        iof = generate_io_function(iob, f_states=[iob.y])
        @test isequal(iof.states, Sym{Real}.([:y, :x]))
    end

    @testset "static function generation" begin
        @parameters t a i1(t) i2(t)
        @variables x(t) y(t)
        @derivatives D'~t
        iob = IOBlock([x~a*i1, y~i2], [i1, i2], [x, y])

        iof = generate_io_function(iob);

        @test isequal(iof.inputs, Sym{Real}.([:i1, :i2]))
        @test isequal(iof.states, Sym{Real}.([:x, :y]))
        @test isequal(iof.params, [Sym{Real}(:a)])

        @test iof.f_oop([4, 1], [-1], 0) == [-4, 1]
        a = Vector{Float64}(undef, 2)
        iof.f_ip(a, [-1,1], [-1], 0)
        @test a == [1, 1]
    end

    @testset "parameter ordering" begin
        using ModelingToolkit: makesym

        @parameters t a b i1(t) i2(t)
        @variables o1(t) o2(t)
        @derivatives D'~t

        iob = IOBlock([D(o1)~a*i1 + b, D(o2)~o1+i2], [i1, i2], [o1, o2])

        allequal(v1, v2) = all([isequal(e1, e2) for (e1, e2) ∈ zip(v1, v2)])

        gen = generate_io_function(iob, f_states=[o1, o2])
        @test allequal(gen.states, makesym.([o1, o2], states=[]))
        gen = generate_io_function(iob, f_states=[o2, o1])
        @test allequal(gen.states, makesym.([o2, o1], states=[]))

        gen = generate_io_function(iob, f_inputs=[i1, i2])
        @test allequal(gen.inputs, makesym.([i1, i2], states=[]))
        gen = generate_io_function(iob, f_inputs=[i2, i1])
        @test allequal(gen.inputs, makesym.([i2, i1], states=[]))

        gen = generate_io_function(iob, f_params=[a, b])
        @test allequal(gen.params, makesym.([a, b], states=[]))
        gen = generate_io_function(iob, f_params=[b, a])
        @test allequal(gen.params, makesym.([b, a], states=[]))
    end
end
