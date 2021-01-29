using Test
using IOSystems
using ModelingToolkit
using ModelingToolkit: value
using SymbolicUtils
using SymbolicUtils: operation

@info "Tests of utils.jl"

@testset "utils.jl" begin
    @testset "check macro" begin
        @variables a b c d
        IOSystems.@check Set([a, b, c]) ⊆ Set([a,b,c,d]) "Shoud be subset"
        try
            IOSystems.@check Set([a, b, c]) ⊆ Set([a,b,d]) "Shoud be subset"
        catch e
            @test e isa ArgumentError
        end
    end

    @testset "remove namespace of symbols" begin
        using IOSystems: remove_namespace
        using ModelingToolkit: renamespace, to_symbolic, rename
        @parameters t a b(t)
        a = to_symbolic(a)
        b = to_symbolic(b)
        an = rename(a, renamespace(:ns, a.name))
        bn = rename(b, renamespace(:ns, operation(b).name))
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

    @testset "eqsubstitute" begin
        using IOSystems: eqsubstitute
        @parameters t p i(t) pn inew(t)
        @variables x(t) y xn(t) yn(t)
        D = Differential(t)
        eq = D(x) ~ x + y + p + i
        @test isequal(eqsubstitute(eq, x=>xn), D(xn) ~ xn + y + p + i)
        @test isequal(eqsubstitute(eq, y=>yn), D(x) ~ x + yn + p + i)
        @test isequal(eqsubstitute(eq, p=>pn), D(x) ~ x + y + pn + i)
        @test isequal(eqsubstitute(eq, i=>inew), D(x) ~ x + y + p + inew)
    end

    @testset "test uniquenames" begin
        using IOSystems: uniquenames

        (t, a, b, mp) = value.(@parameters t a b(t) m)
        (x, y, mv) = value.(@variables x y(t) m)
        @test uniquenames([a, b, x, y, mp])
        @test !uniquenames([a, b, x, y, mp, mv])
    end

    @testset "equation type" begin
        using IOSystems
        using IOSystems: eq_type
        @parameters t
        @variables x(t) y
        D = Differential(t)

        @test eq_type(D(x) ~ 0) == (:diffeq, x.val)
        @test eq_type(D(y) ~ 0) == (:diffeq, y.val)
        @test eq_type(x ~ 0) == (:explicit_algebraic, x.val)
        @test eq_type(y ~ 0) == (:explicit_algebraic, y.val)
        @test eq_type(x ~ x^2) == (:implicit_algebraic, x.val)
        @test eq_type(y ~ y^2) == (:implicit_algebraic, y.val)
        @test eq_type(y ~ x^2) == (:explicit_algebraic, y.val)
        @test eq_type(x ~ y^2) == (:explicit_algebraic, x.val)
        @test eq_type(0 ~ x + y) == (:implicit_algebraic, nothing)
    end

    @testset "recusive subsitute" begin
        using IOSystems: recursive_substitute
        @syms a b c
        rules = Dict([a=>b+2, b=>c])
        @test isequal(recursive_substitute(a, rules), c+2)
    end
end
