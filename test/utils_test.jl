using Test
using IOSystems
using ModelingToolkit

@info "Tests of utils.jl"

@testset "utils.jl" begin
    @testset "test uniquenames" begin
        using IOSystems: uniquenames

        (t, a, b, mp) = value.(@parameters t a b(t) m)
        (x, y, mv) = value.(@variables x y(t) m)
        @test uniquenames([a, b, x, y, mp])
        @test !uniquenames([a, b, x, y, mp, mv])
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

    @testset "eqsubstitute" begin
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

    @testset "check macro" begin
        @variables a b c d
        IOSystems.@check Set([a, b, c]) ⊆ Set([a,b,c,d]) "Shoud be subset"
        try
            IOSystems.@check Set([a, b, c]) ⊆ Set([a,b,d]) "Shoud be subset"
        catch e
            @test e isa ArgumentError
        end
    end
end
