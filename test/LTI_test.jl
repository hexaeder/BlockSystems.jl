using Test
using BlockSystems

@info "Tests of LTI.jl"

@testset "LTI.jl" begin
    @testset "_state_matrix tests" begin
        using BlockSystems: _state_matrix
        @variables t x(t) y(t) z(t)
        @parameters a b c d
        @test all(isequal.(_state_matrix([a*x + b*y, c*x + d*y], [x,y]), [a b; c d]))
        @test_throws ErrorException _state_matrix([a*x*y + b*y, c*x + d*y], [x,y])
    end

    @testset "identify lti" begin
        @variables t x1(t) x2(t) o1(t) o2(t)
        @parameters i1(t) i2(t)
        dt = Differential(t)
        blk = IOBlock([dt(x1) ~ x1 - x2 + 2*i1 + i2,
                       dt(x2) ~ -x1 + x2 + i2,
                       o1 ~ x1 + x2,
                       o2 ~ x2 + i1],
                      [i1, i2], [o1, o2, x1])

        A, B, C, D = identify_lti(blk)

        # TODO: test all the error parts

        @test A == [1 -1; -1 1]
        @test B == [2 1; 0 1]
        @test C == [1 1; 0 1; 1 0]
        @test D == [0 0; 1 0; 0 0]
    end
end
