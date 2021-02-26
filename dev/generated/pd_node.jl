using BlockSystems
using ModelingToolkit

@parameters t
D = Differential(t)
nothing #hide

@parameters τ input(t)
@variables filtered(t)

lpf = IOBlock([D(filtered) ~ 1/τ * (- filtered + input)],
              [input], [filtered])

@parameters ω(t) v(t) τ
@variables u_i(t) u_r(t) A(t)

# explicit algebraic equation for A will be reduced at connect
voltage_source = IOBlock([A ~ 1/τ * (v/√(u_i^2 + u_r^2) - 1),
                          D(u_r) ~ -ω * u_i + A*u_r,
                          D(u_i) ~  ω * u_r + A*u_i],
                         [ω, v], [u_i, u_r])

@parameters K u_ref x_ref x(t)
@variables u(t)

droop_control = IOBlock([
    u ~ - K * (x - x_ref) + u_ref # output is the droop voltage v
    ], [x], [u])

p_filter = IOBlock(lpf, name = :p_filter)
q_filter = IOBlock(lpf, name = :q_filter)
p_droop = IOBlock(droop_control, name = :p_droop)
q_droop = IOBlock(droop_control, name = :q_droop)
v_source = IOBlock(voltage_source, name = :v_source)

nothing #hide

gfi = IOSystem([p_filter.filtered => p_droop.x,
                q_filter.filtered => q_droop.x,
                p_droop.u => v_source.ω,
                q_droop.u => v_source.v],
               [p_filter, q_filter, p_droop, q_droop, v_source],
               name = :GridForming,
               namespace_map = [p_filter.input => :P_in,
                                q_filter.input => :Q_in,
                                p_filter.filtered => :p_filtered,
                                q_filter.filtered => :q_filtered,
                                # parameter names which match VSIVoltagePT1
                                v_source.τ => :τ_v, # time constant voltage control delay
                                p_filter.τ => :τ_P, # time constant active power measurement
                                q_filter.τ => :τ_Q, # time constant reactive power measurement
                                p_droop.K  => :K_P, # droop constant frequency droop
                                q_droop.K  => :K_Q, # droop constant voltage droop
                                q_droop.u_ref => :V_r, # reference/ desired voltage
                                p_droop.u_ref => :ω_r, # reference/ desired frequency
                                p_droop.x_ref => :P, # active (real) power infeed
                                q_droop.x_ref => :Q], # reactive (imag) power infeed                .
               outputs = [v_source.u_i, v_source.u_r])

@parameters u_i(t) u_r(t) i_i(t) i_r(t)
@variables P_in(t) Q_in(t)
pow = IOBlock([P_in ~ u_r*i_r + u_i*i_i,
               Q_in ~ u_i*i_r - u_r*i_i],
              [u_i, u_r, i_i, i_r], [P_in, Q_in], name=:pow)

gfi2 = IOSystem([gfi.u_i => pow.u_i,
                 gfi.u_r => pow.u_r,
                 pow.P_in => gfi.P_in,
                 pow.Q_in => gfi.Q_in],
                [gfi, pow], outputs=[gfi.u_i, gfi.u_r], name=:gfi_i_to_u)
nothing #hide

connected = connect_system(gfi2)
nothing #hide

using PowerDynamics
using NetworkDynamics: ODEVertex
struct IONode{T} <: AbstractNode
    block::IOBlock
    generated::T
    parameters::Vector{Float64}
end

function IONode(blk::IOBlock, parameters::Dict)
    # BlockSpec: blk must by of type (i_r, i_i) ↦ (u_r, u_i)
    spec = BlockSpec([:i_r, :i_i], [:u_r, :u_i])
    @assert spec(blk) "Block has to follow PowerDynamics i/o conventions!"
    # TODO check parameters
    gen = generate_io_function(blk,
                               f_states=[blk.u_r, blk.u_i],
                               f_inputs=[blk.i_r, blk.i_i],
                               f_params=keys(parameters), warn=false);
    IONode(blk, gen, collect(values(parameters)))
end

import PowerDynamics.construct_vertex
function construct_vertex(ion::IONode)
    gen = ion.generated
    function rhs!(dx, x, e_s, e_d, _p, t)
        i = total_current(e_s, e_d)
        gen.f_ip(dx, x, (real(i), imag(i)), ion.parameters, t)
    end
    ODEVertex(f! = rhs!, dim = length(gen.states), mass_matrix = gen.massm, sym = Symbol.(gen.states))
end

import PowerDynamics.symbolsof
PowerDynamics.symbolsof(ionode::IONode) = Symbol.(ionode.generated.states)
PowerDynamics.dimension(ionode::IONode) = length(ionode.generated.states)

@info "connected system:" connected.iparams

para = Dict(
    :τ_v => 0.005,    # time constant voltage control delay
    :τ_P => 0.5,      # time constant active power measurement
    :τ_Q => 0.5,      # time constant reactive power measurement
    :K_P => 0.396,     # droop constant frequency droop
    :K_Q => 0.198,     # droop constant voltage droop
    :V_r => 1.0,   # reference/ desired voltage
    :P   => 0.303, # active (real) power infeed
    :Q   => 0.126, # reactive (imag) power infeed                .
    :ω_r => 0.0)   # refrence/ desired frequency
nothing #hide

using PowerDynamics: SlackAlgebraic, PiModelLine,PowerGrid,find_operationpoint
using OrderedCollections: OrderedDict

buses=OrderedDict(
    "bus1"=> IONode(connected, para),
    "bus2"=> SlackAlgebraic(U=1))

branches=OrderedDict(
    "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=4.999131600798035-1im*15.263086523179553, y_shunt_km=0.0528/2, y_shunt_mk=0.0528/2))

powergrid = PowerGrid(buses,branches)
operationpoint = find_operationpoint(powergrid)
timespan= (0.0,0.1)
nothing #hide

fault1 = ChangeInitialConditions(node="bus1", var=:u_r, f=Inc(0.2))
solution1 = simulate(fault1, powergrid, operationpoint, timespan)
using Plots
plot(solution1.dqsol)

using Test
@testset "Compare IONode vs. PD Node" begin
    for i in 1:100
        para = Dict(
            :τ_v => rand(), # time constant voltage control delay
            :τ_P => rand(), # time constant active power measurement
            :τ_Q => rand(), # time constant reactive power measurement
            :K_P => rand(), # droop constant frequency droop
            :K_Q => rand(), # droop constant voltage droop
            :V_r => rand(), # reference/ desired voltage
            :P   => rand(), # active (real) power infeed
            :Q   => rand(), # reactive (imag) power infeed                .
            :ω_r => 0.0)    # refrence/ desired frequency

        # create IONode
        node_bs = IONode(connected, para)
        f_bs = construct_vertex(node_bs).f!

        # create PDNode
        # the PD node does not have the explicit ω_r parameter
        para_pd = delete!(copy(para), :ω_r)
        node_pd = VSIVoltagePT1(; NamedTuple(para_pd)...)
        f_pd = construct_vertex(node_pd).f!

        # create fake "edge data", 4 incoming, 4 outgooing with 4 states each
        es = [randn(4) for i in 1:4]
        ed = [randn(4) for i in 1:4]

        # select random time
        t = rand()

        # chose random initial state and account for initial ω in PD node
        x_bs = randn(4)
        x_pd = copy(x_bs)
        x_pd[3] = - para[:K_P] * (x_bs[3] - para[:P])

        # create arrays for the results
        dx_bs = similar(x_bs)
        dx_pd = similar(x_pd)

        # call both functions
        f_bs(dx_bs, x_bs, es, ed, nothing, t)
        f_pd(dx_pd, x_pd, es, ed, nothing, t)

        # compare results
        # we have to correct variable 3 of bs implementation to match dω
        @test dx_bs[1] ≈ dx_pd[1]                # u_r
        @test dx_bs[2] ≈ dx_pd[2]                # u_i
        @test - para[:K_P] * dx_bs[3] ≈ dx_pd[3] # ω
        @test dx_bs[4] ≈ dx_pd[4]                # q_filtered
    end
end
nothing #hide

para = Dict(:τ_v=>rand(),:τ_P=>rand(), :τ_Q=>rand(),
            :K_P=>rand(), :K_Q=>rand(), :V_r=>rand(),
            :P=>rand(), :Q=>rand(), :ω_r=>0.0)

node_bs = IONode(connected, para)
f_bs = construct_vertex(node_bs).f!
para_pd = delete!(copy(para), :ω_r)
node_pd = VSIVoltagePT1(; NamedTuple(para_pd)...)
f_pd = construct_vertex(node_pd).f!

es = [randn(4) for i in 1:4]
ed = [randn(4) for i in 1:4]
t = rand()

# chose random initial state and account for initial ω in PD node
x_bs = randn(4)
x_pd = copy(x_bs)
x_pd[3] = - para[:K_P] * (x_bs[3] - para[:P])
dx = similar(x_bs)

using BenchmarkTools
@btime $f_bs($dx, $x_bs, $es, $ed, nothing, $t)
@btime $f_pd($dx, $x_pd, $es, $ed, nothing, $t)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

