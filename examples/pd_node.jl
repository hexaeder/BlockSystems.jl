#=
## Generate `PowerDynamics.jl`-Node with `BlockSystems`
We start by defining some general stuff...
=#
using BlockSystems
using ModelingToolkit

@parameters t
D = Differential(t)

#=
### Defining the common Blocks
We need some common parts such as filters or integrators. The actual Blocks will
be defined by using thes 'blueprints'. Stuff like this should go to a library of
components at some point.
=#
## low pass filter
@parameters τ input(t)
@variables filtered(t)

lpf = IOBlock([D(filtered) ~ 1/τ * (- filtered + input)],
              [input], [filtered])

## integrator
@parameters x(t)
@variables int(t)

integrator = IOBlock([D(int) ~ x], [x], [int])

## voltage source

@parameters ω(t) v(t) τ
@variables u_i(t) u_r(t)

voltage_source = IOBlock([D(u_r) ~ -ω * u_i + (u_r^2 + u_i^2 - v^2)*τ*u_r,
                          D(u_i) ~  ω * u_r + (u_r^2 + u_i^2 - v^2)*τ*u_i],
                         [ω, v], [u_i, u_r])

## Droop control

@parameters K u_ref x_ref x(t)
@variables u(t)

droop_control = IOBlock([
    u ~ - K * (x - x_ref) + u_ref # output is the droop voltage v
    ], [x], [u])

#=
After defining the blueprints we can create the Blocks of the system
based on them:
=#

p_filter = IOBlock(lpf, name = :q_filter)
q_filter = IOBlock(lpf, name = :r_filter)
p_droop = IOBlock(droop_control, name = :p_droop)
q_droop = IOBlock(droop_control, name = :q_droop)
v_source = IOBlock(voltage_source, name = :v_source)

nothing #hide
#=
We can connect the blocks in in order to get a system

```


```
=#

sys = IOSystem([voltage_source.ω => p_droop.u,
                voltage_source.v => q_droop.u,
                p_droop.x => p_filter.filtered,
                q_droop.x => q_filter.filtered],
               [p_filter, q_filter, p_droop, q_droop, voltage_source],
               name = :GridForming,
               namespace_map = [p_filter.input => :p, q_filter.input => :q, q_droop.u => :v],
               outputs = [voltage_source.u_i, voltage_source.u_r])

##
@parameters u_i(t) u_r(t) i_i(t) i_r(t)
@variables P(t) Q(t)
power = IOBlock([P ~ u_r*i_r + u_i*i_i,
                 Q ~ u_i*i_r - u_r*i_i],
                [u_i, u_r, i_i, i_r], [P, Q], name=:power)

gfi2 = IOSystem([power.u_i => gfi.u_i,
                power.u_r => gfi.u_r,
                gfi.p => power.P,
                gfi.q => power.Q],
               [gfi, power], outputs=[gfi.u_i, gfi.u_r])

##

connected_gfi2 = connect_system(gfi2)

##

connected_gfi2.system.eqs
gen = generate_io_function(connected_gfi2,
                           f_states=[connected_gfi2.u_r,connected_gfi2.u_i],
                           f_inputs=[connected_gfi2.i_r,connected_gfi2.i_i]);

##

gen.states
Symbol.(gen.states)

using PowerDynamics
using NetworkDynamics: ODEVertex
struct IONode{T} <: AbstractNode
    block::IOBlock
    generated::T
    parameters::Vector{Float64}
end

function IONode(blk::IOBlock, parameters::Dict)
    spec = BlockSpec([:i_r, :i_i], [:u_r, :u_i])
    @assert spec(blk) "Block has to follow PowerDynamics i/o conventions!"
    # TODO check parameters
    gen = generate_io_function(blk,
                               f_states=[blk.u_r, blk.u_i],
                               f_inputs=[blk.i_r, blk.i_i],
                               f_params=keys(parameters));
    IONode(blk, gen, collect(values(parameters)))
end

import PowerDynamics.construct_vertex
function construct_vertex(ion::IONode)
    gen = ion.generated
    function rhs!(dx, x, e_s, e_d, _p, t)
        i = total_current(e_s, e_d)
        gen.f_ip(dx, x, [real(i), imag(i)], ion.parameters, t)
    end
    ODEVertex(f! = rhs!, dim = length(gen.states), mass_matrix = gen.massm, sym = Symbol.(gen.states))
end

import PowerDynamics.symbolsof
PowerDynamics.symbolsof(ionode::IONode) = Symbol.(ionode.generated.states)
PowerDynamics.dimension(ionode::IONode) = length(ionode.generated.states)


para = Dict(p_filter.τ => 1.0,
            q_filter.τ => 1.0,
            voltage_source.τ => 1.0,
            p_droop.u_ref => 1.0,
            q_droop.u_ref => 1.0,
            p_droop.x_ref => 1.0,
            q_droop.x_ref => 1.0,
            q_droop.K => 1.0,
            p_droop.K => 1.0)
n = IONode(connected_gfi2, para)
n.generated.f_ip

connected_gfi2.iparams

using PowerDynamics: SlackAlgebraic, PiModelLine,PowerGrid,find_operationpoint
using OrderedCollections: OrderedDict

##
# τ_P, τ_Q, K_P, K_Q = rand(1:10, 4)
# P,Q,V_r = rand(1.:10.,3)

buses=OrderedDict(
    "bus1"=> IONode(connected_gfi2, para),
    "bus2"=> SlackAlgebraic(U=1))

branches=OrderedDict(
    "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=4.999131600798035-1im*15.263086523179553, y_shunt_km=0.0528/2, y_shunt_mk=0.0528/2))

##
powergrid = PowerGrid(buses,branches)
operationpoint = find_operationpoint(powergrid)
timespan= (0.0,5.)

operationpoint |> typeof |> fieldnames
operationpoint.vec

# simulating a frequency perturbation at node 1
fault1 = ChangeInitialConditions(node="bus1", var=:u_r, f=Inc(0.2))
solution1 = simulate(fault1, powergrid, operationpoint, timespan)
using Plots
plot(solution1.dqsol)
solutions1.dqsol
